import os
import sys
import gzip
import string
import shutil
import numpy as np
import prody as pr

from datetime import datetime 
from multiprocessing import Pool, Manager
from dataclasses import dataclass
import logomaker
import matplotlib.pyplot as plt
import pandas as pd

from scipy.stats import mode
from scipy.spatial.distance import cdist
from itertools import product, permutations

import qbits

from smallprot import pdbutils, query, cluster_loops, smallprot_config
from smallprot import logger, peputils, constant, plot, extract_master, struct_analysis

@dataclass
class Struct_info:
    trunc_info: str
    loop_info: str
    loop_len: int 
    cent_pdb:str
    clust_num: int   
    min_rmsd: float
    median_rmsd: float
    redundancy:str = "Unknown" 
    clust_num_2nd: int = 0

class Loop_sse:
    """A small protein in the process of generation by MASTER and Qbits.

    Parameters
    ----------
    para : Parameter

    Attributes
    ----------
    pdbs : list
        List of paths to PDB files for the protein at each step of 
        the design process.
    exclusion_pdbs : list
        List of paths to PDB files containing structures defining a region 
        of excluded volume for the designed protein at each step of the 
        design process.
    qbit_reps : list
        List of paths to PDB files for the qbit reps to be added to 
        the designed protein.

    Methods
    -------
    build_protein()
        Build the protein according to the specifications provided to the 
        constructor function.
    loop_seed_structure()
        Add loops to an existing structure that has been passed to the 
        constructor as its seed_pdb argument.
    """

    def __init__(self, seed_pdb, query_pdb, exclusion_pdb,  workdir, para):
        if workdir:
            _workdir = os.path.realpath(workdir)
            if not os.path.exists(_workdir):
                os.mkdir(_workdir)
        else:
            _workdir = os.getcwd() + '/output_' + datetime.now().strftime('%Y-%m-%d-%H-%M-%S')          
            os.mkdir(_workdir)
        self.para = para
        #--------------------------------
        self.log = logger.logger_config(log_path=_workdir + '/log.txt', logging_name='smallprot')
        self.log.info("Creat Smallprot object.")      
        self.log.info("seed_pdb: {}".format(seed_pdb))
        self.log.info("query_pdb: {}".format(query_pdb))
        self.log.info("exclusion_pdb: {}".format(exclusion_pdb))
        smallprot_config.writeConfig(_workdir + '/parameter.ini', self.para)
        #--------------------------------       
        self._prepare_pdbs(seed_pdb, query_pdb, exclusion_pdb,  _workdir)
        self.workdir = _workdir  
        self.loop_workdir = _workdir     
        self.loop_range = [self.para.min_loop_length, self.para.max_loop_length]
        self.targetList = os.path.realpath(self.para.database) + '/pds_list_2p5.txt'
        self.chains_dict = os.path.realpath(self.para.database) + '/db_2p5A_0p3rfree_chains_dictionary.pkl'
        #--------------------------------   
        self.n_truncations = []
        self.c_truncations = []
        self.sat = None
        self.looped_pdbs = []

        #for loop struct infos
        self.infos = []
        self.combs = []
        #--------------------------------

    def _prepare_pdbs(self, seed_pdb, query_pdb, exclusion_pdb, _workdir):
        # if necessary, split query pdb file into chains
        if query_pdb:
            query_chain_dir = _workdir + '/query_chains'
            if not os.path.exists(query_chain_dir):
                os.mkdir(query_chain_dir)
            pdbutils.split_pdb(query_pdb, query_chain_dir, set_bfac=np.log(self.para.min_nbrs))
            self.query_sse_list = [query_chain_dir + '/' + f for f in 
                                   os.listdir(query_chain_dir)
                                   if 'chain_' in f]
            self.query_sse_list.sort()
        else:
            self.query_sse_list = []
        # if necessary, prepare seed pdb files
        _seed_pdb = _workdir + '/seed.pdb'
        if seed_pdb:
            seed_chain_dir = _workdir + '/seed_chains'
            if not os.path.exists(seed_chain_dir):
                os.mkdir(seed_chain_dir)
            pdbutils.split_pdb(seed_pdb, seed_chain_dir)
            self.full_sse_list = [seed_chain_dir + '/' + f for f in 
                                  os.listdir(seed_chain_dir)
                                  if 'chain_' in f]
            self.full_sse_list.sort()
            if query_pdb:
                pdbutils.merge_pdbs([query_pdb, seed_pdb], _seed_pdb)
            else:
                pdbutils.merge_pdbs([seed_pdb], _seed_pdb, set_bfac=np.log(self.para.min_nbrs))      
        elif query_pdb:
            pdbutils.merge_pdbs([query_pdb], _seed_pdb, set_bfac=np.log(self.para.min_nbrs))
            self.full_sse_list = []
        self.seed_pdb = _seed_pdb        
        if not query_pdb and not seed_pdb:
            raise AssertionError('Must provide either query_pdb or seed_pdb.') 
        # if necessary, determine path to exclusion PDB file
        if exclusion_pdb:
            _exclusion_pdb = _workdir + '/exclusion.pdb'
            pdbutils.merge_pdbs([_seed_pdb, exclusion_pdb], _exclusion_pdb, 
                                set_bfac=np.log(self.para.min_nbrs))
            self.orig_exclusion = exclusion_pdb
        else:
            _exclusion_pdb = _seed_pdb
            self.orig_exclusion = None
        self.exclusion_pdb = _exclusion_pdb

    def loop_structure(self, direction=[], n_truncations=[0], c_truncations=[0]):
        if len(self.full_sse_list) == 0:
            raise AssertionError('seed_pdb not provided to constructor.')
        # compute the number of satisfied N- and C-termini
        all_sat = pdbutils.satisfied_termini(self.seed_pdb, self.para.max_nc_dist)
        sat = np.zeros_like(all_sat)
        if len(direction)!=0:          
            for i in range(len(direction)-1):
                j = i+1
                if all_sat[direction[i], direction[j]]:
                    sat[direction[i], direction[j]] = 1
        else:
            sat = all_sat.copy()
        print(sat)
        self.sat = sat
        # set n_truncations and c_truncations for loop generation
        self.n_truncations = n_truncations
        self.c_truncations = c_truncations
        # generate loops
        _full_sse_list = self.full_sse_list.copy()
        # for sse in _full_sse_list:
        #     self.log.info(sse) 
        self._generate_trunc_loops(self.full_sse_list, sat, self.loop_workdir, direction, n_truncations, c_truncations, self.para.cluster_count_cut, self.loop_range)     
        self.log.info('Finish build protein.')

    ### NEW FUNCTIONS FOR GENERATING LOOPS

    def _cal_win_dist(self, full_sse_list, loop_query_win): 
        n_reps = len(full_sse_list)     
        qrep_natoms, qrep_nres, dists = peputils.cal_cdist(full_sse_list)

        #calcualte distances (with window = 7) between two the first sse and another sse, 
        #store the minimum index.
        win_dists = np.zeros((n_reps, qrep_nres[0]- loop_query_win + 1))
        win_dists[0] = np.arange(qrep_nres[0]- loop_query_win + 1)
        for k in range(1, n_reps):
            for x in range(qrep_nres[0]- loop_query_win + 1):
                win_dist = [0]* (qrep_nres[k]- loop_query_win + 1)
                for y in range(qrep_nres[k]- loop_query_win + 1):
                    for z in range(28):
                        win_dist[y] += dists[qrep_natoms[0]+x*4 + z, qrep_natoms[k]+y*4 + z]
                win_dists[k, x]=win_dist.index(min(win_dist))
        return qrep_nres, win_dists  

    def _get_truncs(self, full_sse_list, sat, loop_query_win, n_truncations=[], c_truncations=[]):
        """It may be better to use z index for helix bundle."""
        n_chains = len(full_sse_list)
        qrep_nres, win_dists = self._cal_win_dist(full_sse_list, loop_query_win)           
        truncs = []
        all_local_sats = []
        all_directs = []
        for n in n_truncations:
            if n > qrep_nres[0] - loop_query_win:
                continue
            t = win_dists[:,n]
            truncs.append(t) 
            the_sat = sat.copy()
            for i in range(n_chains):
                if i%2==0:
                    the_sat[i, ] = 0
            all_local_sats.append(the_sat)
            all_directs.append('n' + str(n))

        for c in c_truncations:
            if c > qrep_nres[0] - loop_query_win:
                continue
            t = win_dists[:, -c-1]
            truncs.append(t)
            the_sat = sat.copy()
            for i in range(n_chains):
                if i%2!=0:
                    the_sat[i, ] = 0
            all_local_sats.append(the_sat)
            all_directs.append('c'+ str(c))
        
        return truncs, all_local_sats, all_directs

    def _generate_trunc_loops(self, the_full_sse_list, sat, workdir, direction, n_truncations=[], c_truncations=[], cluster_count_cut=20, loop_range=[3, 20]):
        n_chains = len(sat)
        # find loops for each pair of nearby N- and C-termini
        all_truncs, all_local_sats, all_directs = self._get_truncs(the_full_sse_list, sat, self.para.loop_query_win, n_truncations, c_truncations)
        # for each trunc, find all loops.
        for i in range(len(all_truncs)):
            trunc = all_truncs[i]
            the_sat = all_local_sats[i]
            the_direct = all_directs[i]
            _workdir = workdir +'/trunc_'+ the_direct
            if not os.path.exists(_workdir):
                os.mkdir(_workdir)           
            #Search loops
            self.new_loop_search_fast(the_full_sse_list, the_sat, trunc, _workdir, loop_range)     
            #Summarize loops
            _infos = self._get_top_cluster_summary(_workdir, n_chains, cluster_count_cut, loop_range)
            self.infos.extend(_infos)   

        #Remove redundency. Some loops share the same structure. 
        reduced_order_infos = self._remove_redundancy(cluster_count_cut)
        
        #Write struct info
        self._write_loop_summary(workdir + '/summary_loops.txt', self.infos) 

        self.write_looped_pdb(the_full_sse_list, reduced_order_infos, n_chains, sat, cluster_count_cut)

    def write_looped_pdb(self, the_full_sse_list, reduced_order_infos, n_chains, sat, cluster_count_cut, output_cut = 100):
        #Find loop combine candidates.
        combs, ps, scores = self._extract_loop_combs(reduced_order_infos, n_chains, sat, cluster_count_cut)
        self.combs.extend(combs)
        #print("combs: "+ str(len(combs)))     
        output_cut = output_cut if len(combs) > output_cut else len(combs)
        inds = np.argsort(scores)[::-1][:output_cut]

        looped_pdb_info = []
        for rank in range(len(inds)):
            ind = inds[rank]
            comb = combs[ind]
            p = ps[ind]
            centroids = [c.cent_pdb for c in comb]
            clashing = pdbutils.check_clashes(centroids)
            if clashing:
                continue

            outdir = self.loop_workdir + '/loop_' + '-'.join([str(s) for s in p]) 
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            for c in centroids:
                dst_dir = outdir + '/' + c.split('/')[-1]
                shutil.copy(c, dst_dir)
                shutil.copy(c.split('.')[0] + '_info.png', dst_dir.split('.')[0] + '_info.png')

            structs, slices = self._connect_loops_struct(the_full_sse_list, n_chains, p, centroids, self.para.construct_keep)
            out_pdb = 'output_' + '-'.join(str(v) for v in p) + '_' + str(rank) + '_' + str(scores[ind]) + '.pdb'
            out_pdb_path = outdir + '/' + out_pdb
            pdbutils.merge_save_struct(out_pdb_path, structs, slices) 

            looped_pdb_info.append(out_pdb + '\t' + str(scores[ind]) + '\t' + '\t'.join([os.path.basename(c) for c in centroids]))    
        ##TO DO

        self._write_looped_pdb_summary(self.workdir + '/summary_proteins.txt', looped_pdb_info) 
   
    def new_loop_search_fast(self, _full_sse_list, sat, trunc, workdir, loop_range=[3, 20]):
        """#Find loops for each pair of nearby N- and C-termini. return slice_lengths?"""      
        n_chains = len(sat)
        for j, k in product(range(n_chains), repeat=2):
            # ensure selected SSEs satisfy the distance constraint
            if not sat[j, k] or j == k:
                continue
            print('Generating loops between SSEs {} ' 'and {}'.format(string.ascii_uppercase[j], string.ascii_uppercase[k]))
            loop_workdir = workdir + '/loops_{}_{}'.format(string.ascii_uppercase[j], string.ascii_uppercase[k])
            if not os.path.exists(loop_workdir):
                os.mkdir(loop_workdir)
            loop_query = loop_workdir + '/loop_query.pdb'
            loop_outfile = loop_workdir + '/stdout'
            # calculate how many residues are required for an overlap region 
            # of length 10 Angstroms between the query SSEs and the loops
            inds = [j,k]
            pdbutils.gen_loop_query_win(_full_sse_list, loop_query, inds, trunc, self.para.loop_query_win)
            # find loops with MASTER
            # sort PDBs into directories by loop length
            clusters_exist = self._loop_search_query_search(loop_workdir, loop_query, loop_outfile, loop_range)
            # cluster loops if the clusters do not already exist
            if not clusters_exist:
                cluster_loops.run_cluster(loop_workdir + '/', self.log, self.para.loop_query_win, outfile=loop_outfile)

    def _loop_search_query_search(self, loop_workdir, loop_query, loop_outfile, loop_range):
        """find loops with MASTER"""
        gapLen = str(loop_range[0]) + '-' + str(loop_range[1])
        if not os.path.exists(loop_outfile):
            print('Querying MASTER for loops of length {} to {}.'.format(
                    str(loop_range[0]), str(loop_range[1])))
            query.master_query_loop(loop_query, self.para.loop_target_list, 
                                    rmsdCut=self.para.rmsdCut, topN=self.para.master_query_loop_top,
                                    gapLen=gapLen, outdir=loop_workdir, 
                                    outfile=loop_outfile)
        clusters_exist = True
        loop_workdir_paths = os.listdir(loop_workdir)
        print('Sorting loop PDBs by loop length.')
        # sort PDBs into directories by loop length
        for path in loop_workdir_paths:
            if '.pdb' in path and 'loop_query' not in path:
                with open(loop_workdir + '/' + path, 'r') as f:
                    res_ids = set([int(line[23:26]) for line in 
                                    f.read().split('\n') if 
                                    line[:4] == 'ATOM'])
                    # subtract query ends from loop length
                    l = len(res_ids) - 2*self.para.loop_query_win
                l_dir = loop_workdir + '/' + str(l)
                # create a directory for the loop length if necessary
                if str(l) not in loop_workdir_paths:
                    os.mkdir(l_dir)
                    loop_workdir_paths.append(str(l))
                os.rename(loop_workdir + '/' + path, l_dir + '/' + os.path.basename(path))
            elif os.path.basename(path) in [str(n) for n in range(100)]:
                clusters_path = loop_workdir + '/' + path + '/clusters'
                if not os.path.exists(clusters_path):
                    os.mkdir(clusters_path)
                    clusters_exist = False
        return clusters_exist

    def _get_top_cluster_summary(self, workdir, n_chains, cluster_count_cut, loop_range):
        '''
        After search & cluster loops, calculate the clustered loop info. 
        Copy the centroid or min-rmsd loop into self.workdir/loops_X_Y folder for futher consideration.
        Plot the phi/psi, sequence logo, hydrophobicity.
        '''
        _infos = []
        for p in permutations(range(n_chains), 2):          
            loop_workdir = workdir + '/loops_{}_{}'.format(string.ascii_uppercase[p[0]], string.ascii_uppercase[p[1]])
            for l in range(loop_range[0], loop_range[1] + 1):
                subdir = loop_workdir + '/{}/clusters/1'.format(str(l))                 
                try:
                    loop_pdbs = os.listdir(subdir)
                except:
                    loop_pdbs = []
                if len(loop_pdbs) > 1:       
                    loop_rmsds, loop_seqs, loop_pdss = extract_master._get_pdbs_master_info(loop_workdir + '/match.txt', loop_workdir + '/seq.txt', loop_pdbs)
                    _cent_pdb = ''                  
                    #Copy centroid pdb and plot   
                    if len(loop_pdbs) >= cluster_count_cut:
                        _cent_pdb_workdir = self.loop_workdir + '/loops_{}_{}'.format(string.ascii_uppercase[p[0]], string.ascii_uppercase[p[1]])        
                        _cent = 'loops_{}_{}'.format(string.ascii_uppercase[p[0]], string.ascii_uppercase[p[1]]) + '_' \
                            + workdir.split('/')[-1] + '_rg' + str(l) + '_' + str(len(loop_pdbs))
                        _cent_pdb = _cent_pdb_workdir + '/' + _cent + '.pdb'    
                        if not os.path.exists(_cent_pdb_workdir):
                            os.mkdir(_cent_pdb_workdir)
                        if self.para.select_min_rmsd_pdb:
                            in_pdb = subdir + '/' + loop_pdbs[np.argmin(loop_rmsds)]
                        else:
                            in_pdb = subdir + '/' + [lpdb for lpdb in loop_pdbs if 'centroid' in lpdb][0]

                        with gzip.open(in_pdb, 'rb') as f_in:
                            with open(_cent_pdb, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)                  

                        phi, psi, _sel_seq = struct_analysis.meature_phipsi(_cent_pdb)
                        plot._plot_all(_cent_pdb_workdir + '/' + _cent, loop_seqs, loop_rmsds, l, self.para.loop_query_win, phi, psi, _sel_seq)     

                    #Add summary                
                    info = Struct_info(trunc_info = loop_workdir.split('/')[-2], loop_info = loop_workdir.split('/')[-1], 
                            loop_len = l, clust_num = len(loop_pdbs), min_rmsd = min(loop_rmsds), median_rmsd = np.median(loop_rmsds), cent_pdb = _cent_pdb)              
                    
                    #DO we need to check the size of the 2nd cluster              
                    try:
                        subdir2 = loop_workdir + '/{}/clusters/2'.format(str(l))
                        loop_pdbs2 = os.listdir(subdir2)
                        if len(loop_pdbs2) > 0:
                            info.clust_num_2nd = str(len(loop_pdbs2)) 
                    except:
                        print('Cluster 2 has zero')
              
                    _infos.append(info)      

        return _infos

    def _remove_redundancy(self, cluster_count_cut):
        reduced_order_infos = []
        self.infos.sort(key = lambda x: x.clust_num, reverse = True) 
        redundant_inds = set()
        for i in range(len(self.infos)-1):
            if i in redundant_inds or self.infos[i].cent_pdb=='':             
                continue
            self.infos[i].redundancy = 'New'
            reduced_order_infos.append(self.infos[i])
            for j in range(i+1, len(self.infos)-1):
                if self.infos[j].cent_pdb=='':
                    continue
                if j in redundant_inds:                 
                    continue
                if self.infos[i].loop_len <= self.infos[j].loop_len:
                    #TO DO: The distance calculation is not the real distance.
                    min_dist, min_dist_ind = peputils.cal_sse_dist([self.infos[i].cent_pdb, self.infos[j].cent_pdb])
                else:
                    min_dist, min_dist_ind = peputils.cal_sse_dist([self.infos[j].cent_pdb, self.infos[i].cent_pdb])               
                if min_dist < self.para.rmsdCut: 
                    self.infos[j].redundancy = self.infos[i].trunc_info + " " + self.infos[i].loop_info + " " + str(self.infos[i].loop_len)         
                    redundant_inds.add(j)

        # outdir = self.loop_workdir + '/reduced_loops'
        # if not os.path.exists(outdir):
        #     os.mkdir(outdir)
        # for v in reduced_order_infos:
        #     dst_dir = outdir + '/' + v.cent_pdb.split('/')[-1]
        #     shutil.copy(v.cent_pdb,dst_dir)

        return reduced_order_infos

    def _extract_loop_combs(self, reduced_order_infos, n_chains, sat, cluster_count_cut):
        '''
        With all the loops with >= cluster_count_cut. 
        Try to find the combination so that they can loop the input sse.
        For example, combination (loop_A_B, loop_B_C, loop_C_D) could loop (A, B, C, D) with direction [0, 1, 2, 3]
        '''
        combs = []
        ps = []
        for p in permutations(range(n_chains)):
            if not np.all([sat[p[j], p[j+1]] for j in range(n_chains - 1)]):
                continue
            all_keys = []
            for j in range(n_chains-1):
                k = j+1
                loop_key = 'loops_{}_{}'.format(string.ascii_uppercase[p[j]], string.ascii_uppercase[p[k]])
                keys = [key for key in reduced_order_infos if loop_key in key.loop_info]
                values = [v.clust_num for v in keys]
                all_keys.append([keys[i] for i in np.argsort(values) if values[i] > cluster_count_cut])
            if 0 in [len(v) for v in all_keys]:
                continue
            for comb in product(*all_keys):
                print(comb)
                if self._check_comb_validity(comb, self.para.loop_distance_cut):
                    continue
                combs.append(comb)
                ps.append(p)
        scores = []
        for comb in combs:
            scores.append(sum([c.clust_num for c in comb]))
        return combs, ps, scores

    def _check_comb_validity(self, comb, loop_distance_cut = 15):
        """#Currently, if the distance between mid residues of two loops on the same sides should be smaller than a certain number."""
        n_side = comb[0::2]
        c_side = comb[1::2]
        if len(n_side) >= 2:
            for j, k in product(range(len(n_side)), repeat = 2):
                dist = peputils.cal_loop_mid_dist([n_side[j].cent_pdb, n_side[k].cent_pdb])
                if dist > loop_distance_cut:
                    return False
        if len(c_side) >= 2:
            for j, k in product(range(len(c_side)), repeat = 2):
                dist = peputils.cal_loop_mid_dist([c_side[j].cent_pdb, c_side[k].cent_pdb])
                if dist > loop_distance_cut:
                    return False
        return True
 
    def _connect_loops_struct(self, _full_sse_list, n_chains, permutation, centroids, keep=1):
        '''
        Put seed sses and loop candidates together and calculate the cut position for each.      
        @ _full_sse_list: [string], list of sse path.
        @ permutation: [int], list of direaction. for example [0, 1, 2, 3] means connect sses with order of [A, B C, D].
        @ keep: keep == 0, keep min; keep == -1, keep loop; keep == 1, keep seed. 
        '''
        pdbs_to_combine = [''] * (2 * n_chains - 1)
        pdbs_to_combine[::2] = [_full_sse_list[idx] for idx in permutation]
        pdbs_to_combine[1::2] = centroids

        inds = [[0,0] for i in range(len(pdbs_to_combine))]
        for i in range(len(pdbs_to_combine)-1):
            j = i + 1
            # 
            min_dist, min_dist_ind, qrep_nres = peputils.cal_aa_dist([pdbs_to_combine[i], pdbs_to_combine[j]], i%2==0, self.para.loop_query_win, keep)       
            if i%2==0:          
                inds[i][1] = min_dist_ind[0]
                inds[j][0] = min_dist_ind[1]
            else:
                inds[i][1] = min_dist_ind[0]                
                inds[j][0] = min_dist_ind[1]
                inds[j][1] = qrep_nres[1]
        #print('inds')
        #print(inds)                     
        slices = [] 
        for d in inds:
            slices.append(slice(d[0], d[1]))
        structs = [pdbutils.get_struct('test_' + str(i), pdbs_to_combine[i]) for i in range(len(pdbs_to_combine))]
        return structs, slices
      
    def _write_loop_summary(self, filename, loop_infos):
        '''
        Write information of all loops.
        @ loop_infos: [Struct_info]
        '''
        with open(filename, 'w') as f:
            f.write('trunc_info\tloop_info\tloop_len\tclust_num\tmin_rmsd\tmedian_rmsd\tredundancy\tcent_pdb\tclust_num_2nd\n')
            for r in loop_infos:
                f.write(r.trunc_info + '\t' + r.loop_info + '\t' + str(r.loop_len) + '\t'
                    + str(r.clust_num) + '\t'+ str(r.min_rmsd) + '\t'+ str(r.median_rmsd)+ '\t'+ r.redundancy+ '\t' 
                    + r.cent_pdb + '\t' + str(r.clust_num_2nd) + '\n')       
 
    def _write_looped_pdb_summary(self, filename, looped_pdb_infos):
        '''
        Write the name of the looped_pdb, total score, name of each loop.
        @ looped_pdb_infos: [string], each string contains info of one looped_pdb.
        '''
        with open(filename, 'w') as f:
            f.write('filename\ttotal_score\tloops\n')
            for r in looped_pdb_infos:
                f.write(r + '\n')   