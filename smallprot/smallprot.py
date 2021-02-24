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

from smallprot import pdbutils, query, cluster_loops, smallprot_config, logger, peputils, constant

@dataclass
class Struct_info:
    trunc_info: str
    loop_info: str
    loop_len: int 
    cent_pdb:str
    clust_num: int   
    redundancy:str = "Unknown" 
    clust_num_2nd: int = 0

class SmallProt:
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

    def __init__(self, seed_pdb, query_pdb, exclusion_pdb,  workdir, para_file_path = 'parameter.ini'):
        if os.path.exists(para_file_path):
            self.para = smallprot_config.readConfig(para_file_path)
        else:
            self.para = smallprot_config.Parameter()
        if workdir:
            _workdir = os.path.realpath(workdir)
            if not os.path.exists(_workdir):
                os.mkdir(_workdir)
        else:
            _workdir = os.getcwd() + '/output_' + datetime.now().strftime('%Y-%m-%d-%H-%M-%S')          
            os.mkdir(_workdir)
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
        self.looped_pdbs = []
        #self.output_pdbs = []

        #For job queue, to parallel jobs
        self.queues = []

        #for loop struct infos
        self.infos = []
        self.combs = []
        #--------------------------------
        self.pre_build_pdbs = []
        self.pre_full_sses = []

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

    def build_protein(self, n_truncations=[0], c_truncations=[0]):
        """Iteratively generate a protein using MASTER and Qbits."""
        self.log.info('Start build protein.')
        self.n_truncations = n_truncations
        self.c_truncations = c_truncations
        self._generate_proteins()
        # print('output pdbs :')
        # print('\n'.join(self.output_pdbs))
        self.log.info('Finish build protein.')
  
    def build_protein_parallel(self):
        """!!!Not working due to conflict with qbit."""
        """Iteratively generate a protein using MASTER and Qbits."""
        self._generate_proteins_parallel(self.para.num_iter)
        # print('output pdbs :')
        # print('\n'.join(self.output_pdbs))
        self.log.info('Finish build protein.')

    def loop_structure(self, direction=[], n_truncations=[0], c_truncations=[0]):
        if len(self.full_sse_list) == 0:
            raise AssertionError('seed_pdb not provided to constructor.')
        # compute the number of satisfied N- and C-termini
        all_sat = pdbutils.satisfied_termini(self.seed_pdb, self.para.max_nc_dist)
        print(all_sat)
        sat = np.zeros_like(all_sat)
        print(len(direction))
        if len(direction)!=0:          
            for i in range(len(direction)-1):
                j = i+1
                if all_sat[direction[i], direction[j]]:
                    sat[direction[i], direction[j]] = 1
        else:
            sat = all_sat.copy()
        print(sat)
        # set n_truncations and c_truncations for loop generation
        self.n_truncations = n_truncations
        self.c_truncations = c_truncations
        # generate loops
        _full_sse_list = self.full_sse_list.copy()
        for sse in _full_sse_list:
            self.log.info(sse) 
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
        self._write_file(workdir + '/summary.txt', self.infos) 

        #Find loop combine candidates.
        combs, ps, scores = self._extract_top_hit(reduced_order_infos, n_chains, sat, cluster_count_cut)
        self.combs.extend(combs)

        # #Build whole structure.
        outdir = workdir + '/output'
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        
        output_cut = 100 if len(combs) > 100 else len(combs)
        inds = np.argsort(scores)[::-1][:output_cut]
        for i in inds:
            comb = combs[i]
            p = ps[i]
            centroids = [c.cent_pdb for c in comb]
            clashing = pdbutils.check_clashes(centroids)
            if clashing:
                continue
            structs, slices = self._connect_loops_struct(the_full_sse_list, n_chains, p, centroids)
            out_path = outdir + '/output_' + '-'.join(str(v) for v in p) + '_' + str(i) + '.pdb'
            pdbutils.merge_save_struct(out_path, structs, slices)      
   
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
                cluster_loops.run_cluster(loop_workdir + '/', outfile=loop_outfile)

    def _loop_search_query_search(self, loop_workdir, loop_query, loop_outfile, loop_range):
        """find loops with MASTER"""
        gapLen = str(loop_range[0]) + '-' + str(loop_range[1])
        if not os.path.exists(loop_outfile):
            print('Querying MASTER for loops of length {} to {}.'.format(
                    str(loop_range[0]), str(loop_range[1])))
            query.master_query_loop(loop_query, self.para.loop_target_list, 
                                    rmsdCut=self.para.rmsdCut, topN=200,
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
                    l = len(res_ids) - 14
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
                    #Copy centroid pdb
                    _cent_pdb = ''
                    _cent_pdb_workdir = self.loop_workdir + '/loops_{}_{}'.format(string.ascii_uppercase[p[0]], string.ascii_uppercase[p[1]])
                    if not os.path.exists(_cent_pdb_workdir):
                        os.mkdir(_cent_pdb_workdir)
                    if len(loop_pdbs) >= cluster_count_cut:
                        _cent = _cent_pdb_workdir + '/loops_{}_{}'.format(string.ascii_uppercase[p[0]], string.ascii_uppercase[p[1]]) + '_' \
                            + workdir.split('/')[-1] + '_rg' + str(l) + '_' + str(len(loop_pdbs))
                        _cent_pdb = _cent + '.pdb'
                        # print(_cent_pdb)
                        # print([subdir + '/' + lpdb for lpdb in loop_pdbs if 'centroid' in lpdb][0])
                        in_pdb = [subdir + '/' + lpdb for lpdb in loop_pdbs if 'centroid' in lpdb][0]
                        with gzip.open(in_pdb, 'rb') as f_in:
                            with open(_cent_pdb, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)  
                    
                        self._plot_log(loop_workdir + '/seq.txt', l, _cent)    
                        self._plot_hydro(loop_workdir + '/seq.txt', l, _cent, self.para.loop_query_win)                  
                        phipsi = pdbutils.meaure_phipsi(_cent_pdb)
                        self._plot_phipsi(phipsi, _cent)    
                        self._plot_propensity(loop_workdir + '/seq.txt', l, _cent, self.para.loop_query_win)
                        
                    #Add summary
                    info = Struct_info(trunc_info = loop_workdir.split('/')[-2], loop_info = loop_workdir.split('/')[-1], 
                            loop_len = l, clust_num = len(loop_pdbs), cent_pdb = _cent_pdb)              
                    #DO we need to check the size of the 2nd cluster
                    subdir = loop_workdir + '/{}/clusters/2'.format(str(l))
                    try:
                        loop_pdbs2 = os.listdir(subdir)
                    except:
                        loop_pdbs2 = []
                    if len(loop_pdbs2) > 0:
                        info.clust_num_2nd = str(len(loop_pdbs2))               
                    _infos.append(info)
        return _infos

    def _remove_redundancy(self, cluster_count_cut):
        reduced_order_infos = []
        self.infos.sort(key = lambda x: x.clust_num, reverse = True) 
        redundant_inds = []
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
                    min_dist, min_dist_ind = peputils.cal_sse_dist([self.infos[i].cent_pdb, self.infos[j].cent_pdb])
                else:
                    min_dist, min_dist_ind = peputils.cal_sse_dist([self.infos[j].cent_pdb, self.infos[i].cent_pdb])               
                self.log.info(str(min_dist))
                if min_dist < self.para.rmsdCut: 
                #if min_dist < 100:        
                    self.infos[j].redundancy = self.infos[i].trunc_info + " " + self.infos[i].loop_info + " " + str(self.infos[i].loop_len)         
                    redundant_inds.append(j)

        outdir = self.loop_workdir + '/reduced_loops'
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        for v in reduced_order_infos:
            dst_dir = outdir + '/' + v.cent_pdb.split('/')[-1]
            shutil.copy(v.cent_pdb,dst_dir)

        return reduced_order_infos

    def _extract_top_hit(self, reduced_order_infos, n_chains, sat, cluster_count_cut):
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
 
    def _connect_loops_struct(self, _full_sse_list, n_chains, permutation, centroids):
        """Find the min distance aa pair to connect."""
        pdbs_to_combine = [''] * (2 * n_chains - 1)
        pdbs_to_combine[::2] = [_full_sse_list[idx] for idx in permutation]
        pdbs_to_combine[1::2] = centroids

        inds = [[0,0] for i in range(len(pdbs_to_combine))]
        for i in range(len(pdbs_to_combine)-1):
            j = i + 1
            min_dist, min_dist_ind, qrep_nres = peputils.cal_aa_dist([pdbs_to_combine[i], pdbs_to_combine[j]])
            #print('min_dist_ind')
            #print(min_dist_ind)
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
      
    def _write_file(self, filename, infos):
        with open(filename, 'w') as f:
            f.write('trunc_info\tloop_info\tloop_len\tclust_num\tcent_pdb\tredundancy\tclust_num_2nd\n')
            for r in infos:
                f.write(r.trunc_info + '\t' + r.loop_info + '\t' + str(r.loop_len) + '\t'
                    + str(r.clust_num) + '\t' + r.cent_pdb + '\t'+ r.redundancy+ '\t' + str(r.clust_num_2nd) + '\n')       
    
    def _plot_log(self, seqfile, seqlen, filepath):
        with open(seqfile, 'r') as f:
            lines = f.read().split('\n')
        all_seqs = []
        for line in lines:
            if len(line) > 0:
                seq = ''
                for res in line.split(' '):
                    if len(res) == 3 and res[0].isalpha():
                        seq += qbits.constants.one_letter_code[res]
                    elif len(res) == 4 and res[0] == '[':
                        seq += qbits.constants.one_letter_code[res[1:]]
                    elif len(res) == 4 and res[-1] == ']':
                        seq += qbits.constants.one_letter_code[res[:-1]]
                all_seqs.append(seq)
        seqs = []
        for s in all_seqs:
            if len(s) == 14 + seqlen:            
                seqs.append(s)

        df = logomaker.alignment_to_matrix(sequences=seqs, to_type='counts',
                                               characters_to_ignore='-', pseudocount=0.01)
        logo = logomaker.Logo(df,
                         font_name='Arial',
                         color_scheme='NajafabadiEtAl2017',
                         vpad=.1,
                         width=.8)
        logo.style_xticks(anchor=0, spacing=1)      
        logo.ax.set_ylabel('Count')
        logo.ax.set_xlim([-1, len(df)])
        #logo.fig.savefig(filepath) 
        plt.savefig(filepath + '_logo.png')
        plt.close()

    def _plot_hydro(self, seqfile, seqlen, filepath, loop_query_win):
        with open(seqfile, 'r') as f:
            lines = f.read().split('\n')
        all_hydro = []
        for line in lines:
            if len(line) > 0:
                hydro = []
                for res in line.split(' '):
                    if len(res) == 3 and res[0].isalpha():
                        hydro.append(constant.hydro_dict[res]) 
                    elif len(res) == 4 and res[0] == '[':
                        hydro.append(constant.hydro_dict[res[1:]])
                    elif len(res) == 4 and res[-1] == ']':
                        hydro.append(constant.hydro_dict[res[:-1]])
                if len(hydro) == seqlen + 2*loop_query_win:
                    all_hydro.append(hydro)
        all_hydro_arr = np.array(all_hydro)

        means = np.mean(all_hydro_arr, 0)
        sds = np.std(all_hydro_arr, 0)

        x = list(range(1, seqlen + 2*loop_query_win+1))
        for i in range(len(constant.hydro_scale)):
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111)            
            ax.set_xlabel('AA', fontsize = 18)
            ax.set_ylabel(constant.hydro_scale[i], fontsize = 18)
            #for i in range(len(constant.hydro_scale)):
            #    ax.errorbar(x, means[:, i], yerr = sds[:, i], label = constant.hydro_scale[i])
            ax.errorbar(x, means[:, i], yerr = sds[:, i])
            plt.savefig(filepath+'_hydro_'+str(i)+'.png')
            plt.close()

    def _plot_propensity(self, seqfile, seqlen, filepath, loop_query_win):
        with open(seqfile, 'r') as f:
            lines = f.read().split('\n')
        all_propen = []
        for line in lines:
            if len(line) > 0:
                propen = []
                for res in line.split(' '):
                    if len(res) == 3 and res[0].isalpha():
                        propen.append(constant.propensity_dict[res]) 
                    elif len(res) == 4 and res[0] == '[':
                        propen.append(constant.propensity_dict[res[1:]])
                    elif len(res) == 4 and res[-1] == ']':
                        propen.append(constant.propensity_dict[res[:-1]])
                if len(propen) == seqlen + 2*loop_query_win:
                    all_propen.append(propen)
        all_propen_arr = np.array(all_propen)
        means = np.mean(all_propen_arr, 0)
        sds = np.std(all_propen_arr, 0)

        x = list(range(1, seqlen + 2*loop_query_win+1))
        #for i in range(len(constant.propensity_scale)):
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)            
        ax.set_xlabel('AA', fontsize = 18)
        ax.set_ylabel(constant.propensity_scale[0], fontsize = 18)

        ax.errorbar(x, means[:, 0], yerr = sds[:, 0])
        plt.savefig(filepath+'_helix_propen.png')
        plt.close()

    def _plot_phipsi(self, phipsi, filepath):
        x = list(range(1, int(len(phipsi)/2+1)))
        plt.plot(x, phipsi[::2], label = 'phi')
        plt.plot(x, phipsi[1::2], label = 'psi')
        plt.legend()
        plt.savefig(filepath + '_phipsi.png')
        plt.close()

        plt.plot(phipsi[2:-2:2], phipsi[3:-1:2], 's', color='red', markersize=5, markerfacecolor='white')
        plt.xlim(-180, 180)
        plt.ylim(-180, 180)
        xticks = [-180, -135, -90, -45, 0, 45, 90, 135, 180]
        yticks = [-180, -135, -90, -45, 0, 45, 90, 135, 180]
        plt.xticks(xticks)
        plt.yticks(yticks)
        plt.savefig(filepath + '_ramachandran.png')
        plt.close()
  
    ### FUNCTIONS FOR GENERATING SSEs

    def _generate_proteins(self):     
        self.queues.append([self.seed_pdb, self.exclusion_pdb, self.full_sse_list, self.para.num_iter])
        while len(self.queues) > 0:            
            queue = self.queues.pop(0)
            print('--Get queue--')
            self._const_protein(queue[0], queue[1], queue[2], queue[3]) 

        print('Finish construct sses!')
        if len(self.pre_build_pdbs) > 0:                      
            for i in range(len(self.pre_build_pdbs)):
                loop_build_path = os.path.dirname(self.pre_build_pdbs[i]) + '/' + os.path.splitext(self.pre_build_pdbs[i])[0].split('/')[-1]
                if not os.path.exists(loop_build_path):
                    os.mkdir(loop_build_path)
                self.loop_workdir = loop_build_path
                self._const_prot_loop(self.pre_build_pdbs[i], self.pre_full_sses[i], loop_build_path)
            print('Construct sses found.')
        else:
            print('No construct sses found.')
            return
        
    def _cut_pdbs(self, full_sse_list, loop_query_win):
        n_reps = len(full_sse_list)     
        qrep_natoms, qrep_nres, dists = peputils.cal_cdist(full_sse_list)

        win_dists = np.zeros((n_reps, 2), dtype = int)
        win_dists[0] = [0, qrep_nres[0]- loop_query_win]
        for k in range(1, n_reps):
            for x in range(2):
                win_dist = [0]* (qrep_nres[k]- loop_query_win + 1)
                for y in range(qrep_nres[k]- loop_query_win + 1):
                    for z in range(28):
                        win_dist[y] += dists[qrep_natoms[0]+win_dists[0][x]*4 + z, qrep_natoms[k]+y*4 + z]
                win_dists[k, x]=win_dist.index(min(win_dist))

        slices = []
        for i in range(n_reps):
            #if i%2 == 0:
            if win_dists[i][0] < win_dists[i][1]:
                slices.append(slice(win_dists[i][0], win_dists[i][1] + loop_query_win))
            else:
                slices.append(slice(win_dists[i][1], win_dists[i][0] + loop_query_win))

        structs = [pdbutils.get_struct(str(i), full_sse_list[i]) for i in range(n_reps)]
        struct = pdbutils.merge_structs(structs, slices)
        return struct

    def _const_protein(self, pdb, exclusion_pdb, full_sse_list, recursion_order):
        outdir = os.path.dirname(pdb)
        #Construct final protein.
        if recursion_order == 0:
            folders = outdir.split('/')
            filename = 'full'
            for i in range(self.para.num_iter, 0, -1):
                filename += '_' + folders[-i]
            filename += '.pdb'
            _full_pdb = outdir + '/' + filename
            pdbutils.merge_pdbs(full_sse_list, _full_pdb, min_nbrs=self.para.min_nbrs)
            #struct = self._cut_pdbs(full_sse_list, self.para.loop_query_win)
            #pdbutils.save_struct(struct, _full_pdb)
            pre_build_workdir = self.workdir + '/pre_build_pdbs'
            if not os.path.exists(pre_build_workdir):
                os.mkdir(pre_build_workdir)
            _full_pdb_new = pre_build_workdir + '/' + filename
            shutil.copyfile(_full_pdb, _full_pdb_new)

            self.pre_build_pdbs.append(_full_pdb_new)
            self.pre_full_sses.append(full_sse_list)
            #self._const_prot_loop(_full_pdb, exclusion_pdb, full_sse_list, outdir)    
            return
        #Generate queries for next recursion.
        origin_qreps = self._generate_qreps(pdb, exclusion_pdb, recursion_order, outdir)
        qreps = self._resize_qreps(origin_qreps, full_sse_list, self.para.loop_query_win)
        if qreps == None:
            return
        for i, qrep in enumerate(qreps):
            seed_sse_lists = self._prepare_seed_sses(qrep, full_sse_list)
            for j, seed_sse_list in enumerate(seed_sse_lists):
                self._add_seed_sse(i, qrep, j, seed_sse_list, full_sse_list, exclusion_pdb, recursion_order, outdir)
        return
   
    def _const_prot_loop(self, pdb, full_sse_list, outdir):
        sat = pdbutils.satisfied_termini(pdb, self.para.max_nc_dist)
        n_sat = np.sum(sat)
        # if self.para.num_iter - recursion_order - 2 > n_sat:
        #     # if it is impossible to satisfy all N- or C- termini within 
        #     # the remaining number of iterations, exit the branch early
        #     return
        if n_sat >= self.para.num_iter:
            try_loopgen = False
            n_chains = len(sat)
            for p in permutations(range(n_chains)):
                if np.all([sat[p[k], p[k+1]] for k in range(n_chains - 1)]):
                    try_loopgen = True
            if try_loopgen:
                #I don't think this could avoid build same proteins. The chains of the protein have different order. 
                with open(pdb, 'r') as f0:
                    f0_read = f0.read()
                    for pdb in self.looped_pdbs:
                        with open(pdb, 'r') as f1:
                            if f0_read == f1.read():
                                try_loopgen = False 
            # if necessary, ensure the compactness criterion is met
            if self.para.screen_compactness:
                compactness = pdbutils.calc_compactness(pdb)
                try_loopgen = try_loopgen and (compactness > 0.1)
            if try_loopgen:
                self.looped_pdbs.append(pdb)                             
                self._generate_trunc_loops(full_sse_list, sat, outdir, None, self.n_truncations, self.c_truncations, self.para.cluster_count_cut, self.loop_range)

    def _generate_qreps(self, pdb, exclusion_pdb, recursion_order, outdir):
        print('Adding a qbit rep.')
        # search for a contiguous secondary structural element to add
        outfile = outdir + '/stdout'
        print('Querying MASTER')
        query.master_query(pdb, self.targetList, self.para.rmsdCut, 
            topN=None, outfile=outfile, clobber=False)
        print('Searching with Qbits')
        if not os.path.exists(outdir + '/qbit_reps/'):
            try:
                # ensure the second SSE is antiparallel to the first
                #query_exists = int(bool(len(self.query_sse_list)))
                #first_recursion = (recursion_order == self.para.num_iter - query_exists)

                first_recursion = (self.para.num_iter == recursion_order)
                print('first_recursion: ' + str(self.para.num_iter)  +  '-' + str(recursion_order) + str(first_recursion))
                query.qbits_search(pdb, exclusion_pdb, 
                                   self.chains_dict, outdir, 
                                   self.para.qbits_window, self.para.qbits_rmsd, 
                                   top=5, sec_struct=self.para.secstruct,
                                   antiparallel=first_recursion,
                                   min_nbrs=self.para.min_nbrs, contiguous=True)
            except:
                pass
        qreps = None
        if os.path.exists(outdir + '/qbit_reps/'):
            qreps = [outdir + '/qbit_reps/' + pdb_path for pdb_path in os.listdir(outdir + '/qbit_reps/')]
        return qreps

    def _resize_qreps(self, qreps, full_sse_list, loop_query_win):
        resize_qreps = []
        for qrep in qreps:
            two_sses = [full_sse_list[0], qrep]
            struct = self._cut_pdbs(two_sses, loop_query_win)
            chain = list(struct.get_chains())[-1]     
            re_qrep_path = os.path.dirname(qrep) + '_resize'
            if not os.path.exists(re_qrep_path):    
                os.mkdir(re_qrep_path)
            re_qrep = re_qrep_path + '/resize_' + os.path.basename(qrep)
            pdbutils.save_struct(chain, re_qrep)
            resize_qreps.append(re_qrep)
        return resize_qreps


    def _prepare_seed_sses(self, qrep, full_sse_list):
        # if using an external query structure, include it in the Qbits 
        # search until the protein under construction has at least two SSEs
        # if len(full_sse_list) < 2:
        #     all_reps = self.query_sse_list + full_sse_list + [qrep]
        # else:
        #     all_reps = full_sse_list + [qrep]


        all_reps = self.query_sse_list + full_sse_list + [qrep]
        return [all_reps]
        # print('---------------------------')
        # print(qrep)
        # print(all_reps)
        # print(len(all_reps))
        # print('---------------------------')
        if len(all_reps) < 3:
            seed_sse_lists = [all_reps]
        else:
            # extract atomic coordinates of each SSE
            seed_sse_lists = []
            qrep_xyz = []
            qrep_natoms = [0]
            for j, pair_qrep in enumerate(all_reps):
                title = 'struct{}'.format(str(j))
                this_struct = pdbutils.get_struct(title, pair_qrep, self.para.min_nbrs)
                atoms = this_struct.get_atoms()
                qrep_xyz.append(np.array([atom.get_coord() for atom in atoms]))
                qrep_natoms.append(qrep_natoms[-1] + len(qrep_xyz[-1]))
            qrep_xyz = np.vstack(qrep_xyz)
            # compute minimum interatomic distance between SSE pairs
            dists = cdist(qrep_xyz, qrep_xyz)
            n_reps = len(all_reps)
            for j in range(0, n_reps - 1):
                for k in range(j + 1, n_reps):
                    min_dist = np.min(dists[qrep_natoms[j]:qrep_natoms[j+1], qrep_natoms[k]:qrep_natoms[k+1]])
                    # add pairs of SSEs if they are adjacent in space
                    if min_dist < 5.:
                        seed_sse_lists.append([all_reps[j], all_reps[k]])
        return seed_sse_lists


    def _add_seed_sse(self, i, qrep, j, seed_sse_list, full_sse_list, exclusion_pdb, recursion_order, outdir):
        """Add new queries into the self.queues."""
        _workdir = '{}/{}'.format(outdir, str(i) + string.ascii_lowercase[j])
        if not os.path.exists(_workdir):
            os.mkdir(_workdir)
        _seed_pdb = _workdir + '/seed.pdb'
        pdbutils.merge_pdbs(seed_sse_list, _seed_pdb, min_nbrs=self.para.min_nbrs)
        _full_sse_list = full_sse_list.copy()
        _full_sse_list.append(qrep)
        print('SSE List:')
        print('\n'.join(_full_sse_list))
        # compute the number of satisfied N- and C-termini
        #_the_pdb = ''
        # if len(_full_sse_list) > 2:
        #     _full_pdb = _workdir + '/full.pdb'
        #     pdbutils.merge_pdbs(_full_sse_list, _full_pdb, min_nbrs=self.para.min_nbrs)
        #     sat = pdbutils.satisfied_termini(_full_pdb, self.para.max_nc_dist)
        #     _the_pdb = _full_pdb
        # else:
        #     sat = pdbutils.satisfied_termini(_seed_pdb, self.para.max_nc_dist)
        #     _the_pdb = _seed_pdb
        sat = pdbutils.satisfied_termini(_seed_pdb, self.para.max_nc_dist)
        n_sat = np.sum(sat)
        if self.para.num_iter - recursion_order + 1 > n_sat:
            # if it is impossible to satisfy all N- or C- termini within 
            # the remaining number of iterations, exit the branch early
            return

        if recursion_order > 0:
            _exclusion_pdb = _workdir + '/exclusion.pdb'
            pdbutils.merge_pdbs([exclusion_pdb, qrep], _exclusion_pdb, min_nbrs=self.para.min_nbrs)
            self.queues.append([_seed_pdb, _exclusion_pdb, _full_sse_list, recursion_order - 1])
        return    
            
    ### NEW FUNCTIONS FOR GENERATING SSEs IN PARALLEL

    def _generate_proteins_parallel(self, recursion_order):
        manager = Manager()
        pool = Pool(6)
        queues = manager.Queue()

        counter = 1
        item = [self.seed_pdb, self.exclusion_pdb, self.full_sse_list, recursion_order]
        pool.apply_async(self._const_protein_parallel, (item, queues))

        while counter > 0:           
            _item = queues.get(block = True)
            print(_item)
            if _item == 'empty':
                counter -= 1
            else:
                counter += 1
                pool.apply_async(self._const_protein_parallel, (_item, queues))
        print('---out while loop---')
        pool.close()
        pool.join()

    def _const_protein_parallel(self, item, queues):
        """
        The current parallel is not working due to NULL 'qreps = self._generate_qreps(pdb, exclusion_pdb, recursion_order, outdir)'
        query.py -> qbits_search() -> 'p.parse(outdir=outdir, show_parsing_progress=False)'
        p.parse also apply parallel, which may cause a conflict.
        """   
        pdb = item[0]
        exclusion_pdb = item[1]
        full_sse_list = item[2]
        recursion_order = item[3]
        #print('Run task (%s)...' % (os.getpid()))
        outdir = os.path.dirname(pdb)
        #Construct final protein.
        if recursion_order == 0:
            self._const_prot_loop(pdb, exclusion_pdb, full_sse_list, recursion_order, outdir) 
            queues.put('empty', block = True)      
            print('---_const_protein_parallel: first---')
            return
        #Generate queries for next recursion.
        qreps = self._generate_qreps(pdb, exclusion_pdb, recursion_order, outdir)
        if qreps == None:
            print('---_const_protein_parallel: second---')
            queues.put('empty', block = True)    
            return
        _queues = []
        for i, qrep in enumerate(qreps):
            seed_sse_lists = self._prepare_seed_sses(qrep, full_sse_list)
            for j, seed_sse_list in enumerate(seed_sse_lists):
                _queues = self._add_seed_sse(i, qrep, j, seed_sse_list, full_sse_list, recursion_order, outdir)
        for _item in _queues:
            queues.put(_item, block = True)        
        queues.put('empty', block = True)    
        print('---_const_protein_parallel: add more---')
        return
