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

from scipy.stats import mode
from scipy.spatial.distance import cdist
from itertools import product, permutations

import qbits

from smallprot import pdbutils, query, cluster_loops, smallprot_config, logger

@dataclass
class Struct_info:
    trunc_info: str
    loop_info: str
    loop_len: int
    clust_num: int
    cent_pdb:str
    cluster_key_res : []
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

    def __init__(self, para_file_path = 'parameter.ini'):
        if os.path.exists(para_file_path):
            self.para = smallprot_config.readConfig(para_file_path)
        else:
            self.para = Parameter()
        if self.para.workdir:
            _workdir = os.path.realpath(self.para.workdir)
            if not os.path.exists(_workdir):
                os.mkdir(_workdir)
        else:
            _workdir = os.getcwd() + '/output_' + datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
            self.para.workdir = _workdir
            os.mkdir(_workdir)

        self.log = logger.logger_config(log_path=_workdir + '/log.txt', logging_name='smallprot')
        self.log.info("Creat Smallprot object.")

        _seed_pdb = _workdir + '/seed.pdb'
        # if necessary, split query pdb file into chains
        if self.para.query_pdb:
            query_chain_dir = _workdir + '/query_chains'
            if not os.path.exists(query_chain_dir):
                os.mkdir(query_chain_dir)
            pdbutils.split_pdb(self.para.query_pdb, query_chain_dir, set_bfac=np.log(self.para.min_nbrs))
            self.query_sse_list = [query_chain_dir + '/' + f for f in 
                                   os.listdir(query_chain_dir)
                                   if 'chain_' in f]
            self.query_sse_list.sort()
        else:
            self.query_sse_list = []
        # if necessary, prepare seed pdb files
        if self.para.seed_pdb:
            if self.para.query_pdb:
                pdbutils.merge_pdbs([self.para.query_pdb, self.para.seed_pdb], _seed_pdb)
            else:
                pdbutils.merge_pdbs([self.para.seed_pdb], _seed_pdb, set_bfac=np.log(self.para.min_nbrs))
            seed_chain_dir = _workdir + '/seed_chains'
            if not os.path.exists(seed_chain_dir):
                os.mkdir(seed_chain_dir)
            pdbutils.split_pdb(self.para.seed_pdb, seed_chain_dir)
            self.full_sse_list = [seed_chain_dir + '/' + f for f in 
                                  os.listdir(seed_chain_dir)
                                  if 'chain_' in f]
            self.full_sse_list.sort()
        elif self.para.query_pdb:
            pdbutils.merge_pdbs([self.para.query_pdb], _seed_pdb, set_bfac=np.log(self.para.min_nbrs))
            self.full_sse_list = []
        else:
            raise AssertionError('Must provide either query_pdb or seed_pdb.')            
        self.pdbs = [_seed_pdb]
        # if necessary, determine path to exclusion PDB file
        if self.para.exclusion_pdb:
            _exclusion_pdb = _workdir + '/exclusion.pdb'
            pdbutils.merge_pdbs([_seed_pdb, self.para.exclusion_pdb], _exclusion_pdb, 
                                set_bfac=np.log(self.para.min_nbrs))
            self.orig_exclusion = self.para.exclusion_pdb
        else:
            _exclusion_pdb = _seed_pdb
            self.orig_exclusion = None
        # set remaining attributes
        self.workdir = _workdir
        self.exclusion_pdbs = [_exclusion_pdb]
        self.loop_range = [self.para.min_loop_length, self.para.max_loop_length]
        self.targetList = os.path.realpath(self.para.database) + '/pds_list_2p5.txt'
        self.chains_dict = os.path.realpath(self.para.database) + '/db_2p5A_0p3rfree_chains_dictionary.pkl'
        
        self.n_truncations = []
        self.c_truncations = []
        self.chain_key_res = []
        self.looped_pdbs = []
        self.output_pdbs = []

        #For job queue, to parallel jobs
        self.queues = []

        #for loop struct infos
        self.infos = []
        self.combs = []


    def build_protein(self):
        """Iteratively generate a protein using MASTER and Qbits."""
        self.log.info('Start build protein.')
        self._generate_proteins(self.para.num_iter)
        print('output pdbs :')
        print('\n'.join(self.output_pdbs))
        self.log.info('Finish build protein.')

    def build_protein_parallel(self):
        """Iteratively generate a protein using MASTER and Qbits."""
        self._generate_proteins_parallel(self.para.num_iter)
        print('output pdbs :')
        print('\n'.join(self.output_pdbs))

    def build_protein_deprecate(self):
        """Iteratively generate a protein using MASTER and Qbits."""
        self._generate_recursive(self.para.num_iter)
        print('output pdbs :')
        print('\n'.join(self.output_pdbs))
 
    def loop_seed_structure(self, n_truncations=[], c_truncations=[], chain_key_res=[]):
        """Treat the seed PDB as a complete structure and generate loops.

        Parameters
        ----------
        n_truncations : list
           List of integer counts of residues to truncate from the N-terminus 
           of each chain in the seed structure prior to loop generation.
        c_truncations : list
           List of integer counts of residues to truncate from the C-terminus 
           of each chain in the seed structure prior to loop generation.
        chain_key_res : list
           List of lists, one for each chain in the seed structure, denoting 
           the indices (counted up from 0 at the N-terminus of each chain, 
           regardless of the residue indices in the PDB file) of residues to 
           be retained in steric clash calculations.
        """
        # ensure a seed structure was provided to the class constructor
        if len(self.full_sse_list) == 0:
            raise AssertionError('seed_pdb not provided to constructor.')
        if len(self.pdbs) > 1:
            raise AssertionError('build_protein() has already been run.')
        # compute the number of satisfied N- and C-termini
        sat = pdbutils.satisfied_termini(self.pdbs[0], self.para.max_nc_dist)
        # set n_truncations and c_truncations for loop generation
        # self.n_truncations = n_truncations
        # self.c_truncations = c_truncations
        self.chain_key_res = chain_key_res
        # generate loops
        _full_sse_list = self.full_sse_list.copy()
        self._generate_loops(_full_sse_list, sat, self.workdir, self.loop_range) 
        print('output pdbs :')
        print('\n'.join(self.output_pdbs))

    def loop_seed_single_structure(self, direction=[], n_truncations=[0], c_truncations=[0], chain_key_res=[]):
        if len(self.full_sse_list) == 0:
            raise AssertionError('seed_pdb not provided to constructor.')
        if len(self.pdbs) > 1:
            raise AssertionError('build_protein() has already been run.')
        # compute the number of satisfied N- and C-termini
        all_sat = pdbutils.satisfied_termini(self.pdbs[0], self.para.max_nc_dist)
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
        # self.n_truncations = n_truncations
        # self.c_truncations = c_truncations
        self.chain_key_res = chain_key_res
        # generate loops
        _full_sse_list = self.full_sse_list.copy()
        for sse in _full_sse_list:
            self.log.info(sse)
        #self._generate_loops(_full_sse_list, sat, self.workdir, self.loop_range)   
        self._generate_trunc_loops(direction, sat, self.workdir, n_truncations, c_truncations, self.para.cluster_count_cut, self.loop_range)     
        print('output pdbs :')
        print('\n'.join(self.output_pdbs))

    ### NEW FUNCTIONS FOR GENERATING LOOPS
    
    def _cal_win_dist(self, full_sse_list):
        qrep_xyz = []
        qrep_natoms = [0]
        qrep_nres = []
        backbone = ['N', 'CA', 'C', 'O']
        for j, pair_qrep in enumerate(full_sse_list):
            title = 'struct{}'.format(str(j))
            this_struct = pdbutils.get_struct(title, pair_qrep, self.para.min_nbrs)
            atoms = this_struct.get_atoms()
            qrep_xyz.append(np.array([atom.get_coord() for atom in atoms if atom.get_name() in backbone]))            
            qrep_nres.append(int(len(qrep_xyz[-1])/4))
            qrep_natoms.append(qrep_natoms[-1] + len(qrep_xyz[-1]))

        qrep_xyz = np.vstack(qrep_xyz)
        # compute interatomic distance between SSE pairs
        dists = cdist(qrep_xyz, qrep_xyz)

        n_reps = len(full_sse_list)

        #calcualte distances (with window = 7) between two the first sse and another sse, 
        #store the minimum index.
        win_dists = np.zeros((4, qrep_nres[0]- 6))
        win_dists[0] = np.arange(qrep_nres[0]- 6)
        for k in range(1, n_reps):
            for x in range(qrep_nres[0]- 6):
                win_dist = [0]* (qrep_nres[k]- 6)
                for y in range(qrep_nres[k]- 6):
                    for z in range(28):
                        win_dist[y] += dists[qrep_natoms[0]+x*4 + z, qrep_natoms[k]+y*4 + z]
                win_dists[k, x]=win_dist.index(min(win_dist))
        return qrep_nres, win_dists  

    # It may be better to use z index for helix bundle.
    def _get_truncs(self, full_sse_list, sat, n_truncations=[], c_truncations=[]):
        qrep_nres, win_dists = self._cal_win_dist(full_sse_list)             

        n_chains = len(full_sse_list)
        all_n_truncs = []
        all_c_truncs = []
        all_local_sats = []
        all_directs = []

        direction = list(range(n_chains))

        for n in n_truncations:
            ns = [0]*n_chains
            cs = [0]*n_chains
            the_sat = sat.copy()
            for i in range(n_chains):
                if i%2==0:
                    the_sat[i, ] = 0
                    ns[direction[i]] = win_dists[direction[i]][n]
                else:
                    cs[direction[i]] = qrep_nres[direction[i]] - 7 - win_dists[direction[i]][n]
            all_n_truncs.append(ns)
            all_c_truncs.append(cs)
            all_local_sats.append(the_sat)
            all_directs.append('n')

        for c in c_truncations:
            cn = qrep_nres[0] - c -7
            ns = [0]*n_chains
            cs = [0]*n_chains
            the_sat = sat.copy()
            for i in range(n_chains):
                if i%2==0:
                    cs[direction[i]] = qrep_nres[direction[i]] - 7 - win_dists[direction[i]][cn]
                else:
                    the_sat[i, ] = 0
                    ns[direction[i]] = win_dists[direction[i]][cn]
            all_n_truncs.append(ns)
            all_c_truncs.append(cs)
            all_local_sats.append(the_sat)
            all_directs.append('c')
        return all_n_truncs, all_c_truncs, all_local_sats, all_directs

    def _get_truncs_depre(self, direction, n_truncations, c_truncations):
        all_n_truncs = []
        all_c_truncs = []
        #assume the sse are ordered in alpha beta seq
        if len(direction)==0:
            direction = list(range(len(direction)+1))

        for n in n_truncations:
            ns = [0]*len(direction)
            cs = [0]*len(direction)
            for i in range(len(direction)):
                if i%2==0:
                    ns[direction[i]] = n
                else:
                    cs[direction[i]] = n
            all_n_truncs.append(ns)
            all_c_truncs.append(cs)
        for c in c_truncations:
            ns = [0]*len(direction)
            cs = [0]*len(direction)
            for i in range(len(direction)):
                if i%2==0:
                    cs[direction[i]] = c
                else:
                    ns[direction[i]] = c
            all_n_truncs.append(ns)
            all_c_truncs.append(cs)
        return all_n_truncs, all_c_truncs

    def _generate_trunc_loops(self, direction, sat, workdir, n_truncations=[], c_truncations=[], cluster_count_cut=20, loop_range=[3, 20]):
        n_chains = len(sat)
        # find loops for each pair of nearby N- and C-termini
        all_n_truncs, all_c_truncs, all_local_sats, all_directs = self._get_truncs(self.full_sse_list, sat, n_truncations, c_truncations)
        # for each trunc, find all loops.
        for i in range(len(all_n_truncs)):
            ns = all_n_truncs[i]
            cs = all_c_truncs[i]
            the_sat = all_local_sats[i]
            the_direct = all_directs[i]
            trunc_len = int(ns[0]) if the_direct == 'n' else int(cs[0])
            _workdir = workdir +'/'+ the_direct +'_trunc_{}'.format(str(trunc_len))
            if not os.path.exists(_workdir):
                os.mkdir(_workdir)           
            pdbutils.split_pdb(self.pdbs[-1], _workdir, self.para.min_nbrs, None, ns, cs)
            _the_full_sse_list = [_workdir + '/' + path for path in os.listdir(_workdir) if 'chain_' in path]
            #Search loops
            slice_lengths = self._loop_search_fast(_the_full_sse_list, the_sat, _workdir, loop_range)
            #loop_success = self._get_loop_success(the_sat, _workdir, loop_range)      
            #Summarize loops
            _infos = self._get_top_cluster_summary(_workdir, n_chains, slice_lengths, cluster_count_cut, loop_range)
            self.infos.extend(_infos)
        #Write struct info
        self._write_file(workdir + '/summary.txt', self.infos)  
        #Find loop combine candidates.
        combs, ps, scores = self._extract_top_hit(n_chains, sat, cluster_count_cut)
        self.combs.extend(combs)
        # #Build whole structure.
        outdir = workdir + '/output'
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        
        output_cut = 30 if len(combs) > 30 else len(combs)
        inds = np.argsort(scores)[::-1][:output_cut]
        for i in inds:
            comb = combs[i]
            p = ps[i]
            centroids = [c.cent_pdb for c in comb]
            clashing = pdbutils.check_clashes(centroids)
            if clashing:
                continue
            structs, slices = self._connect_loops_struct(self.full_sse_list, n_chains, p, centroids)
            out_path = outdir + '/output_' + '-'.join(str(v) for v in p) + '_' + str(i) + '.pdb'
            pdbutils.merge_save_struct(out_path, structs, slices)      

    def _get_top_cluster_summary(self, workdir, n_chains, slice_lengths, cluster_count_cut, loop_range):
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
                    sl = slice_lengths[p[0], p[1]]
                    cluster_key_res = self._key_residues(loop_workdir + '/seq.txt', loop_pdbs, sl)                
                    #Copy centroid pdb
                    _cent_pdb = ''
                    _cent_pdb_workdir = self.para.workdir + '/loops_{}_{}'.format(string.ascii_uppercase[p[0]], string.ascii_uppercase[p[1]])
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
                    
                        _cent_pdb_log = _cent + '.png'
                        self._plot_log(loop_workdir + '/seq.txt', l, _cent_pdb_log)               

                        _cent_pdb_phipsi = _cent + '_phipsi.png'
                        phipsi = pdbutils.meaure_phipsi(_cent_pdb)
                        self._plot_phipsi(phipsi, _cent_pdb_phipsi)    
                        
                    #Add summary
                    info = Struct_info(trunc_info = loop_workdir.split('/')[-2], loop_info = loop_workdir.split('/')[-1], 
                            loop_len = l, clust_num = len(loop_pdbs), cent_pdb = _cent_pdb, cluster_key_res = cluster_key_res)              
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

    def _extract_top_hit(self, n_chains, sat, cluster_count_cut):
        combs = []
        ps = []
        for p in permutations(range(n_chains)):
            if not np.all([sat[p[j], p[j+1]] for j in range(n_chains - 1)]):
                continue
            all_keys = []
            for j in range(n_chains-1):
                k = j+1
                loop_key = 'loops_{}_{}'.format(string.ascii_uppercase[p[j]], string.ascii_uppercase[p[k]])
                keys = [key for key in self.infos if loop_key in key.loop_info]
                values = [v.clust_num for v in keys]
                all_keys.append([keys[i] for i in np.argsort(values) if values[i] > cluster_count_cut])
            if 0 in [len(v) for v in all_keys]:
                continue
            for comb in product(*all_keys):
                combs.append(comb)
                ps.append(p)
        scores = []
        for comb in combs:
            scores.append(sum([c.clust_num for c in comb]))
        return combs, ps, scores

    def _connect_loops_struct(self, _full_sse_list, n_chains, permutation, centroids):
        backbone = ['N', 'CA', 'C', 'O']
        pdbs_to_combine = [''] * (2 * n_chains - 1)
        pdbs_to_combine[::2] = [_full_sse_list[idx] for idx in permutation]
        pdbs_to_combine[1::2] = centroids

        inds = [[0,0] for i in range(len(pdbs_to_combine))]
        for i in range(len(pdbs_to_combine)-1):
            j = i + 1
            if i%2==0:
                qrep_xyz = []
                title2 = 'struct{}'.format(str(j))
                struct2 = pdbutils.get_struct(title2, pdbs_to_combine[j], self.para.min_nbrs)
                pos2 = [atom.get_coord() for atom in list(struct2.get_residues())[0].get_atoms() if atom.get_name() in backbone]   
                qrep_xyz.append(np.array(pos2))

                title1 = 'struct{}'.format(str(i))
                struct1 = pdbutils.get_struct(title1, pdbs_to_combine[i], self.para.min_nbrs)
                atoms1 = struct1.get_atoms()        
                qrep_xyz.append(np.array([atom.get_coord() for atom in atoms1 if atom.get_name() in backbone])) 
                
                inds[i][1] = len(list(struct1.get_residues()))
                inds[j][1] = len(list(struct2.get_residues()))
                qrep_xyz = np.vstack(qrep_xyz)
                # compute interatomic distance between SSE pairs
                dists = cdist(qrep_xyz, qrep_xyz)

                #store the minimum index.
                pdb1_len = len(list(struct1.get_residues()))
                win_dist = [0]*pdb1_len
                for x in range(0, pdb1_len):                  
                    for z in range(4):
                        win_dist[x] += dists[z, x*4 + z]
                ind=win_dist.index(min(win_dist[1:]))-1               
                inds[i][1]=ind
            else:
                qrep_xyz = []
                title2 = 'struct{}'.format(str(i))
                struct2 = pdbutils.get_struct(title2, pdbs_to_combine[i], self.para.min_nbrs)
                pos2 = [atom.get_coord() for atom in list(struct2.get_residues())[-1].get_atoms() if atom.get_name() in backbone]   
                qrep_xyz.append(np.array(pos2))

                title1 = 'struct{}'.format(str(j))
                struct1 = pdbutils.get_struct(title1, pdbs_to_combine[j], self.para.min_nbrs)
                atoms1 = struct1.get_atoms()        
                qrep_xyz.append(np.array([atom.get_coord() for atom in atoms1 if atom.get_name() in backbone])) 

                inds[j][1] = len(list(struct1.get_residues()))
                inds[i][1] = len(list(struct2.get_residues()))
                qrep_xyz = np.vstack(qrep_xyz)
                # compute interatomic distance between SSE pairs
                dists = cdist(qrep_xyz, qrep_xyz)

                #store the minimum index.
                pdb1_len = len(list(struct1.get_residues()))
                win_dist = [0]*pdb1_len
                for x in range(0, pdb1_len):                  
                    for z in range(4):
                        win_dist[x] += dists[z, x*4 + z]
                ind=win_dist.index(min(win_dist[1:]))
                inds[j][0] = ind
                               
        slices = [] 
        for d in inds:
            slices.append(slice(d[0], d[1]))
        structs = [pdbutils.get_struct('test_' + str(i), pdbs_to_combine[i]) for i in range(len(pdbs_to_combine))]
        return structs, slices

    #check_clash for one combination of chains in top cluster. 
    def _check_single_clash_depre(self, chain_pdbs, permutation, loop_comb):
        centroids = []
        res_ids_to_keep = []
        for centroid in loop_comb:         
            # check if loop clashes with exclusion PDB
            # or SSEs it does not connect
            loop_clashes = False
            if self.orig_exclusion is not None:
                loop_clashes = loop_clashes or \
                    pdbutils.check_clashes([centroid.cent_pdb, self.orig_exclusion])
            unlooped_sses = [chain_pdbs[sse_idx] for sse_idx in permutation 
                                if sse_idx not in permutation[j:j+2]]
            if len(self.chain_key_res) == len(unlooped_sses):
                sse_key_res = [self.chain_key_res[p_idx] for p_idx in permutation]
            else:
                sse_key_res = [[]] * len(unlooped_sses)
            loop_clashes = loop_clashes or \
                pdbutils.check_clashes([centroid.cent_pdb] + unlooped_sses, [centroid.cluster_key_res] + sse_key_res)
            if loop_clashes or pdbutils.check_gaps(centroid.cent_pdb):
                break
            centroids.append(centroid.cent_pdb)
            res_ids_to_keep.append(centroid.cluster_key_res)
        return centroids, res_ids_to_keep
      
    def _write_file(self, filename, infos):
        with open(filename, 'w') as f:
            f.write('trunc_info\tloop_info\tloop_len\tclust_num\tcent_pdb\tcluster_key_res\tclust_num_2nd\n')
            for r in infos:
                f.write(r.trunc_info + '\t' + r.loop_info + '\t' + str(r.loop_len) + '\t'
                    + str(r.clust_num) + '\t' + r.cent_pdb + '\t' + str(r.clust_num_2nd) + '\n')       
    
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
                         font_name='Stencil Std',
                         color_scheme='NajafabadiEtAl2017',
                         vpad=.1,
                         width=.8)
        logo.style_xticks(anchor=0, spacing=1)      
        logo.ax.set_ylabel('Count')
        logo.ax.set_xlim([-1, len(df)])
        #logo.fig.savefig(filepath) 
        plt.savefig(filepath)
        plt.close()

    def _plot_phipsi(self, phipsi, filepath):
        x = list(range(1, len(phipsi)+1))
        plt.plot(x, phipsi)
        plt.savefig(filepath)
        plt.close()

    ### FUNCTIONS FOR GENERATING LOOPS

    def _generate_loops(self, _full_sse_list, sat, workdir, loop_range=[3, 20]):
        n_chains = len(sat)
        # find loops for each pair of nearby N- and C-termini
        slice_lengths = self._loop_search_fast(_full_sse_list, sat, workdir, loop_range)
        loop_success = self._get_loop_success(sat, workdir, loop_range)
        outfiles = []
        counter = 0
        for p in permutations(range(n_chains)):
            # if loops were built between all successive SSEs in the 
            # permutation of SSE order, continue on to add in the loops
            if np.all([loop_success[p[j], p[j+1]] for j in range(n_chains - 1)]):
                all_centroids, num_clusters, cluster_key_res, no_clusters = self._get_top_clusters(workdir, n_chains, slice_lengths, p, loop_range)
                if no_clusters:
                    continue
                # test whether any selection of loops avoids clashing
                some_outfiles, counter = \
                    self._test_topologies(_full_sse_list, workdir, p, all_centroids, 
                                          cluster_key_res, slice_lengths, 
                                          num_clusters, n_chains, counter)
                outfiles += some_outfiles
        self.output_pdbs += outfiles
    
    def _loop_search(self, _full_sse_list, sat, workdir, loop_range=[3, 20]):
        n_chains = len(sat)
        loop_success = np.zeros_like(sat)
        slice_lengths = np.zeros((sat.shape[0], sat.shape[1], 2), dtype=int)
        for j, k in product(range(n_chains), repeat=2):
            # ensure selected SSEs satisfy the distance constraint
            if not sat[j, k] or j == k:
                continue
            print('Generating loops between SSEs {} '
                  'and {}'.format(string.ascii_uppercase[j], 
                                  string.ascii_uppercase[k]))
            loop_workdir = workdir + \
                           '/loops_{}_{}'.format(string.ascii_uppercase[j], 
                                                 string.ascii_uppercase[k])
            if not os.path.exists(loop_workdir):
                os.mkdir(loop_workdir)
            loop_query = loop_workdir + '/loop_query.pdb'
            loop_outfile = loop_workdir + '/stdout'
            # calculate how many residues are required for an overlap region 
            # of length 10 Angstroms between the query SSEs and the loops
            slice_lengths[j, k] = \
                pdbutils.gen_loop_query([_full_sse_list[j], 
                                         _full_sse_list[k]], 
                                        loop_query, min_nbrs=self.para.min_nbrs)
            # find loops with MASTER
            clusters_exist = True
            for l in range(loop_range[0], loop_range[1] + 1):
                subdir = loop_workdir + '/{}'.format(str(l))
                if not os.path.exists(subdir):
                    os.mkdir(subdir)
                    print('Querying MASTER for loops '
                          'of length {}'.format(str(l)))
                    query.master_query_loop(loop_query, self.para.loop_target_list, 
                                            rmsdCut=self.para.rmsdCut, topN=200,
                                            gapLen=l, outdir=subdir, 
                                            outfile=loop_outfile)
                if not os.path.exists(subdir + '/clusters'):
                    clusters_exist = False
            # cluster loops if the clusters do not already exist
            if not clusters_exist:
                cluster_loops.run_cluster(loop_workdir + '/', 
                                          outfile=loop_outfile)
            # determine whether clustering succeeded for any loop length
            for l in range(loop_range[0], loop_range[1] + 1):
                subdir = loop_workdir + '/{}'.format(str(l))
                if len(os.listdir(subdir + '/clusters')) > 0:
                    loop_success[j, k] = 1
        return loop_success, slice_lengths

    #Find loops for each pair of nearby N- and C-termini. return slice_lengths?
    def _loop_search_fast(self, _full_sse_list, sat, workdir, loop_range=[3, 20]):
        n_chains = len(sat)
        slice_lengths = np.zeros((sat.shape[0], sat.shape[1], 2), dtype=int)
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
            slice_lengths[j, k] = pdbutils.gen_loop_query([_full_sse_list[j], 
                                                          _full_sse_list[k]], 
                                                          loop_query, min_nbrs=self.para.min_nbrs)
            # find loops with MASTER
            # sort PDBs into directories by loop length
            clusters_exist = self._loop_search_query_search(loop_workdir, loop_query, loop_outfile, loop_range)
            # cluster loops if the clusters do not already exist
            if not clusters_exist:
                cluster_loops.run_cluster(loop_workdir + '/', outfile=loop_outfile)
        return slice_lengths

    def _loop_search_query_search(self, loop_workdir, loop_query, loop_outfile, loop_range):
        # find loops with MASTER
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

    #The function is to find if there is loop found between the two query chains. 
    def _get_loop_success(self, sat, workdir, loop_range=[3, 20]):
        n_chains = len(sat)
        loop_success = np.zeros_like(sat)
        # determine whether clustering succeeded for any loop length
        for j, k in product(range(n_chains), repeat=2):
            loop_workdir = workdir + '/loops_{}_{}'.format(string.ascii_uppercase[j], string.ascii_uppercase[k])
            for l in range(loop_range[0], loop_range[1] + 1):
                subdir = loop_workdir + '/{}'.format(str(l))
                if os.path.exists(subdir):
                    if len(os.listdir(subdir + '/clusters')) > 0:
                        loop_success[j, k] = 1
        return loop_success

    def _get_top_clusters(self, workdir, n_chains, slice_lengths, p, loop_range=[3, 20]):
        all_centroids = []
        num_clusters = []
        cluster_key_res = []
        no_clusters = False
        if not os.path.exists(workdir + '/loop_centroids'):
            os.mkdir(workdir + '/loop_centroids')
        for j in range(n_chains - 1):          
            loop_workdir = workdir + '/loops_{}_{}'.format(string.ascii_uppercase[p[j]], string.ascii_uppercase[p[j+1]])
            l_centroids = []
            cluster_sizes = []
            l_cluster_key_res = []
            for l in range(loop_range[0], loop_range[1] + 1):
                subdir = loop_workdir + '/{}/clusters/1'.format(str(l))                 
                try:
                    loop_pdbs = os.listdir(subdir)
                except:
                    loop_pdbs = []
                if len(loop_pdbs) > 0:
                    if self.para.lowest_rmsd_loop:
                        l_centroids.append(subdir + '/' + self._lowest_rmsd_loop(loop_workdir + '/match.txt', loop_pdbs))
                    else:
                        l_centroids.append([subdir + '/' + lpdb for lpdb in loop_pdbs if 'centroid' in lpdb][0])                   
                    cluster_sizes.append(len(loop_pdbs))                                                       
                    sl = slice_lengths[p[j], p[j+1]]
                    # find "key positions" along the loop that are 
                    # statistically enriched in one residue, so their  
                    # side chains can be included in clash checks
                    l_cluster_key_res.append(self._key_residues(loop_workdir + '/seq.txt', loop_pdbs, sl))                     
            # sort clusters by size
            if len(cluster_sizes) > 0:
                idxsort = np.argsort(cluster_sizes)[::-1]
                all_centroids.append([l_centroids[idx] for idx in idxsort])
                num_clusters.append(len(cluster_sizes))
                cluster_key_res.append([l_cluster_key_res[idx] for idx in idxsort])
            else:
                no_clusters = True
        return all_centroids, num_clusters, cluster_key_res, no_clusters, summaries
                   
    def _lowest_rmsd_loop(self, matchfile, loop_pdbs):
        with open(matchfile, 'r') as f:
            rmsds = [float([val for val in line.split(' ') if val != ''][0]) 
                     for line in f.read().split('\n') if len(line) > 0]
        loop_rmsds = []
        for loop_pdb in loop_pdbs:
            if 'wgap' in loop_pdb:
                idx = int(loop_pdb.split('_')[-1][4:-7]) - 1
            else:
                idx = int(loop_pdb.split('_')[-1][5:-7]) - 1
            loop_rmsds.append(rmsds[idx])
        return loop_pdbs[np.argmin(loop_rmsds)]

    def _key_residues(self, seqfile, loop_pdbs, slice_lengths):
        with open(seqfile, 'r') as f:
            lines = f.read().split('\n')
        all_seqs = []
        for line in lines:
            if len(line) > 0:
                reslist = []
                for res in line.split(' '):
                    if len(res) == 3 and res[0].isalpha():
                        reslist.append(res)
                    elif len(res) == 4 and res[0] == '[':
                        reslist.append(res[1:])
                    elif len(res) == 4 and res[-1] == ']':
                        reslist.append(res[:-1])
                all_seqs.append(reslist)
        seqs = []
        for loop_pdb in loop_pdbs:
            if 'wgap' in loop_pdb:
                idx = int(loop_pdb.split('_')[-1][4:-7]) - 1
            else:
                idx = int(loop_pdb.split('_')[-1][5:-7]) - 1
            seqs.append(all_seqs[idx])
        seqs = np.array(seqs, dtype=str)
        try:
            modes = mode(seqs, axis=0)
        except:
            return []
        # return which positions have one residue across most cluster members
        idxs = np.argwhere(modes.count > 0.7 * len(seqs)).reshape(-1)
        idx_min = slice_lengths[0]
        idx_max = seqs.shape[1] - slice_lengths[1]
        return [idx for idx in idxs if idx > idx_min and idx < idx_max]

    def _test_topologies(self, _full_sse_list, workdir, permutation, all_centroids, 
                         cluster_key_res, slice_lengths, num_clusters, 
                         n_chains, counter):
        some_outfiles = []
        # iterate through loop structures until one is found without 
        # clashes and (optionally) satisfying a compactness criterion
        pdb_dir = os.path.dirname(self.pdbs[-1])
        pdbutils.split_pdb(self.pdbs[-1], pdb_dir, self.para.min_nbrs, None, 
                           self.n_truncations, self.c_truncations)
        chain_pdbs = [pdb_dir + '/' + path for path 
                      in os.listdir(pdb_dir) if 'chain_' in path]
        chain_pdbs.sort()
        #get all combinations of matches in top clusters from each chain.
        loop_idx_sets = np.array(list(product(*[range(num) for num in num_clusters])))
        loop_idx_sets = loop_idx_sets[np.argsort(loop_idx_sets.sum(axis=1))]
        #forbidden is a list of list of matches for each chain that find clash. No need to check clash of it again.
        forbidden = [[]] * len(all_centroids)
        for idxs in loop_idx_sets:
            skip_idxs = False
            for j in range(len(all_centroids)):
                if idxs[j] in forbidden[j]:
                    skip_idxs = True
            if skip_idxs:
                continue
            centroids_gz = [all_centroids[j][idxs[j]] for j in range(len(all_centroids))]
            filenames, centroids, res_ids_to_keep = self._check_clash(workdir, centroids_gz, forbidden, chain_pdbs, permutation, cluster_key_res, idxs)           
            # ensure enough loops were found and check clashes between them
            if len(centroids) != n_chains - 1 or pdbutils.check_clashes(centroids, res_ids_to_keep):
                [os.remove(filename) for filename in filenames]
                continue
            # permute SSEs and connect with loops
            clashing, outfile_path = self._connect_loops(_full_sse_list, workdir, n_chains, permutation, filenames, centroids, counter, slice_lengths)
            if clashing:
                [os.remove(filename) for filename in filenames]
                continue
            if self.para.screen_compactness:
                compactness = pdbutils.calc_compactness(outfile_path)
                if compactness < 0.138:
                    [os.remove(filename) for filename in filenames]
                    continue
            print('full protein output :', outfile_path)
            #print('pdbs_to_combine :', pdbs_to_combine)
            some_outfiles.append(outfile_path)
            counter += 1
            break
        return some_outfiles, counter

    #check_clash for one combination of chains in top cluster. 
    def _check_clash(self, workdir, centroids_gz, forbidden, chain_pdbs, permutation, cluster_key_res, idxs):
        filenames = []
        centroids = []
        res_ids_to_keep = []
        for j, centroid in enumerate(centroids_gz):
            dest_name = workdir + '/loop_centroids/' + os.path.basename(centroid)[:-3]
            if not os.path.exists(dest_name):
                with gzip.open(centroid, 'rb') as f_in:
                    with open(dest_name, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                filenames.append(dest_name)
                # check if loop clashes with exclusion PDB
                # or SSEs it does not connect
                loop_clashes = False
                if self.orig_exclusion is not None:
                    loop_clashes = loop_clashes or \
                        pdbutils.check_clashes([dest_name, self.orig_exclusion])
                unlooped_sses = [chain_pdbs[sse_idx] for sse_idx in permutation 
                                    if sse_idx not in permutation[j:j+2]]
                if len(self.chain_key_res) == len(unlooped_sses):
                    sse_key_res = [self.chain_key_res[p_idx] for p_idx in permutation]
                else:
                    sse_key_res = [[]] * len(unlooped_sses)
                loop_clashes = loop_clashes or \
                    pdbutils.check_clashes([dest_name] + unlooped_sses, [cluster_key_res[j][idxs[j]]] + sse_key_res)
                if loop_clashes or pdbutils.check_gaps(dest_name):
                    # remove loops that clash with SSEs or that have gaps 
                    # from future consideration
                    forbidden[j].append(idxs[j])
                    break
            centroids.append(dest_name)
            res_ids_to_keep.append(cluster_key_res[j][idxs[j]])
        return filenames, centroids, res_ids_to_keep
        
    # permute SSEs and connect with loops      
    def _connect_loops(self, _full_sse_list, workdir, n_chains, permutation, filenames, centroids, counter, slice_lengths):
        pdbs_to_combine = [''] * (2 * n_chains - 1)
        pdbs_to_combine[::2] = [_full_sse_list[idx] for idx in permutation]
        pdbs_to_combine[1::2] = centroids
        outfile_path = workdir + '/output_{}.pdb'.format(str(counter))
        filenames.append(outfile_path)
        overlaps = []
        for j in range(n_chains - 1):
            overlaps.append(slice_lengths[permutation[j], permutation[j+1], 0])
            overlaps.append(slice_lengths[permutation[j], permutation[j+1], 1])
        clashing = pdbutils.stitch(pdbs_to_combine, outfile_path, 
                                    overlaps=overlaps, 
                                    min_nbrs=self.para.min_nbrs, 
                                    from_closest=False)
        return clashing, outfile_path
    
    ### FUNCTIONS FOR GENERATING SSEs

    def _generate_recursive(self, recursion_order):
        print('Adding a qbit rep.')
        # search for a contiguous secondary structural element to add
        outdir = os.path.dirname(self.pdbs[-1])
        outfile = outdir + '/stdout'
        print('Querying MASTER')
        query.master_query(self.pdbs[-1], self.targetList, self.para.rmsdCut, 
            topN=None, outfile=outfile, clobber=False)
        print('Searching with Qbits')
        if not os.path.exists(outdir + '/qbit_reps/'):
            try:
                # ensure the second SSE is antiparallel to the first
                query_exists = int(bool(len(self.query_sse_list)))
                first_recursion = (recursion_order == self.para.num_iter - query_exists)
                query.qbits_search(self.pdbs[-1], self.exclusion_pdbs[-1], 
                                   self.chains_dict, outdir, 
                                   self.para.qbits_window, self.para.qbits_rmsd, 
                                   top=10, sec_struct=self.para.secstruct,
                                   antiparallel=first_recursion,
                                   min_nbrs=self.para.min_nbrs, contiguous=True)
            except:
                pass
        if os.path.exists(outdir + '/qbit_reps/'):
            qreps = [outdir + '/qbit_reps/' + pdb_path for pdb_path in os.listdir(outdir + '/qbit_reps/')]
            # iterate over qbit reps to attempt adding more structure to each 
            # until a suitable small protein is found
            print('Testing Qbit reps')
            for i, qrep in enumerate(qreps):
                self._add_qbit_rep(qrep, i, recursion_order, outdir)

    def _add_qbit_rep(self, qrep, i, recursion_order, outdir):
        # determine permissible sets of seed SSEs for Qbits searches
        seed_sse_lists = self._prepare_seed_sse_lists(qrep)
        # iterate over permissible sets of seed SSEs to add more Qbit reps
        for j, seed_sse_list in enumerate(seed_sse_lists):
            _workdir = '{}/{}'.format(outdir, str(i) + string.ascii_lowercase[j])
            if not os.path.exists(_workdir):
                os.mkdir(_workdir)
            _seed_pdb = _workdir + '/seed.pdb'
            pdbutils.merge_pdbs(seed_sse_list, _seed_pdb, min_nbrs=self.para.min_nbrs)
            self.full_sse_list.append(qrep)
            print('SSE List:')
            print('\n'.join(self.full_sse_list))
            # compute the number of satisfied N- and C-termini
            if len(self.full_sse_list) > 2:
                _full_pdb = _workdir + '/full.pdb'
                pdbutils.merge_pdbs(self.full_sse_list, _full_pdb, min_nbrs=self.para.min_nbrs)
                sat = pdbutils.satisfied_termini(_full_pdb, self.para.max_nc_dist)
                self.pdbs.append(_full_pdb)
            else:
                sat = pdbutils.satisfied_termini(_seed_pdb, self.para.max_nc_dist)
                self.pdbs.append(_seed_pdb)
            n_sat = np.sum(sat)
            if self.para.num_iter - recursion_order - 1 > n_sat:
                # if it is impossible to satisfy all N- or C- termini within 
                # the remaining number of iterations, exit the branch early
                self.pdbs = self.pdbs[:-1]
                self.full_sse_list = self.full_sse_list[:-1]
                continue
            # if recursion_order is not 1, continue adding qbit reps 
            if recursion_order > 1:
                _exclusion_pdb = _workdir + '/exclusion.pdb'
                pdbutils.merge_pdbs([self.exclusion_pdbs[-1], qrep], _exclusion_pdb, min_nbrs=self.para.min_nbrs)
                self.exclusion_pdbs.append(_exclusion_pdb)
                self._generate_recursive(recursion_order - 1)
            # if recursion order is 1 and there are enough N/C termini 
            # satisfied, try building loops
            elif n_sat >= self.para.num_iter:
                try_loopgen = False
                n_chains = len(sat)
                for p in permutations(range(n_chains)):
                    if np.all([sat[p[k], p[k+1]] for k in range(n_chains - 1)]):
                        try_loopgen = True
                # check to make sure self.pdbs[-1] hasn't been looped before
                if try_loopgen:
                    with open(self.pdbs[-1], 'r') as f0:
                        f0_read = f0.read()
                        for pdb in self.looped_pdbs:
                            with open(pdb, 'r') as f1:
                                if f0_read == f1.read():
                                    try_loopgen = False
                # if necessary, ensure the compactness criterion is met
                if self.para.screen_compactness:
                    compactness = pdbutils.calc_compactness(self.pdbs[-1])
                    try_loopgen = try_loopgen and (compactness > 0.1)
                if try_loopgen:
                    self.looped_pdbs.append(self.pdbs[-1])
                    self._generate_loops(sat, _workdir, self.loop_range)
            # if unsuccessful, remove the PDB from the running lists
            self.pdbs = self.pdbs[:-1]
            if len(self.exclusion_pdbs) == len(self.pdbs) + 1:
                self.exclusion_pdbs = self.exclusion_pdbs[:-1]
            self.full_sse_list = self.full_sse_list[:-1]
            # do not iterate over other SSE lists if recursion is complete
            if recursion_order == 1:
                break

    def _prepare_seed_sse_lists(self, qrep):
        # if using an external query structure, include it in the Qbits 
        # search until the protein under construction has at least two SSEs
        if len(self.full_sse_list) < 2:
            all_reps = self.query_sse_list + self.full_sse_list + [qrep]
        else:
            all_reps = self.full_sse_list + [qrep]
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

    ### NEW FUNCTIONS FOR GENERATING SSEs

    def _generate_proteins(self, recursion_order):
        self.queues.append([self.pdbs[-1], self.exclusion_pdbs[-1], self.full_sse_list, recursion_order])
        while len(self.queues) > 0:            
            queue = self.queues.pop(0)
            print('--Get queue--')
            self._const_protein(queue[0], queue[1], queue[2], queue[3])             
        print('---out while loop---')

    def _const_protein(self, pdb, exclusion_pdb, full_sse_list, recursion_order):
        outdir = os.path.dirname(pdb)
        #Construct final protein.
        if recursion_order == 0:
            self._const_prot_loop(pdb, exclusion_pdb, full_sse_list, recursion_order, outdir)    
            return
        #Generate queries for next recursion.
        qreps = self._generate_qreps(pdb, exclusion_pdb, recursion_order, outdir)
        if qreps == None:
            return
        for i, qrep in enumerate(qreps):
            seed_sse_lists = self._prepare_seed_sses(qrep, full_sse_list)
            for j, seed_sse_list in enumerate(seed_sse_lists):
                self._add_seed_sse(i, qrep, j, seed_sse_list, full_sse_list, exclusion_pdb, recursion_order, outdir)
                #_queues = self._add_seed_sse(i, qrep, j, seed_sse_list, full_sse_list, exclusion_pdb, recursion_order, outdir)
        # for _item in _queues:
        #     self.queues.append(_item)   
        #print('---add more into queue---')
        return

    #The function will be called when the recursion_order == 0 to construct the final small protein.
    def _const_prot_loop(self, pdb, exclusion_pdb, full_sse_list, recursion_order, outdir):
        sat = pdbutils.satisfied_termini(pdb, self.para.max_nc_dist)
        n_sat = np.sum(sat)
        if self.para.num_iter - recursion_order - 2 > n_sat:
            # if it is impossible to satisfy all N- or C- termini within 
            # the remaining number of iterations, exit the branch early
            return
        if n_sat >= self.para.num_iter:
            try_loopgen = False
            n_chains = len(sat)
            for p in permutations(range(n_chains)):
                if np.all([sat[p[k], p[k+1]] for k in range(n_chains - 1)]):
                    try_loopgen = True
            # check to make sure self.pdbs[-1] hasn't been looped before
            if try_loopgen:
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
                self._generate_loops(full_sse_list, sat, outdir, self.loop_range)

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
                query_exists = int(bool(len(self.query_sse_list)))
                first_recursion = (recursion_order == self.para.num_iter - query_exists)
                query.qbits_search(pdb, exclusion_pdb, 
                                   self.chains_dict, outdir, 
                                   self.para.qbits_window, self.para.qbits_rmsd, 
                                   top=10, sec_struct=self.para.secstruct,
                                   antiparallel=first_recursion,
                                   min_nbrs=self.para.min_nbrs, contiguous=True)
            except:
                pass
        qreps = None
        if os.path.exists(outdir + '/qbit_reps/'):
            qreps = [outdir + '/qbit_reps/' + pdb_path for pdb_path in os.listdir(outdir + '/qbit_reps/')]
        return qreps

    def _prepare_seed_sses(self, qrep, full_sse_list):
        # if using an external query structure, include it in the Qbits 
        # search until the protein under construction has at least two SSEs
        if len(full_sse_list) < 2:
            all_reps = self.query_sse_list + full_sse_list + [qrep]
        else:
            all_reps = full_sse_list + [qrep]
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

    #Add new queries into the self.queues.
    def _add_seed_sse(self, i, qrep, j, seed_sse_list, full_sse_list, exclusion_pdb, recursion_order, outdir):
        #_queues = []
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
        _the_pdb = ''
        if len(_full_sse_list) > 2:
            _full_pdb = _workdir + '/full.pdb'
            pdbutils.merge_pdbs(_full_sse_list, _full_pdb, min_nbrs=self.para.min_nbrs)
            sat = pdbutils.satisfied_termini(_full_pdb, self.para.max_nc_dist)
            _the_pdb = _full_pdb
        else:
            sat = pdbutils.satisfied_termini(_seed_pdb, self.para.max_nc_dist)
            _the_pdb = _seed_pdb
        n_sat = np.sum(sat)
        if self.para.num_iter - recursion_order - 1 > n_sat:
            # if it is impossible to satisfy all N- or C- termini within 
            # the remaining number of iterations, exit the branch early
            #return _queues
            return
        # if recursion_order is not 1, continue adding qbit reps 
        if recursion_order > 1:
            _exclusion_pdb = _workdir + '/exclusion.pdb'
            pdbutils.merge_pdbs([exclusion_pdb, qrep], _exclusion_pdb, min_nbrs=self.para.min_nbrs)
            #_queues.append([_the_pdb, _exclusion_pdb, _full_sse_list, recursion_order - 1])
            self.queues.append([_the_pdb, _exclusion_pdb, _full_sse_list, recursion_order - 1])
        #return _queues
        return    
            
    ### NEW FUNCTIONS FOR GENERATING SSEs IN PARALLEL

    def _generate_proteins_parallel(self, recursion_order):
        manager = Manager()
        pool = Pool(6)
        queues = manager.Queue()

        counter = 1
        item = [self.pdbs[-1], self.exclusion_pdbs[-1], self.full_sse_list, recursion_order]
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

    #The current parallel is not working due to NULL 'qreps = self._generate_qreps(pdb, exclusion_pdb, recursion_order, outdir)'
    #query.py -> qbits_search() -> 'p.parse(outdir=outdir, show_parsing_progress=False)'
    #p.parse also apply parallel, which may cause a conflict.
    def _const_protein_parallel(self, item, queues):
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

                