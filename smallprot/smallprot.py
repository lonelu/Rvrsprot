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
class SSE_info:
    ahull_in_ratio: float
    median_min_dists: np.ndarray

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
        #For job queue, to parallel jobs
        self.queues = []
        #--------------------------------
        self.pre_build_pdbs = []
        self.pre_full_sses = []
        self.pre_build_pdbs_summary = []


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

    ### FUNCTIONS

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

    ### FUNCTIONS FOR GENERATING SSEs

    def _generate_proteins(self):     
        self.queues.append([self.seed_pdb, self.exclusion_pdb, self.full_sse_list, self.para.num_iter])
        while len(self.queues) > 0:            
            queue = self.queues.pop(0)
            print('--Get queue--')
            self._const_protein(queue[0], queue[1], queue[2], queue[3]) 

        print('Finish construct sses!')
        if len(self.pre_build_pdbs) > 0:   
            filename = self.workdir + '/protein_build_summary.txt'                  
            #self._write_protein_summary(filename, self.pre_build_pdbs, self.pre_build_pdbs_summary)

            print('Construct sses found.')
        else:
            print('No construct sses found.')
            return

    def _write_protein_summary(self, filename, pdb_paths, summaries):
        with open(filename, 'w') as f:
            f.write('pdb_path\tca_in_ahull_ratio\tmedian_min_dist_diff\n')
            for i in range(len(pdb_paths)):
                median_min_dist = '\t'.join([ str(abs(x[0] -x[1])) for x in summaries[i][1]])
                f.write(pdb_paths[i] + '\t' + str(summaries[i][0]) +  '\t' + '\n')           

    def _const_protein(self, pdb, exclusion_pdb, full_sse_list, recursion_order):
        outdir = os.path.dirname(pdb)
        #Construct final protein.
        if recursion_order == 0:
            self._const_final_protein(outdir, full_sse_list)  
            return
        #Generate queries for next recursion.
        qreps = self._generate_qreps(pdb, exclusion_pdb, recursion_order, outdir)
        if qreps == None:
            return
        if True:
            qreps = self._resize_qreps(qreps, full_sse_list, self.para.loop_query_win)
        for i, qrep in enumerate(qreps):
            seed_sse_lists = self._prepare_seed_sses(qrep, full_sse_list)
            for j, seed_sse_list in enumerate(seed_sse_lists):
                self._add_seed_sse(i, qrep, j, seed_sse_list, full_sse_list, exclusion_pdb, recursion_order, outdir)
        return
   
    def _const_final_protein(self, outdir, full_sse_list):
        '''
        The smallprot building process is at the last step as the recursion_order == 0. 
        Try to construct the protein and check if it satisfy ahull and distance limitations.
        outdir : output path of the final_protein
        full_sse_list : all sses for the final protein
        '''
        #Merge the sses to generate the pdb.
        folders = outdir.split('/')
        filename = 'full'
        for i in range(self.para.num_iter, 0, -1):
            filename += '_' + folders[-i]
        filename += '.pdb'
        _full_pdb = outdir + '/' + filename
        pdbutils.merge_pdbs(full_sse_list, _full_pdb, min_nbrs=self.para.min_nbrs)

        #Filter final pdbs with limitations.
        #COMMENT: Note that to use alpha_hull for the helix bundles, the helix must be cutted properly. 
        #         Otherwise, the alpha_hull will be out of control.
        ahull_in_ratio = struct_analysis.get_in_ahull_ratio(_full_pdb)
        # if ahull_in_ratio < 0.5:
        #     return
        distances = peputils.get_neighbor_dists(full_sse_list, list(range(len(full_sse_list))))
        
        #Copy out pre_build_pdbs into one folder, so it is easy to analyze them together.
        pre_build_workdir = self.workdir + '/pre_build_pdbs'
        if not os.path.exists(pre_build_workdir):
            os.mkdir(pre_build_workdir)
        _full_pdb_new = pre_build_workdir + '/' + filename
        shutil.copyfile(_full_pdb, _full_pdb_new)
        
        self.pre_build_pdbs.append(_full_pdb_new)
        self.pre_full_sses.append(full_sse_list)  
        self.pre_build_pdbs_summary.append([ahull_in_ratio, distances])

    def _const_prot_loop(self, pdb, full_sse_list, outdir):
        '''
        #TO DO: This function requires changes.
        '''
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
            topN=self.para.master_query_top, outfile=outfile, clobber=False)
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
                                   self.para.top, sec_struct=self.para.secstruct,
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
            try:
                struct = pdbutils._cut_pdbs(two_sses, loop_query_win)
                chain = list(struct.get_chains())[-1]     
                re_qrep_path = os.path.dirname(qrep) + '_resize'
                if not os.path.exists(re_qrep_path):    
                    os.mkdir(re_qrep_path)
                re_qrep = re_qrep_path + '/resize_' + os.path.basename(qrep)
                pdbutils.save_struct(chain, re_qrep)
                resize_qreps.append(re_qrep)
            except:
                #There is a bug here in the self.cut_pdbs. Try prody???
                print('resize crash!')
                resize_qreps.append(qrep)
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
