import os
import sys
import gzip
import string
import shutil
import numpy as np
import prody as pr
from smallprot import pdbutils, query, cluster_loops, smallprot_config


class SSE:
    """

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
    """

    def __init__(self, seed_pdb, query_pdb, exclusion_pdb,  workdir, para):
        self._prepare_pdbs(seed_pdb, query_pdb, exclusion_pdb,  workdir, para.min_nbrs)
        #--------------------------------   


    def _prepare_pdbs(self, seed_pdb, query_pdb, exclusion_pdb, _workdir, min_nbrs):
        # if necessary, split query pdb file into chains
        if query_pdb:
            query_chain_dir = _workdir + '/query_chains'
            if not os.path.exists(query_chain_dir):
                os.mkdir(query_chain_dir)
            pdbutils.split_pdb(query_pdb, query_chain_dir, set_bfac=np.log(min_nbrs))
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
                pdbutils.merge_pdbs([seed_pdb], _seed_pdb, set_bfac=np.log(min_nbrs))      
        elif query_pdb:
            pdbutils.merge_pdbs([query_pdb], _seed_pdb, set_bfac=np.log(min_nbrs))
            self.full_sse_list = []
        self.seed_pdb = _seed_pdb        
        if not query_pdb and not seed_pdb:
            raise AssertionError('Must provide either query_pdb or seed_pdb.') 
        # if necessary, determine path to exclusion PDB file
        if exclusion_pdb:
            _exclusion_pdb = _workdir + '/exclusion.pdb'
            pdbutils.merge_pdbs([_seed_pdb, exclusion_pdb], _exclusion_pdb, 
                                set_bfac=np.log(min_nbrs))
            self.orig_exclusion = exclusion_pdb
        else:
            _exclusion_pdb = _seed_pdb
            self.orig_exclusion = None
        self.exclusion_pdb = _exclusion_pdb