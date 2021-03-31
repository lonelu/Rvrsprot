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

from smallprot import pdbutils, query, cluster_loops, smallprot_config, sse
from smallprot import logger, peputils, constant, plot, extract_master, struct_analysis


class ExtendSSE:
    '''
    To extend sses of exsiting structures.
    '''
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
        self.SSE = sse.SSE(seed_pdb, query_pdb, exclusion_pdb, _workdir, para)
        self.targetList = os.path.realpath(self.para.database) + '/pds_list_2p5.txt'
        self.chains_dict = os.path.realpath(self.para.database) + '/db_2p5A_0p3rfree_chains_dictionary.pkl'
        #--------------------------------
        self.workdir = _workdir

    def _generate_qreps(self, workdir,  pdb, exclusion_pdb):
        # search for a contiguous secondary structural element to add
        outfile = workdir + '/stdout'
        print('Querying MASTER')
        query.master_query(pdb, self.targetList, self.para.rmsdCut, 
            topN=self.para.master_query_top, outfile=outfile, clobber=False)

    def get_extention(self, extension_lengths = [(12, 0), (0, 12)]):
        self._generate_qreps(self.workdir, self.SSE.seed_pdb, self.SSE.exclusion_pdb)

        match = self.workdir + '/match.txt'
        seq = self.workdir + '/seq.txt'

        outdir = self.workdir + '/exsearch/'
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        qbits.exsearch.extend_search(outdir, None, self.SSE.seed_pdb, match, seq, extension_lengths)

        



