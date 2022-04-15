


import os
import numpy as np
import prody as pr


workdir = '/mnt/e/DesignData/tm/'

pdbpath = 'hp1_75_start_model_renumbered_2.pdb'

pdb = pr.parsePDB(workdir + pdbpath)

tm_sel = ''

center = pr.calcCenter(pdb.select('name CA'))

pdb_cr = pdb.copy()

pr.calcTransformation(center.reshape(1, 3), np.zeros((1,3))).apply(pdb_cr)

pr.calcCenter(pdb_cr.select('name CA'))

pr.writePDB(workdir + pdb.getTitle() + '_cr.pdb', pdb_cr)