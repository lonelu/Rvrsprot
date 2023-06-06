'''
If you want to manually check all vdms.
'''

import os
import sys
import prody as pr
import pandas as pd
import numpy as np

from rvrsprot.basic import constant
from rvrsprot.combs import gvdm_helper


'''
python /mnt/e/GitHub_Design/Rvrsprot/scripts/sample_vdm/manual_all_vdm.py

'''


path_to_vdm_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'


workdir = '/mnt/e/DesignData/Chemodrugs/HB_RUC_1st_vdm/'
target_name = 'bb_prep_599.pdb'
target = pr.parsePDB(workdir + target_name)

outdir = workdir + 'vdm_all_bbcnh_ASP/'
os.makedirs(outdir, exist_ok= True)

gvdm_helper.write_vdm(path_to_vdm_database, 'bb_cnh', 'ASP', outdir, total = None)


tf = pr.calcTransformation(target.select('resnum 131 and name N CA C'), constant.ideal_ala_coords)
    
tf.apply(target)

pr.writePDB('_' + target_name, target_name)