'''
If you want to manually check all/selected vdms.
'''

import os
import sys
import prody as pr
import pandas as pd
import numpy as np

from rvrsprot.basic import constant
from rvrsprot.combs import gvdm_helper


'''
python /Users/lonelu/GitHub_Design/Rvrsprot/scripts/sample_vdm/manual_vdm_cluster.py

'''


path_to_vdm_database='/Users/lonelu/DesignData/Combs/vdMs/'


workdir = '/Users/lonelu/DesignData/Chemodrugs/HB_RUC_1st_vdm/'
target_name = 'bb_prep_599.pdb'
target = pr.parsePDB(workdir + target_name)

outdir = workdir + 'vdm_ASP131_clux_manual/'
os.makedirs(outdir, exist_ok= True)
cg = 'hid'
aa = 'ASP'
df_vdm = gvdm_helper.load_old_vdm(path_to_vdm_database, cg, aa)

#df_vdm[['CG', 'rota', 'probe_name']].iloc[0]
#hie
#2_2_3VGH_biomol_1_A_A
#3_2_3WQ0_biomol_1_A_A
#hid
#2_2_4GB5_biomol_1_A_A
#2_1_4NVR_biomol_1_A_A
i_cg =  df_vdm[(df_vdm['CG'] == 2) 
                & (df_vdm['rota'] == 1) 
                & (df_vdm['probe_name'] == '4NVR_biomol_1_A_A')]

cluster_num = i_cg['cluster_number'].iloc[0]
print(cluster_num)

gvdm_helper.write_vdm_cluster(df_vdm, cg, aa, outdir, cluster_num)


tf = pr.calcTransformation(target.select('resnum 131 and name N CA C'), constant.ideal_ala_coords)
    
tf.apply(target)

pr.writePDB(outdir + target_name, target)


