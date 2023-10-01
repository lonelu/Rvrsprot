'''
When the ligand is positioned.
We sample the backbone with rvdms.
The backbone could then be used for RFDiffuse.
'''

import os
import sys
import prody as pr
import pandas as pd
import numpy as np

from rvrsprot.combs import gvdm_helper
from rvrsprot.combs import __supper_cg_rvdms

pd.set_option("display.max_columns", None)


########################################################################################

vdm_cg_aa_atommap_dict_bbcco= {
    (0, 0):{
        'cg' : 'bb_cco',
        'lgd_sel' : ['C13', 'C8', 'O3'],
        'represent_name' : 'O',
        'correspond_resname' : 'GLY',
        'correspond_names' : ['CA', 'C', 'O'],
    },  
}

########################################################################################

workdir = '/Users/lonelu/DesignData/Combs/vdM_paper/'
outdir = workdir + 'out_rvdm/'
os.makedirs(outdir, exist_ok= True)

target_name = '6w70_mono.pdb'
target = pr.parsePDB(workdir + target_name)

#lig_path = workdir + 'ligand.pdb' 
ligand = pr.parsePDB(workdir + target_name).select('resname GG2')

aa = 'HIS'

########################################################################################

vdm_cg_aa_atommap_dict = vdm_cg_aa_atommap_dict_bbcco
path_to_rvdm_database='/Users/lonelu/DesignData/Database/rvdMs/'

########################################################################################

def run_rvdm_sample(outdir, ligand, vdm_cg_aa_atommap_dict, path_to_rvdm_database, aa, target = None):

    cg_ids = vdm_cg_aa_atommap_dict.keys()

    for cg_id in cg_ids:
        #Load vdm database.
        df_vdm = gvdm_helper.load_old_vdm(path_to_rvdm_database, vdm_cg_aa_atommap_dict[cg_id]['cg'], aa)
        
        #Transformation.
        lgd_sel_atomnames = vdm_cg_aa_atommap_dict[cg_id]['lgd_sel']
        cg_sel_atomnames = vdm_cg_aa_atommap_dict[cg_id]['correspond_names']
        __supper_cg_rvdms.rvdm_transform(ligand, lgd_sel_atomnames, cg_sel_atomnames, df_vdm)

        #Write top 10 rvdms.
        labels_cgs = df_vdm.sort_values(by='cluster_size', ascending=False)[['CG', 'rota', 'probe_name']].drop_duplicates()
        #print(df_vdm['cluster_size'].value_counts())

        for i in range(10):
            x = labels_cgs.iloc[i]
            v = df_vdm[(df_vdm['CG'] == x['CG']) & (df_vdm['rota'] == x['rota']) & (df_vdm['probe_name'] == x['probe_name'])]
            print(v['cluster_number'].iloc[0])
            ag = gvdm_helper.df2ag(v, 'test_rvdm_' + str(i))
            pr.writePDB(outdir + ag.getTitle(), ag)
    return

######################################################################################
if __name__ == '__main__':
    
    run_rvdm_sample(outdir, ligand, vdm_cg_aa_atommap_dict, path_to_rvdm_database, aa, target = None)