'''
When the ligand is anchored in a specific position. 
We first generate all positions and conformations of the ligands on the anchored pocket.
We then sample vdMs on the backbone with the given ligand ensembles.
'''

import os
import sys
import prody as pr
import pandas as pd
import numpy as np

from rvrsprot.basic import constant
from rvrsprot.combs import gvdm_helper
from rvrsprot.combs import __search_cg_vdms as search_cg_vdms

pd.set_option("display.max_columns", None)


'''
python /mnt/e/GitHub_Design/Rvrsprot/scripts/sample_vdm/sample_vdm_with_anchor_ligand_yuda.py
'''

########################################################################################


vdm_cg_aa_atommap_dict = {
    (0, 0):{
        'cg' : 'coo',
        'lgd_sel' : ['O2', 'C26', 'O27'],
        'represent_name' : 'OD1',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['OD1', 'CG', 'OD2'],
    },    
        (0, 1):{
        'cg' : 'coo',
        'lgd_sel' : ['O2', 'C26', 'O27'],
        'represent_name' : 'OD1',
        'correspond_resname' : 'ASP',
        'correspond_names' : ['OD2', 'CG', 'OD1'],
    },    
        (0, 2):{
        'cg' : 'coo',
        'lgd_sel' : ['O2', 'C26', 'O27'],
        'represent_name' : 'OE1',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['OE1', 'CD', 'OE2'],
    },    
        (0, 3):{
        'cg' : 'coo',
        'lgd_sel' : ['O2', 'C26', 'O27'],
        'represent_name' : 'OE1',
        'correspond_resname' : 'GLU',
        'correspond_names' : ['OE2', 'CG', 'OE1'],
    },    
}


########################################################################################

workdir = '/mnt/e/DesignData/kaipeng_yuda/'
target_name = '5kau_carboxy.pdb'
target = pr.parsePDB(workdir + target_name)


########################################################################################

# def load_ligs(lig_path):
#     ligs = []
#     for lig_name in os.listdir(lig_path):
#         if not '.pdb' in lig_name:
#             continue
#         lig = pr.parsePDB(lig_path + lig_name)
#         ligs.append(lig)
#     return ligs

#ligands = load_ligs('/mnt/e/DesignData/Chemodrugs/HB_RUC__ligs/ligs_1st_tf2middle/')

ligands = [target.select('resname RHQ')]

path_to_vdm_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'

########################################################################################

def run_vdm_sample(target, resnum, ligands, path_to_vdm_database, vdm_cg_aa_atommap_dict, outdir):

    labels_cgs = {}
    df_cgs = {}
    dist_ind_cgs = {}

    cg_ids = vdm_cg_aa_atommap_dict.keys()
    for cg_id in cg_ids:
        #Load vdm database.
        #Here we specify to check ASP vdMs. Note this can be programed to be easier.
        df_vdm = gvdm_helper.load_old_vdm(path_to_vdm_database, vdm_cg_aa_atommap_dict[cg_id]['cg'], aa)
        
        #Transformation.
        pos = target.select('name N CA C and resnum ' + str(resnum))
        tf = pr.calcTransformation(constant.ideal_ala_coords, pos.getCoords())
        df_vdm.loc[:, ['c_x', 'c_y', 'c_z']] = tf.apply(df_vdm[['c_x', 'c_y', 'c_z']].to_numpy())

        #Test if the first vdm's coords are changed.
        #x = labels_cgs[cg_id].iloc[0]
        #v = df_vdm[(df_vdm['CG'] == x['CG']) & (df_vdm['rota'] == x['rota']) & (df_vdm['probe_name'] == x['probe_name'])]
        #ag = gvdm_helper.df2ag(v, 'test_position_of_1st_vdm')
        #pr.writePDB(outdir + ag.getTitle(), ag)

        #Search vdms.
        search_cg_vdms.search_vdm(df_vdm, ligands, cg_id, vdm_cg_aa_atommap_dict, labels_cgs, df_cgs, dist_ind_cgs, rmsd = 0.75)
        

    CgCombInfoDict = search_cg_vdms.construct_vdm_write(outdir, ligands, labels_cgs, vdm_cg_aa_atommap_dict, df_cgs, dist_ind_cgs, clash_radius = 2.3)

    search_cg_vdms.write_summary(outdir, CgCombInfoDict, name = '_summary.tsv')
    return

######################################################################################

'''
#Run one residue.
outdir = workdir + 'vdm_sampling_res90_ARG_coo/'
resnum = 90
aa = 'ARG'
os.makedirs(outdir, exist_ok= True)

run_vdm_sample(target, resnum, ligands, path_to_vdm_database, vdm_cg_aa_atommap_dict, outdir)
'''


resnums = [104, 104]
aas = ['ARG', 'LYS']

       
for i in range(len(resnums)):
    resnum = resnums[i]
    aa = aas[i]
    outdir = workdir + 'vdm_sampling_res' + str(resnum) + '_' + aa + '/'
    os.makedirs(outdir, exist_ok=True)
    run_vdm_sample(target, resnum, ligands, path_to_vdm_database, vdm_cg_aa_atommap_dict, outdir)