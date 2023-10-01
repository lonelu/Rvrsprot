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
python /Users/lonelu/GitHub_Design/Rvrsprot/scripts/sample_vdm/sample_vdm_with_anchor_ligand.py
'''

########################################################################################

vdm_cg_aa_atommap_dict_bbcnh = {
    (0, 0):{
        'cg' : 'bb_cnh',
        'lgd_sel' : ['C18', 'N3', 'H15'],
        'represent_name' : 'N',
        'correspond_resname' : 'GLY',
        'correspond_names' : ['CA', 'N', 'H'],
    },
    (0, 1):{
        'cg' : 'bb_cnh',
        'lgd_sel' : ['C19', 'N3', 'H15'],
        'represent_name' : 'N',
        'correspond_resname' : 'GLY',
        'correspond_names' : ['CA', 'N', 'H']
    },    
    (0, 2):{
        'cg' : 'bb_cnh',
        'lgd_sel' : ['C18', 'N3', 'H15'],
        'represent_name' : 'N',
        'correspond_resname' : 'ALA',
        'correspond_names' : ['CA', 'N', 'H']
    },
    (0, 3):{
        'cg' : 'bb_cnh',
        'lgd_sel' : ['C19', 'N3', 'H15'],
        'represent_name' : 'N',
        'correspond_resname' : 'ALA',
        'correspond_names' : ['CA', 'N', 'H']
    },    
}


vdm_cg_aa_atommap_dict_his = {
    (0, 0):{
        'cg' : 'hid',
        'lgd_sel' : ['C4', 'N1', 'C5'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CG', 'ND1', 'CE1'],
    },
    (0, 1):{
        'cg' : 'hid',
        'lgd_sel' : ['C5', 'N1', 'C4'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CG', 'ND1', 'CE1'],
    },    
    (0, 2):{
        'cg' : 'hie',
        'lgd_sel' : ['C4', 'N1', 'C5'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CE1', 'NE2', 'CD2'],
    },
    (0, 3):{
        'cg' : 'hie',
        'lgd_sel' : ['C5', 'N1', 'C4'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CE1', 'NE2', 'CD2'],
    },    
}


vdm_cg_aa_atommap_dict_indole = {
    (0, 0):{
        'cg' : 'indole',
        'lgd_sel' : ['N1', 'C2', 'C12'],
        'represent_name' : 'NE1',
        'correspond_resname' : 'TRP',
        'correspond_names' : ['NE1', 'CD2', 'CH2'],
    },    
}


vdm_cg_aa_atommap_dict_ph ={
    (0, 0):{
        'cg' : 'ph',
        'lgd_sel' : ['C5', 'C6', 'C12'],
        'represent_name' : 'CZ',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CD1', 'CD2', 'CZ'],
    },    
    (0, 1):{
        'cg' : 'ph',
        'lgd_sel' : ['C5', 'C6', 'C12'],
        'represent_name' : 'CZ',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CD2', 'CD1', 'CZ'],
    },    
    (0, 2):{
        'cg' : 'ph',
        'lgd_sel' : ['C5', 'C6', 'C12'],
        'represent_name' : 'CZ',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CZ', 'CD1', 'CD2'],
    },  
    (0, 3):{
        'cg' : 'ph',
        'lgd_sel' : ['C5', 'C6', 'C12'],
        'represent_name' : 'CZ',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CZ', 'CD2', 'CD1'],
    },  
    (0, 4):{
        'cg' : 'ph',
        'lgd_sel' : ['C5', 'C6', 'C12'],
        'represent_name' : 'CZ',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CD1', 'CZ', 'CD2'],
    },  
    (0, 5):{
        'cg' : 'ph',
        'lgd_sel' : ['C5', 'C6', 'C12'],
        'represent_name' : 'CZ',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CD2', 'CZ', 'CD1'],
    },  

    # (1, 0):{
    #     'cg' : 'ph',
    #     'lgd_sel' : ['C5', 'C6', 'C12'],
    #     'represent_name' : 'CG',
    #     'correspond_resname' : 'PHE',
    #     'correspond_names' : ['CE1', 'CE2', 'CG'],
    # },    
    # (1, 1):{
    #     'cg' : 'ph',
    #     'lgd_sel' : ['C5', 'C6', 'C12'],
    #     'represent_name' : 'CG',
    #     'correspond_resname' : 'PHE',
    #     'correspond_names' : ['CE2', 'CE1', 'CG'],
    # },    
    # (1, 2):{
    #     'cg' : 'ph',
    #     'lgd_sel' : ['C5', 'C6', 'C12'],
    #     'represent_name' : 'CG',
    #     'correspond_resname' : 'PHE',
    #     'correspond_names' : ['CE1', 'CG', 'CE2'],
    # },    
    # (1, 3):{
    #     'cg' : 'ph',
    #     'lgd_sel' : ['C5', 'C6', 'C12'],
    #     'represent_name' : 'CG',
    #     'correspond_resname' : 'PHE',
    #     'correspond_names' : ['CE2', 'CG', 'CE1'],
    # },    

    (6, 0):{
        'cg' : 'phenol',
        'lgd_sel' : ['C5', 'C6', 'C12'],
        'represent_name' : 'CZ',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CD1', 'CD2', 'CZ'],
    },    
    (6, 1):{
        'cg' : 'phenol',
        'lgd_sel' : ['C5', 'C6', 'C12'],
        'represent_name' : 'CZ',
        'correspond_resname' : 'PHE',
        'correspond_names' : ['CD2', 'CD1', 'CZ'],
    },    

}


########################################################################################

#workdir = '/Users/lonelu/DesignData/Chemodrugs/HB_RUC_1st_vdm/'
#target_name = 'bb_prep_599.pdb'
workdir = '/Users/lonelu/DesignData/Chemodrugs/RUC_vdm/HB_RUC_rucA_vdm/'
target_name = 'ruc_A.pdb'
target = pr.parsePDB(workdir + target_name)

outdir = workdir + 'vdm_sampling_res29_ASP_hidhie2/'
resnum = 29
aa = 'ASP'
vdm_cg_aa_atommap_dict = vdm_cg_aa_atommap_dict_bbcnh

# outdir = workdir + 'vdm_sampling_res90_LEU_ph/'
# resnum = 90
# aa = 'LEU'
# vdm_cg_aa_atommap_dict = vdm_cg_aa_atommap_dict_ph

# outdir = workdir + 'vdm_sampling_res131_ASP_hidhie2/'
# resnum = 131
# aa = 'ASP'
# vdm_cg_aa_atommap_dict = vdm_cg_aa_atommap_dict_his

os.makedirs(outdir, exist_ok= True)
########################################################################################

def load_ligs(lig_path):
    ligs = []
    for lig_name in os.listdir(lig_path):
        if not '.pdb' in lig_name:
            continue
        lig = pr.parsePDB(lig_path + lig_name)
        ligs.append(lig)
    return ligs

#ligands = load_ligs('/Users/lonelu/DesignData/Chemodrugs/HB_RUC__ligs/ligs_1st_tf2middle/')
ligands = load_ligs('/Users/lonelu/DesignData/Chemodrugs/HB_RUC__ligs/ligs_rucA_tf2middle/')
path_to_vdm_database='/Users/lonelu/DesignData/Combs/vdMs/'

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

run_vdm_sample(target, resnum, ligands, path_to_vdm_database, vdm_cg_aa_atommap_dict, outdir)