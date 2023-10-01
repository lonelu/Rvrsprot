import os
import sys
import prody as pr
import pandas as pd
import numpy as np

from rvrsprot.basic import constant
from rvrsprot.combs import gvdm_helper
from rvrsprot.combs import __search_cg_vdms

workdir = '/Users/lonelu/DesignData/Chemodrugs/_structure_split/'
target_name = 'ruc_A.pdb'
target = pr.parsePDB(workdir + target_name)

workdir = '/Users/lonelu/DesignData/Database/vdMh2o/'
outdir = workdir + 'vdm_map_599_D131_indole/'
os.makedirs(outdir, exist_ok= True)

tf = pr.calcTransformation(target.select('resnum 131 and name N CA C'), constant.ideal_ala_coords)
tf.apply(target)
pr.writePDB(outdir + 'tf_' + target_name, target)


h2ovdm_db = workdir + 'vdm_pdb_ASP_indole/'
vdms = []
no_his_vdms = []
for v in os.listdir(h2ovdm_db):
    if not '.pdb' in v:
        continue
    vdm = pr.parsePDB(h2ovdm_db + v)
    if not vdm.select('resname TRP and name CD1 NE1 CE2'):
        no_his_vdms.append(v)
        continue
    vdms.append(vdm)


# def load_ligs(lig_path):
#     ligs = []
#     for lig_name in os.listdir(lig_path):
#         if not '.pdb' in lig_name:
#             continue
#         lig = pr.parsePDB(lig_path + lig_name)
#         tf.apply(lig)
#         ligs.append(lig)
#     return ligs

# ligands = load_ligs('/Users/lonelu/DesignData/Chemodrugs/HB_RUC__ligs/ligs_1st_tf2middle/')

ligands = [target.select('resname RPW')]

vdm_cg_aa_atommap_dict_his = {
    (0, 0):{
        'cg' : 'hid',
        'lgd_sel' : ['C8', 'N23', 'C25'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CG', 'ND1', 'CE1'],
    },
    (0, 1):{
        'cg' : 'hid',
        'lgd_sel' : ['C25', 'N23', 'C8'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CG', 'ND1', 'CE1'],
    },    
    (0, 2):{
        'cg' : 'hie',
        'lgd_sel' : ['C8', 'N23', 'C25'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CE1', 'NE2', 'CD2'],
    },
    (0, 3):{
        'cg' : 'hie',
        'lgd_sel' : ['C25', 'N23', 'C8'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CE1', 'NE2', 'CD2'],
    },    
}



vdm_cg_aa_atommap_dict_indole = {
    (1, 0):{
        'cg' : 'indole',
        'lgd_sel' : ['C25', 'N23', 'C8'],
        'represent_name' : 'NE1',
        'correspond_resname' : 'TRP',
        'correspond_names' : ['CD1', 'NE1', 'CE2'],
    },
    (1, 1):{
        'cg' : 'indole',
        'lgd_sel' : ['C25', 'N23', 'C8'],
        'represent_name' : 'NE1',
        'correspond_resname' : 'TRP',
        'correspond_names' : ['CE2', 'NE1', 'CD1'],
    },
}


vdm_cg_aa_atommap_dict_his_2atom = {
    (0, 0):{
        'cg' : 'hid',
        'lgd_sel' : ['C8', 'N23'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CG', 'ND1'],
    },
    (0, 1):{
        'cg' : 'hid',
        'lgd_sel' : ['C25', 'N23'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CG', 'ND1'],
    },    
    (0, 2):{
        'cg' : 'hie',
        'lgd_sel' : ['C8', 'N23'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CE1', 'NE2'],
    },
    (0, 3):{
        'cg' : 'hie',
        'lgd_sel' : ['C25', 'N23'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CE1', 'NE2'],
    },    
    (1, 0):{
        'cg' : 'hid',
        'lgd_sel' : ['C8', 'N23'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CE1', 'ND1'],
    },
    (1, 1):{
        'cg' : 'hid',
        'lgd_sel' : ['C25', 'N23'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CE1', 'ND1'],
    },    
    (1, 2):{
        'cg' : 'hie',
        'lgd_sel' : ['C8', 'N23'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CD2', 'NE2'],
    },
    (1, 3):{
        'cg' : 'hie',
        'lgd_sel' : ['C25', 'N23'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CD2', 'NE2'],
    },    
}

vdm_cg_aa_atommap_dict = vdm_cg_aa_atommap_dict_his
vdm_cg_aa_atommap_dict = vdm_cg_aa_atommap_dict_indole

rmsd = 0.75
num_cg_atoms = 2
radius = np.sqrt(num_cg_atoms) * rmsd
matched = []
for cg_id in vdm_cg_aa_atommap_dict.keys():
    ligand_coords = __search_cg_vdms.get_ligand_coords(ligands, vdm_cg_aa_atommap_dict[cg_id]['lgd_sel'])
    coords_a = []
    for vdm in vdms:
        coords_a.append(np.array([vdm.select('resname TRP and name ' + ss).getCoords()[0] for ss in vdm_cg_aa_atommap_dict[cg_id]['correspond_names']]))
    coords_a = np.array(coords_a).reshape((-1, 9))

    dists, inds = __search_cg_vdms.get_nearest_vdms_rmsd(ligand_coords, coords_a, radius = radius)
    for i in range(len(inds)):
        if len(inds[i]) <= 0:
            continue
        matched.append((vdms[i], dists[i]))


for m in matched:
    pr.writePDB(outdir + m[0].getTitle() + '_' + str(round(min(m[1])/np.sqrt(num_cg_atoms), 2)) + '.pdb', m[0])



