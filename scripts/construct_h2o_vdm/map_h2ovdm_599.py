import os
import sys
import prody as pr
import pandas as pd
import numpy as np

from rvrsprot.basic import constant
from rvrsprot.combs import gvdm_helper
from rvrsprot.combs import __search_cg_vdms

workdir = '/mnt/e/DesignData/Chemodrugs/HB_RUC_1st_vdm/'
target_name = 'bb_prep_599.pdb'
target = pr.parsePDB(workdir + target_name)

workdir = '/mnt/e/DesignData/Database/vdMh2o/'
outdir = workdir + 'vdm_map_599_D131/'
os.makedirs(outdir, exist_ok= True)

tf = pr.calcTransformation(target.select('resnum 131 and name N CA C'), constant.ideal_ala_coords)
tf.apply(target)
pr.writePDB(outdir + 'tf_' + target_name, target)


h2ovdm_db = workdir + 'vdm_pdb_ASP-his/'
vdms = []
no_his_vdms = []
for v in os.listdir(h2ovdm_db):
    if not '.pdb' in v:
        continue
    vdm = pr.parsePDB(h2ovdm_db + v)
    if not vdm.select('resname HIS and name CG ND1 CE1'):
        no_his_vdms.append(v)
        continue
    vdms.append(vdm)



'''
# Get the ligand from the protein
coords_lig = []
x = np.array([target.select('resname RUC and name ' + ss).getCoords()[0] for ss in ['C4', 'N1', 'C5']])
coords_lig.append(x)
x = np.array([target.select('resname RUC and name ' + ss).getCoords()[0] for ss in ['C5', 'N1', 'C4']])
coords_lig.append(x)
coords_lig = np.array(coords_lig).reshape((-1, 9))
'''

def load_ligs(lig_path):
    ligs = []
    for lig_name in os.listdir(lig_path):
        if not '.pdb' in lig_name:
            continue
        lig = pr.parsePDB(lig_path + lig_name)
        tf.apply(lig)
        ligs.append(lig)
    return ligs

ligands = load_ligs('/mnt/e/DesignData/Chemodrugs/HB_RUC__ligs/ligs_1st_tf2middle/')


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



vdm_cg_aa_atommap_dict_his = {
    (0, 0):{
        'cg' : 'hid',
        'lgd_sel' : ['C4', 'N1'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CG', 'ND1'],
    },
    (0, 1):{
        'cg' : 'hid',
        'lgd_sel' : ['C5', 'N1'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CG', 'ND1'],
    },    
    (0, 2):{
        'cg' : 'hie',
        'lgd_sel' : ['C4', 'N1'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CE1', 'NE2'],
    },
    (0, 3):{
        'cg' : 'hie',
        'lgd_sel' : ['C5', 'N1'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CE1', 'NE2'],
    },    
    (1, 0):{
        'cg' : 'hid',
        'lgd_sel' : ['C4', 'N1'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CE1', 'ND1'],
    },
    (1, 1):{
        'cg' : 'hid',
        'lgd_sel' : ['C5', 'N1'],
        'represent_name' : 'ND1',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CE1', 'ND1'],
    },    
    (1, 2):{
        'cg' : 'hie',
        'lgd_sel' : ['C4', 'N1'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CD2', 'NE2'],
    },
    (1, 3):{
        'cg' : 'hie',
        'lgd_sel' : ['C5', 'N1'],
        'represent_name' : 'NE2',
        'correspond_resname' : 'HIS',
        'correspond_names' : ['CD2', 'NE2'],
    },    
}


rmsd = 0.75
num_cg_atoms = 2
radius = np.sqrt(num_cg_atoms) * rmsd
matched = []
for cg_id in vdm_cg_aa_atommap_dict_his.keys():
    ligand_coords = __search_cg_vdms.get_ligand_coords(ligands, vdm_cg_aa_atommap_dict_his[cg_id]['lgd_sel'])
    coords_a = []
    for vdm in vdms:
        coords_a.append(np.array([vdm.select('resname HIS and name ' + ss).getCoords()[0] for ss in vdm_cg_aa_atommap_dict_his[cg_id]['correspond_names']]))
    coords_a = np.array(coords_a).reshape((-1, 6))

    dists, inds = __search_cg_vdms.get_nearest_vdms_rmsd(ligand_coords, coords_a, radius = radius)
    for i in range(len(inds)):
        if len(inds[i]) <= 0:
            continue
        matched.append((vdms[i], dists[i]))


for m in matched:
    pr.writePDB(outdir + m[0].getTitle() + '_' + str(round(min(m[1])/np.sqrt(num_cg_atoms), 2)) + '.pdb', m[0])



