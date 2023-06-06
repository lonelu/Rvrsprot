import os
import sys
import prody as pr
import pandas as pd
import numpy as np

from rvrsprot.basic import constant
from rvrsprot.combs import gvdm_helper
from rvrsprot.combs import __search_cg_vdms

workdir = '/mnt/e/DesignData/lab_designs/vdM_paper/'
target_name = '6w70_h2o.pdb'
target = pr.parsePDB(workdir + target_name)

workdir = '/mnt/e/DesignData/Database/vdMh2o/'
outdir = workdir + 'vdm_map_6W70_Y46/'
os.makedirs(outdir, exist_ok= True)

tf = pr.calcTransformation(target.select('resnum 46 and name N CA C'), constant.ideal_ala_coords)
tf.apply(target)
pr.writePDB(outdir + 'tf_'+ target_name , target)


h2ovdm_db = workdir + 'vdm_pdb_TYR-pro/'
vdms = []
no_gly_vdms = []
for v in os.listdir(h2ovdm_db):
    if not '.pdb' in v:
        continue
    vdm = pr.parsePDB(h2ovdm_db + v)
    if not vdm.select('resname PRO and name CA C O'):
        no_gly_vdms.append(v)
        continue
    vdms.append(vdm)



# Get the ligand from the protein
ligand_coords = []
#x = np.array([target.select('resname GG2 and name ' + ss).getCoords()[0] for ss in ['C23', 'C19', 'O2']])
x = np.array([target.select('resname GG2 and name ' + ss).getCoords()[0] for ss in ['C19', 'O2']])

ligand_coords.append(x)
ligand_coords = np.array(ligand_coords).reshape((-1, 6))



vdm_cg_aa_atommap_dict_gly = {
    (0, 0):{
        'cg' : 'hid',
        'lgd_sel' : ['C23', 'C19', 'O2'],
        'represent_name' : 'O',
        'correspond_resname' : 'ALA',
        'correspond_names' : ['CA', 'C', 'O'],
    },
}

vdm_cg_aa_atommap_dict_gly = {
    (0, 0):{
        'cg' : 'hid',
        'lgd_sel' : ['C19', 'O2'],
        'represent_name' : 'O',
        'correspond_resname' : 'PRO',
        'correspond_names' : ['C', 'O'],
    },
}

rmsd = 0.75
num_cg_atoms = 2
radius = np.sqrt(num_cg_atoms) * rmsd
matched = []
for cg_id in vdm_cg_aa_atommap_dict_gly.keys():
    coords_a = []
    for vdm in vdms:
        coords_a.append(np.array([vdm.select('resname PRO and name ' + ss).getCoords()[0] for ss in vdm_cg_aa_atommap_dict_gly[cg_id]['correspond_names']]))
    coords_a = np.array(coords_a).reshape((-1, 6))

    dists, inds = __search_cg_vdms.get_nearest_vdms_rmsd(ligand_coords, coords_a, radius = radius)
    for i in range(len(inds)):
        if len(inds[i]) <= 0:
            continue
        matched.append((vdms[i], dists[i]))


for m in matched:
    pr.writePDB(outdir + m[0].getTitle() + '_' + str(round(min(m[1])/np.sqrt(num_cg_atoms), 2)) + '.pdb', m[0])



