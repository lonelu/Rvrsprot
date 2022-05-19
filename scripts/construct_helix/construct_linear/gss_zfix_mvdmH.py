'''
Add helix to metal-vdms or mvdms, and superimpose on the std motif rotations.
Keep the mvdmH is z-fixed. 
z-fix is the helix is in the same direction with z-axis.
'''

import prody as pr
import os
from scipy.spatial.transform import Rotation
from numpy.linalg import norm
import math
from itertools import permutations
import numpy as np
from sklearn.neighbors import NearestNeighbors
from metalprot.database import database_cluster
from metalprot.basic import prody_ext
from metalprot.basic import vdmer
import shutil
import sys
sys.path.append(r'/mnt/e/GitHub_Design/Rvrsprot/scripts/construct_helix/construct_linear')
from gss_rvdmH_basics import clash, clash_rvs, cal_angle, calc_z_direction, rot_rvdmH_C2


### For C2, we only want to rotate around z-axis when the helix direction is also in agree with z-axis.

def generate_z_constrained_helix(helix_std, vdm, imido_rot_dir, helix_std_outdir):
    zn = pr.AtomGroup()
    zn.setCoords(np.zeros((1, 3), dtype=float))
    zn.setNames(['ZN'])
    zn.setResnums([1000])
    zn.setResnames(['ZN'])
    for x in os.listdir(imido_rot_dir):
        if not '.pdb' in x or not 'A' in x:
            continue
        std_cp = pr.parsePDB(imido_rot_dir + x)
        std_coords = [std_cp.select('serial ' + str(i)).getCoords()[0] for i in [4, 5, 2, 3, 1]]
        vdm_cp = vdm.copy()
        helix_std_cp = helix_std.copy()
        pr.calcTransformation(vdm_cp.select('name CG ND1 CD2 CE1 NE2').getCoords(), np.array(std_coords)).apply(vdm_cp)

        pr.calcTransformation(helix_std_cp.select('resindex 7 and name N C CA').getCoords(), vdm_cp.select('resindex 1 and name N C CA')).apply(helix_std_cp)
        
        if clash(vdm_cp.select('name CG ND1 CD2 CE1 NE2').getCoords(), helix_std_cp.select('name N C CA O').getCoords(), 2.5):
            continue
        angle, angle_zy, angle_zx = calc_z_direction(helix_std_cp)
        print('angle: {}, angle_zy: {}'.format(angle, angle_zy))
        if ( angle < 30 and (angle_zy <20 or angle_zy > 160)) or (angle > 150 and (angle_zy <20 or angle_zy > 160)):
            ind_vdm_dict = {7:vdm_cp}
            title = vdm.getTitle() + '_' + x.split('.')[0] + '_a_' + str(round(angle, 0)) + '_azy_' + str(round(angle_zy, 0))
            # There is two ways to add Metal, we can add the original ZN or add ZN at [0, 0, 0]
            #new_pdb = prody_ext.mutate_vdm_target_into_ag2(helix_std_cp, ind_vdm_dict, title, vdm_sel='resindex 1 or name ZN')
            new_pdb = prody_ext.mutate_vdm_target_into_ag2(helix_std_cp, ind_vdm_dict, title, vdm_sel='resindex 1')
            new_pdb = prody_ext.combine_ags([new_pdb, zn], new_pdb.getTitle(), ['A', 'M'])
            pr.writePDB(helix_std_outdir + title, new_pdb)
    return


def generate_z_constrained_all_helixs(helix_std, vdm_dir, imido_rot_dir, helix_std_outdir):
    '''
    Note that the current method is specified for NE2 contacting Metal. 
    There are 133 ND1 contacting vdMs. Such as AAA_cluster_518_sc_-1.0, AAA_cluster_520_sc_-1.0, AAA_cluster_821_sc_-4.83.
    '''
    vdm_ND1s = []
    os.makedirs(helix_std_outdir, exist_ok=True)
    for x in os.listdir(vdm_dir):
        if not '.pdb' in x:
            continue
        vdm = pr.parsePDB(vdm_dir + x)
        if vdmer.get_contact_atom(vdm).getName() == 'ND1':
            vdm_ND1s.append(vdm)
            continue
        generate_z_constrained_helix(helix_std, vdm, imido_rot_dir, helix_std_outdir)
    return vdm_ND1s


workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2_coil/'
helix_std = pr.parsePDB(workdir + 'CC_15-Mer.pdb')

vdm_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211015_AAext3/AAA_H_vdms/'
imido_rot_dir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2_coil/std_rots/'
helix_std_outdir = workdir + '_c2_z_fix/'

vdm_ND1s = generate_z_constrained_all_helixs(helix_std, vdm_dir, imido_rot_dir, helix_std_outdir)

'''
vdm = pr.parsePDB(workdir + 'AAA_cluster_0_mem_71_centroid_2007_2j7j_ZN_3_mem0_HIS_A_80.pdb')

generate_z_constrained_helix(helix_std, vdm, imido_rot_dir, helix_std_outdir)

points = np.array([[0, 0, 0], [0, 0, 2], [0, 2, 0], [2, 0, 0]], dtype=float)
prody_ext.write2pymol(points, workdir, filename='axis',  names=['O', 'H', 'H', 'H'])
'''