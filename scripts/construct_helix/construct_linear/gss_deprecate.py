'''
Some deprecated functions used similar as construct tetrahedral geo metal binding proteins.
The idea is first rotate the two sc of std his binding motif. 
Then superimpose the best vdm on the sc. 
The method is limited with vdms used for construction and is not well suited for the problem.
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
from .gss_rvdmH_basics import clash, clash_rvs, cal_angle, calc_z_direction, rot_rvdmH_C2


### Generate vdm-helix rotate with different angles.

def addhelix_tetrahedral_sc_rots(workdir, outdir, outdir_name, helix_std_name):
    '''
    Add vdm-helix into each component sc of std tetrahydral geo. 
    '''
    helix_std_outdir = workdir + outdir_name
    os.makedirs(helix_std_outdir, exist_ok=True)

    helix_std = pr.parsePDB(workdir + helix_std_name)

    for x in os.listdir(outdir):
        if not '.pdb' in x:
            continue
        std_cp = pr.parsePDB(outdir + x)
        std_coords = [std_cp.select('serial ' + str(i)).getCoords()[0] for i in [4, 5, 2, 3, 1]]

        helix_std_cp = helix_std.copy()
        pr.calcTransformation(helix_std_cp.select('name CG ND1 CD2 CE1 NE2').getCoords(), np.array(std_coords)).apply(helix_std_cp)
        pr.writePDB(helix_std_outdir + x, helix_std_cp)
    return

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2/'
outdir = workdir + 'std_rots/'
os.makedirs(outdir, exist_ok = True)

addhelix_tetrahedral_sc_rots(workdir, outdir, 'helix_std_rots_vdmclu0/', 'AAA_cluster_0_sc_0.64.pdb')
addhelix_tetrahedral_sc_rots(workdir, outdir, 'helix_std_rots_vdmclu1/', 'AAA_cluster_83_sc_0.58.pdb')


### Generate C2 by rotating 180 at z. 

def generate_rvdm_C2(workdir, outdir):
    for file in os.listdir(workdir):
        if not '.pdb' in file or not 'A' in file:
            continue
        rvdmH = pr.parsePDB(workdir + file)
        ag = rot_rvdmH_C2(rvdmH)
        pr.writePDB(outdir + ag.getTitle() +'.pdb', ag)
    return

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2/'
indir = workdir + '_z_fix_oriZN/'
outdir = workdir + '_z_fix_rot_oriZN/'
os.makedirs(outdir, exist_ok=True)
generate_rvdm_C2(indir, outdir)

'''
workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2/helix_std_rots_vdmclu0/'
outdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2/rvdm_clu0_c2/'
os.makedirs(outdir, exist_ok=True)
generate_rvdm_C2(workdir, outdir)
'''

### generate 2vdm combinations. (This can be run on gpu)

def generate_2vdm_combination_in_lineargeo(helix_noclash_outdir, vh_dir1, vh_dir2):
    A_s = []
    B_s = []
    for x in os.listdir(vh_dir1):
        if not '.pdb' in x:
            continue
        if 'A' in x:
            A_s.append(pr.parsePDB(vh_dir1 + x))
    for x in os.listdir(vh_dir2):
        if not '.pdb' in x:
            continue
        if 'B' in x:
            B_s.append(pr.parsePDB(vh_dir2 + x))

    os.makedirs(helix_noclash_outdir, exist_ok=True)

    for i in range(120):
        A = A_s[i].select('protein and heavy')
        B = B_s[0].select('protein and heavy')

        if clash(A.getCoords(), B.getCoords()):
            continue

        title = A.getTitle() + '_' + B.getTitle() 
        AB = pr.AtomGroup(title)
        aa = A_s[i].select('resnum 5 6 7 8 9 10 11').toAtomGroup() 
        ba = B_s[0].select('resnum 5 6 7 8 9 10 11').toAtomGroup()
        ba.setChids(['B' for i in range(len(ba))])

        AB = aa + ba
        
        pr.writePDB(helix_noclash_outdir + title + '.pdb', AB) 
    return


workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2/'
#workdir = '/home/gpu/Lei/DesignData/_reverse_design/'
vh_dir1 = workdir + 'helix_std_rots_vdmclu0/'
vh_dir2 = workdir + 'helix_std_rots_vdmclu1/'
generate_2vdm_combination_in_lineargeo(workdir + 'helix_noclash_0-0/', vh_dir1, vh_dir1)
generate_2vdm_combination_in_lineargeo(workdir + 'helix_noclash_0-1/', vh_dir1, vh_dir2)
generate_2vdm_combination_in_lineargeo(workdir + 'helix_noclash_1-1/', vh_dir2, vh_dir2)


