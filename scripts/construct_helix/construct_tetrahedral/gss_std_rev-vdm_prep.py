'''
Generate helix bundles with the best zn vdM from helix bundle. 
The vdMs are inversed for helix bundle master searching.
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
import shutil

def clash(Acoords, Bcoords, vdm_vdm_clash_dist = 2.6):
    neigh_y = NearestNeighbors(radius= vdm_vdm_clash_dist)
    neigh_y.fit(Acoords)
    x_in_y = neigh_y.radius_neighbors(Bcoords)
    x_has_y = any([True if len(a) >0 else False for a in x_in_y[1]])
    if x_has_y:
        return True
    return False 


### Generate vdm rotate with different angles.

def _generate_tetrahedral_sc_rots(outdir, std, std_cp_sel, tf_sel, tag_pre):
    dist = pr.calcDistance(std.select('serial 1'), std.select('serial 2'))[0]
    xcoord = np.array([[0, 0, 0], [dist , 0, 0]])
    tf = pr.calcTransformation(std.select(tf_sel).getCoords(), xcoord)
    tf_rv = pr.calcTransformation(xcoord, std.select(tf_sel).getCoords())
    
    for i in range(120):
        std_cp = std.select(std_cp_sel).copy()
        tf.apply(std_cp)
        v = std_cp.getCoords()
        theta = 2*math.pi/120 * i
        axis = xcoord[1]/norm(xcoord[1])
        rot = Rotation.from_rotvec(theta * axis)
        new_v = rot.apply(v) 
        new_v = tf_rv.apply(new_v) 
        std_cp.setCoords(new_v)
        pr.writePDB(outdir + tag_pre + str(i) +'.pdb', std_cp)


def generate_tetrahedral_sc_rots(std, outdir):
    '''
    It turned out the original function has a bug. 
    The rotation of the molecula will change the molecula as the rotation is change the shape of the 3d object.
    One need to transform the molecula and align along the axis of the rotation vector.

    # Testing rotation
    std_cp = std.copy()
    v = std_cp.select('serial 3 4 5 6 7').getCoords()
    axis = std_cp.select('serial 1').getCoords() - std_cp.select('serial 3').getCoords()
    theta = 1.57

    axis = axis / norm(axis)
    rot = Rotation.from_rotvec(theta * axis)
    new_v = rot.apply(v)  

    std_cp.select('serial 3 4 5 6 7').setCoords(new_v)
    pr.writePDB(workdir + 'test.pdb', std_cp)
    '''
    _generate_tetrahedral_sc_rots(outdir, std.copy(), 'serial 2 3 4 5 6', 'serial 1 2', 'A_')

    _generate_tetrahedral_sc_rots(outdir, std.copy(), 'serial 7 8 9 10 11', 'serial 1 7', 'B_')

    _generate_tetrahedral_sc_rots(outdir, std.copy(), 'serial 12 13 14 15 16', 'serial 1 12', 'C_')

    _generate_tetrahedral_sc_rots(outdir, std.copy(), 'serial 17 18 19 20 21', 'serial 1 17', 'D_')

    return


workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/'
outdir = workdir + 'std_rots_fix/'
os.makedirs(outdir, exist_ok = True)

std = pr.parsePDB(workdir + 'stdt.pdb')
#std.getSerials()
generate_tetrahedral_sc_rots(std, outdir)


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

addhelix_tetrahedral_sc_rots(workdir, outdir, 'helix_std_rots_vdmclu0/', 'AAA_cluster_0_sc_0.64.pdb')
addhelix_tetrahedral_sc_rots(workdir, outdir, 'helix_std_rots_vdmclu1/', 'AAA_cluster_83_sc_0.58.pdb')


### Create C3 bundles. 

def generate_c3_by_rotation(outdir, vdm_helix_dir, std):
    '''
    Select the rotation A. Rotate A 120 degree to create C3 helix. 
    '''

    zn_coord = std.select('serial 1').getCoords()[0]
    cent = pr.calcCenter(std.select('serial 2 7 12'))
    dist = pr.calcDistance(zn_coord, cent)

    origin_coord = np.array((zn_coord, cent))
    xcoord = np.array([[0, 0, 0], [dist , 0, 0]])

    tf = pr.calcTransformation(origin_coord, xcoord)
    tf_rv = pr.calcTransformation(xcoord, origin_coord)

    #imax = 0
    v_hs = []
    for v_h in os.listdir(vdm_helix_dir):
        # imax+=1
        # if imax >= 3:
        #     continue
        if not '.pdb' in v_h or not 'A' in v_h:
            continue
        v_hs.append(tf.apply(pr.parsePDB(vdm_helix_dir + v_h)))
        
        
    for v_h in v_hs:
        metal = v_h.select('name ZN').copy().toAtomGroup()
        v_h_a = v_h.select('protein').copy().toAtomGroup()
        v_h_b = v_h.select('protein').copy().toAtomGroup()
        v_h_b.setChids('B')
        v_h_c = v_h.select('protein').copy().toAtomGroup()
        v_h_c.setChids('C')
        #tf_rv.apply(v_h_a) 

        v = v_h_a.getCoords()
        theta = 2*math.pi/3*1
        axis = xcoord[1]/norm(xcoord[1])
        rot = Rotation.from_rotvec(theta * axis)
        new_v = rot.apply(v) 
        #new_v = tf_rv.apply(new_v) 
        v_h_b.setCoords(new_v)

        #if clash(v_h_a.getCoords(), v_h_b.getCoords(), 2.3):
        #    continue

        v = v_h_a.getCoords()
        theta = 2*math.pi/3*2
        axis = xcoord[1]/norm(xcoord[1])
        rot = Rotation.from_rotvec(theta * axis)
        new_v = rot.apply(v) 
        #new_v = tf_rv.apply(new_v) 
        v_h_c.setCoords(new_v)

        v_all = v_h_a + v_h_b + v_h_c + metal

        pr.writePDB(outdir + v_h.getTitle() +'_C3.pdb', v_all)

    return

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/'
outdir = workdir + 'c3_vc0/'
os.makedirs(outdir, exist_ok = True)

generate_c3_by_rotation(outdir, workdir + 'helix_std_rots_vdmclu0/', std)

outdir = workdir + 'c3_vc1/'
os.makedirs(outdir, exist_ok = True)
generate_c3_by_rotation(outdir, workdir + 'helix_std_rots_vdmclu1/', std)

### generate 3vdm combinations. (This can be run on gpu)

def generate_3vdm_combination_in_tetrahedralgeo(workdir, helix_noclash_outdir, vh_dir1, vh_dir2, vh_dir3):
    A_s = []
    B_s = []
    C_s = []
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
    for x in os.listdir(vh_dir3):
        if not '.pdb' in x:
            continue
        if 'C' in x:
            C_s.append(pr.parsePDB(vh_dir3 + x))   

    os.makedirs(helix_noclash_outdir, exist_ok=True)

    for i, j, k in permutations(range(120), 3):
        A = A_s[i].select('protein and heavy')
        B = B_s[j].select('protein and heavy')
        C = C_s[k].select('protein and heavy')

        if clash(A.getCoords(), B.getCoords()):
            continue
        if clash(A.getCoords(), C.getCoords()):
            continue
        if clash(B.getCoords(), C.getCoords()):
            continue
        title = A.getTitle() + '_' + B.getTitle() + '_' + C.getTitle()
        ABC = pr.AtomGroup(title)
        aa = A_s[i].select('resnum 5 6 7 8 9 10 11').toAtomGroup() 
        ba = B_s[j].select('resnum 5 6 7 8 9 10 11').toAtomGroup()
        ba.setChids(['B' for i in range(len(ba))])
        ca = C_s[k].select('resnum 5 6 7 8 9 10 11').toAtomGroup() 
        ca.setChids(['C' for i in range(len(ca))])
        ABC = aa + ba + ca
        
        pr.writePDB(helix_noclash_outdir + title + '.pdb', ABC) 
    return

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/'
#workdir = '/home/gpu/Lei/DesignData/_reverse_design/'
vh_dir1 = workdir + 'helix_std_rots_vdmclu0/'
vh_dir2 = workdir + 'helix_std_rots_vdmclu1/'
generate_3vdm_combination_in_tetrahedralgeo(workdir, workdir + 'helix_noclash_0-0-0/', vh_dir1, vh_dir1, vh_dir1)
generate_3vdm_combination_in_tetrahedralgeo(workdir, workdir + 'helix_noclash_0-0-1/', vh_dir1, vh_dir1, vh_dir2)
generate_3vdm_combination_in_tetrahedralgeo(workdir, workdir + 'helix_noclash_0-1-1/', vh_dir1, vh_dir2, vh_dir2)
generate_3vdm_combination_in_tetrahedralgeo(workdir, workdir + 'helix_noclash_1-1-1/', vh_dir2, vh_dir2, vh_dir2)

### cluster 3vdm combinations. (this can be run on gpu.)

def cluster_and_extract_centroid(workdir_clu, outdir_clu, outdir_cent, pretag, rmsd = 0.75, len_sel = 21, align_sel = 'name CA'):
    '''
    We don't want to master search every 3vdm or 2vdm, 
    so we cluster them first,
    then extract the centroid and we only search the centroid. 
    '''
    _pdbs = [pr.parsePDB(workdir_clu + x) for x in os.listdir(workdir_clu)]

    clu = database_cluster.superimpose_aa_core(_pdbs, rmsd, len_sel, align_sel, min_cluster_size = 0)

    database_cluster.print_cluster_pdbs(clu, outdir_clu, rmsd, tag = pretag)

    os.makedirs(outdir_cent, exist_ok=True)

    for x in os.listdir(outdir_clu):
        if not os.path.isdir(outdir_clu + x):
            continue
        for y in os.listdir(outdir_clu + x + '/'):
            if 'centroid' not in y:
                continue
            shutil.copy(outdir_clu + x + '/' + y, outdir_cent + y)
    return 

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/'

cluster_and_extract_centroid(workdir + 'helix_noclash_0-0-0/', workdir + 'helix_noclash_0-0-0_clu/', workdir + 'helix_noclash_0-0-0_cent/', 'v0-0-0_', rmsd = 0.75, len_sel = 21, align_sel = 'name CA')
cluster_and_extract_centroid(workdir + 'helix_noclash_0-0-1/', workdir + 'helix_noclash_0-0-1_clu/', workdir + 'helix_noclash_0-0-1_cent/', 'v0-0-1_', rmsd = 0.75, len_sel = 21, align_sel = 'name CA')
cluster_and_extract_centroid(workdir + 'helix_noclash_0-1-1/', workdir + 'helix_noclash_0-1-1_clu/', workdir + 'helix_noclash_0-1-1_cent/', 'v0-1-1_', rmsd = 0.75, len_sel = 21, align_sel = 'name CA')
cluster_and_extract_centroid(workdir + 'helix_noclash_1-1-1/', workdir + 'helix_noclash_1-1-1_clu/', workdir + 'helix_noclash_1-1-1_cent/', 'v1-1-1_', rmsd = 0.75, len_sel = 21, align_sel = 'name CA')

### generate 2vdm combinations. 

def generate_2vdm_combination_in_tetrahedralgeo(workdir, outdir_name, helix_std_dir1, helix_std_dir2):
    A_s = []
    B_s = []
    for x in os.listdir(helix_std_dir1):
        if not '.pdb' in x:
            continue
        if 'A' in x:
            A_s.append(pr.parsePDB(helix_std_dir1 + x))
    for x in os.listdir(helix_std_dir2):
        if not '.pdb' in x:
            continue
        if 'B' in x:
            B_s.append(pr.parsePDB(helix_std_dir2 + x))
   
    helix_noclash_outdir = workdir + outdir_name
    os.makedirs(helix_noclash_outdir, exist_ok=True)

    for i, j in permutations(range(120), 2):
        A = A_s[i].select('heavy')
        B = B_s[j].select('protein and heavy')

        if clash(A.select('protein').getCoords(), B.getCoords()):
            continue

        title = 'A_'  + str(i) + '_B_' + str(j)
        #title = A.getTitle() + '_' + B.getTitle()
        ABC = pr.AtomGroup(title)
        aa = A_s[i].select('resnum 5 6 7 8 9 10 11 or name ZN').toAtomGroup() 
        ba = B_s[j].select('resnum 5 6 7 8 9 10 11').toAtomGroup()
        ba.setChids(['B' for i in range(len(ba))])

        ABC = aa + ba
        
        pr.writePDB(helix_noclash_outdir + title + '.pdb', ABC) 
    return

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c3/'

generate_2vdm_combination_in_tetrahedralgeo(workdir, 'c3_pair/helix_noclash_0-0_vdm/', helix_std_dir1 = workdir + 'helix_std_rots_vdmclu0/', helix_std_dir2 = workdir + 'helix_std_rots_vdmclu0/')

generate_2vdm_combination_in_tetrahedralgeo(workdir, 'c3_pair/helix_noclash_0-1_vdm/', helix_std_dir1 = workdir + 'helix_std_rots_vdmclu0/', helix_std_dir2 = workdir + 'helix_std_rots_vdmclu1/')

generate_2vdm_combination_in_tetrahedralgeo(workdir, 'c3_pair/helix_noclash_1-1_vdm/', helix_std_dir1 = workdir + 'helix_std_rots_vdmclu1/', helix_std_dir2 = workdir + 'helix_std_rots_vdmclu1/')

### Cluster the 2 vdm_combination and use the centroid for master search.


workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c3/c3_pair/'

cluster_and_extract_centroid(workdir + 'helix_noclash_0-0_vdm/', workdir + 'helix_noclash_0-0_vdm_clu/', workdir + 'helix_noclash_0-0_vdm_cent/', 'v0-0_', rmsd = 0.75, len_sel = 14, align_sel = 'name CA')
cluster_and_extract_centroid(workdir + 'helix_noclash_0-1_vdm/', workdir + 'helix_noclash_0-1_vdm_clu/', workdir + 'helix_noclash_0-1_vdm_cent/', 'v0-1_', rmsd = 0.75, len_sel = 14, align_sel = 'name CA')
cluster_and_extract_centroid(workdir + 'helix_noclash_1-1_vdm/', workdir + 'helix_noclash_1-1_vdm_clu/', workdir + 'helix_noclash_1-1_vdm_cent/', 'v1-1_', rmsd = 0.75, len_sel = 14, align_sel = 'name CA')


### using 2vdm_combination to generate 3vdm_combination

def generate_3vdm_combination_with2vdm(workdir, outdir_name, file_name, helix_std_dir1):
    protein = pr.parsePDB(workdir + file_name)
    C_s = []
    for x in os.listdir(helix_std_dir1):
        if not '.pdb' in x:
            continue
        if 'C' in x:
            C_s.append(pr.parsePDB(helix_std_dir1 + x))
    print('there are {} candidates of C.'.format(len(C_s)))
    helix_noclash_outdir = workdir + outdir_name
    os.makedirs(helix_noclash_outdir, exist_ok=True)

    for i in range(120):
        C = C_s[i]
        if clash(protein.select('protein and heavy').getCoords(), C.select('protein and heavy').getCoords()):
            continue
        title = file_name + '_' + C.getTitle()
        #title = A.getTitle() + '_' + B.getTitle()
        ABC = pr.AtomGroup(title)
        aa = protein.toAtomGroup() 
        ca = C.copy().select('resnum 5 6 7 8 9 10 11').toAtomGroup() 
        ca.setChids(['C' for i in range(len(ca))])
        ABC = aa + ca
        
        pr.writePDB(helix_noclash_outdir + title + '.pdb', ABC) 
    return    

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/helix_noclash_1-1_vdm_search/'

file_name = 'A_42_B_44.pdb'
generate_3vdm_combination_with2vdm(workdir, 'A_42_B44_1-1-0_noclash/', file_name, '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/helix_std_rots/')
generate_3vdm_combination_with2vdm(workdir, 'A_42_B44_1-1-1_noclash/', file_name, '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/helix_std_rots_vdmclu1/')

### Statistics on gpu
import os
import shutil
import statistics

workdir = '/home/gpu/Lei/DesignData/_reverse_design/'
helix_noclash_outdir = workdir + 'helix_noclash_1-1_vdm/'
helix_noclash_masterfilter_outdir = workdir + 'helix_noclash_1-1_vdm_master/'

_summary = []
for folder in os.listdir(helix_noclash_masterfilter_outdir):
    _dir = helix_noclash_masterfilter_outdir + folder + '/'
    if not os.path.isdir(_dir):
        continue
    
    rmsds = []
    with open(_dir + 'seq.txt') as f:
        lines = f.readlines()
        if len(lines) <= 0:
            continue
        for line in lines:
            rmsds.append(float(line.strip(' ').split(' ')[0]))
    _summary.append((folder, statistics.mean(rmsds), min(rmsds), len(lines)))
    #shutil.copy(helix_noclash_outdir + folder + '.pdb', _dir + folder + '.pdb')

with open(workdir + '_summary_1-1.txt', 'w') as f:
    f.write('name\tmean_rmsd\tmin_rmsd\tcount\n')
    for s in _summary:
        f.write('\t'.join([str(x) for x in s]) + '\n')


### select and copy top hits into a new folder.

helix_noclash_sel_outdir = workdir + 'helix_noclash_sel/'
os.makedirs(helix_noclash_sel_outdir)
_summary_sort = sorted(_summary, key=lambda tup: tup[3], reverse=True)

for s in _summary_sort[0:74]:
    folder = s[0]
    shutil.copy(helix_noclash_outdir + folder + '.pdb', helix_noclash_sel_outdir + folder + '.pdb')


for s in _summary_sort:
    folder = s[0]
    if s[3] >= 3:
        shutil.copy(helix_noclash_outdir + folder + '.pdb', helix_noclash_sel_outdir + folder + '.pdb')