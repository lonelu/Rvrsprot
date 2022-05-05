import os
import sys
import prody as pr
from scipy.spatial.transform import Rotation
from numpy.linalg import norm
import math
from itertools import permutations
import numpy as np
from sklearn.neighbors import NearestNeighbors
from metalprot.database import database_cluster

import shlex
import subprocess
import multiprocessing as mp
from distutils.dir_util import copy_tree
import shutil

def clash(Acoords, Bcoords, vdm_vdm_clash_dist = 2.6):
    neigh_y = NearestNeighbors(radius= vdm_vdm_clash_dist)
    neigh_y.fit(Acoords)
    x_in_y = neigh_y.radius_neighbors(Bcoords)
    x_has_y = any([True if len(a) >0 else False for a in x_in_y[1]])
    if x_has_y:
        return True
    return False 

### generate 3vdm combinations. (This can be run on gpu)

def generate_3vdm_combination_in_tetrahedralgeo(workdir, helix_noclash_outdir, vh_dir1, vh_dir2, vh_dir3):
    print('run: ' + helix_noclash_outdir)
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


### cluster 3vdm combinations. (this can be run on gpu.)

def cluster_centroid(workdir_clu, outdir_clu, outdir_cent, pretag, rmsd = 0.5, len_sel = 21, align_sel = 'name CA'):
    '''
    We don't want to master search every 3vdm or 2vdm, 
    so we cluster them first,
    then extract the centroid and we only search the centroid. 
    '''
    print('run: ' + outdir_clu)
    _pdbs = [pr.parsePDB(workdir_clu + x) for x in os.listdir(workdir_clu)]

    clu = database_cluster.superimpose_aa_core(_pdbs, rmsd, len_sel, align_sel, min_cluster_size = 0)

    database_cluster.print_cluster_pdbs(clu, outdir_clu, rmsd, tag = pretag)

    return 


def extract_centroid(outdir_clu, outdir_cent):
    os.makedirs(outdir_cent, exist_ok = True)
    for x in os.listdir(outdir_clu):
        if not os.path.isdir(outdir_clu + x):
            continue
        for y in os.listdir(outdir_clu + x):
            if 'centroid' not in y:
                continue
            shutil.copy(outdir_clu + x + '/' + y, outdir_cent + y)
    return


def main_wynton():
    ind = int(sys.argv[1]) -1
    workdir = '/wynton/home/degradolab/lonelu/DesignData/_reverse_design/'

    vh_dir1 = workdir + 'helix_std_rots_vdmclu0/'
    vh_dir2 = workdir + 'helix_std_rots_vdmclu1/'
    parameters1 = [
        [workdir, workdir + 'helix_noclash_0-0-0/', vh_dir1, vh_dir1, vh_dir1], 
        [workdir, workdir + 'helix_noclash_0-0-1/', vh_dir1, vh_dir1, vh_dir2],
        [workdir, workdir + 'helix_noclash_0-1-1/', vh_dir1, vh_dir2, vh_dir2],
        [workdir, workdir + 'helix_noclash_1-1-1/', vh_dir2, vh_dir2, vh_dir2]
    ]

    p1 = parameters1[ind]
    generate_3vdm_combination_in_tetrahedralgeo(p1[0], p1[1], p1[2], p1[3], p1[4])

    parameters2 = [
        [workdir + 'helix_noclash_0-0-0/', workdir + 'helix_noclash_0-0-0_clu/', workdir + 'helix_noclash_0-0-0_cent/', 'v0-0-0_'], 
        [workdir + 'helix_noclash_0-0-1/', workdir + 'helix_noclash_0-0-1_clu/', workdir + 'helix_noclash_0-0-1_cent/', 'v0-0-1_'],
        [workdir + 'helix_noclash_0-1-1/', workdir + 'helix_noclash_0-1-1_clu/', workdir + 'helix_noclash_0-1-1_cent/', 'v0-1-1_'],
        [workdir + 'helix_noclash_1-1-1/', workdir + 'helix_noclash_1-1-1_clu/', workdir + 'helix_noclash_1-1-1_cent/', 'v1-1-1_']
    ]
    p2 = parameters2[ind]
    cluster_centroid(p2[0], p2[1], p2[2], p2[3])
    extract_centroid(p2[1], p2[2])
    return


def main_local():
    workdir = '/home/gpu/Lei/DesignData/_reverse_design/'
    cluster_centroid(workdir + 'helix_noclash_0-0-1/', workdir + 'helix_noclash_0-0-1_clu/', workdir + 'helix_noclash_0-0-1_cent/', 'v0-0-1_')
    extract_centroid(workdir + 'helix_noclash_0-0-1_clu/', workdir + 'helix_noclash_0-0-1_cent/')
    return

def main_gpu():
    #workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/'
    workdir = '/home/gpu/Lei/DesignData/_reverse_design/'

    vh_dir1 = workdir + 'helix_std_rots_vdmclu0/'
    vh_dir2 = workdir + 'helix_std_rots_vdmclu1/'
    parameters = [
        [workdir, workdir + 'helix_noclash_0-0-0/', vh_dir1, vh_dir1, vh_dir1], 
        [workdir, workdir + 'helix_noclash_0-0-1/', vh_dir1, vh_dir1, vh_dir2],
        [workdir, workdir + 'helix_noclash_0-1-1/', vh_dir1, vh_dir2, vh_dir2],
        [workdir, workdir + 'helix_noclash_1-1-1/', vh_dir2, vh_dir2, vh_dir2]
    ]

    
    # generate_3vdm_combination_in_tetrahedralgeo(workdir, workdir + 'helix_noclash_0-0-0/', vh_dir1, vh_dir1, vh_dir1)
    # generate_3vdm_combination_in_tetrahedralgeo(workdir, workdir + 'helix_noclash_0-0-1/', vh_dir1, vh_dir1, vh_dir2)
    # generate_3vdm_combination_in_tetrahedralgeo(workdir, workdir + 'helix_noclash_0-1-1/', vh_dir1, vh_dir2, vh_dir2)
    # generate_3vdm_combination_in_tetrahedralgeo(workdir, workdir + 'helix_noclash_1-1-1/', vh_dir2, vh_dir2, vh_dir2)

    pool = mp.Pool(4)
    [pool.apply_async(generate_3vdm_combination_in_tetrahedralgeo, args=(para[0], para[1], para[2], para[3], para[4])) for para in parameters]
    pool.close()
    pool.join()

    parameters = [
        [workdir + 'helix_noclash_0-0-0/', workdir + 'helix_noclash_0-0-0_clu/', workdir + 'helix_noclash_0-0-0_cent/', 'v0-0-0_'], 
        [workdir + 'helix_noclash_0-0-1/', workdir + 'helix_noclash_0-0-1_clu/', workdir + 'helix_noclash_0-0-1_cent/', 'v0-0-1_'],
        [workdir + 'helix_noclash_0-1-1/', workdir + 'helix_noclash_0-1-1_clu/', workdir + 'helix_noclash_0-1-1_cent/', 'v0-1-1_'],
        [workdir + 'helix_noclash_1-1-1/', workdir + 'helix_noclash_1-1-1_clu/', workdir + 'helix_noclash_1-1-1_cent/', 'v1-1-1_']
    ]

    # cluster_and_extract_centroid(workdir + 'helix_noclash_0-0-0/', workdir + 'helix_noclash_0-0-0_clu/', workdir + 'helix_noclash_0-0-0_cent/', 'v0-0-0_', rmsd = 0.75, len_sel = 21, align_sel = 'name CA')
    # cluster_and_extract_centroid(workdir + 'helix_noclash_0-0-1/', workdir + 'helix_noclash_0-0-1_clu/', workdir + 'helix_noclash_0-0-1_cent/', 'v0-0-1_', rmsd = 0.75, len_sel = 21, align_sel = 'name CA')
    # cluster_and_extract_centroid(workdir + 'helix_noclash_0-1-1/', workdir + 'helix_noclash_0-1-1_clu/', workdir + 'helix_noclash_0-1-1_cent/', 'v0-1-1_', rmsd = 0.75, len_sel = 21, align_sel = 'name CA')
    # cluster_and_extract_centroid(workdir + 'helix_noclash_1-1-1/', workdir + 'helix_noclash_1-1-1_clu/', workdir + 'helix_noclash_1-1-1_cent/', 'v1-1-1_', rmsd = 0.75, len_sel = 21, align_sel = 'name CA')
    
    pool = mp.Pool(4)
    [pool.apply_async(cluster_and_extract_centroid, args=(para[0], para[1], para[2], para[3])) for para in parameters]
    pool.close()
    pool.join()

    return

if __name__=='__main__':
    main3()