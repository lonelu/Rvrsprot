import numpy as np
from scipy.spatial.distance import cdist
from smallprot import pdbutils, query, cluster_loops, smallprot_config, logger

backbone = ['N', 'CA', 'C', 'O']

def cal_cdist(sse_list):
    qrep_xyz = []
    qrep_natoms = [0]
    qrep_nres = []

    for j, pair_qrep in enumerate(sse_list):
        title = 'struct{}'.format(str(j))
        this_struct = pdbutils.get_struct(title, pair_qrep, 1)
        atoms = this_struct.get_atoms()
        qrep_xyz.append(np.array([atom.get_coord() for atom in atoms if atom.get_name() in backbone]))            
        qrep_nres.append(int(len(qrep_xyz[-1])/4))
        qrep_natoms.append(qrep_natoms[-1] + len(qrep_xyz[-1]))

    qrep_xyz = np.vstack(qrep_xyz)
    # compute interatomic distance between SSE pairs
    dists = cdist(qrep_xyz, qrep_xyz)
    return qrep_natoms, qrep_nres, dists

#Calculate minimun distance between two sses.
#sse_list = [a, b], where a <= b.
def cal_sse_dist(sse_list):
    n_reps = len(sse_list)
    if n_reps != 2:
        raise AssertionError('sse_list does not contain 2 sses.')

    qrep_natoms, qrep_nres, dists = cal_cdist(sse_list)

    #calcualte distances between two the first sse and another sse, 
    #store the minimum index.
    win_dists = np.zeros(qrep_nres[1] - qrep_nres[0] + 1)
    for y in range(qrep_nres[1] - qrep_nres[0] +1):      
        for z in range(qrep_nres[0]*4):
            win_dists[y] += dists[qrep_natoms[0] + z, qrep_natoms[1] + y*4 + z]
    min_dist_ind = np.argmin(win_dists)
    min_dist = win_dists[min_dist_ind]/(qrep_nres[0]*4)

    return min_dist, min_dist_ind

def cal_aa_dist(sse_list):
    n_reps = len(sse_list)
    if n_reps != 2:
        raise AssertionError('sse_list does not contain 2 sses.')

    qrep_natoms, qrep_nres, dists = cal_cdist(sse_list)

    # if qrep_nres[0] != qrep_nres[1]:
    #     raise AssertionError('sse_list 2 sses does not have equal length.')
    win_dists = np.zeros((qrep_nres[0], qrep_nres[1]))
    for i in range(qrep_nres[0]):
        for j in range(qrep_nres[1]):
            for z in range(4):
                win_dists[i, j] += dists[qrep_natoms[0] + i*4 + z, qrep_natoms[1] + j*4 + z]      
    min_dist_ind = np.unravel_index(win_dists.argmin(), win_dists.shape)
    min_dist = win_dists[min_dist_ind]
    return min_dist, min_dist_ind, qrep_nres

def cal_loop_mid_dist(sse_list):
    n_reps = len(sse_list)
    if n_reps != 2:
        raise AssertionError('sse_list does not contain 2 sses.')

    qrep_natoms, qrep_nres, dists = cal_cdist(sse_list)

    dist = 0
    for z in range(4):
        dist += dists[qrep_natoms[0] + int(qrep_nres[0]/2)*4 + z, qrep_natoms[1] + int(qrep_nres[1]/2)*4 + z]
    return dist   


