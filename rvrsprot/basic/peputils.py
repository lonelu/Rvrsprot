import numpy as np
from scipy.spatial.distance import cdist
from . import pdbutils, cluster_loops, logger
from ..external import query
import prody as pr

backbone = ['N', 'CA', 'C', 'O']

def cal_cdist(sse_list):
    '''
    compute interatomic distance between SSE pairs
    '''
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
    dists = cdist(qrep_xyz, qrep_xyz)
    return qrep_natoms, qrep_nres, dists

def cal_sse_dist(sse_list):
    '''
    Calculate minimun distance between two sses.
    sse_list = [a, b], where a <= b.
    '''
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

def cal_aa_dist(sse_list, seedfirst, loop_query_win, keep =0):
    '''
    Calculate minimum distance between [seed and loop] or [loop and seed]
    return index of the two aa postions of miniumn distant.
    '''
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
    if keep ==0:
        min_dist_ind = np.unravel_index(win_dists.argmin(), win_dists.shape)
    elif keep == -1 and seedfirst:
        ind = win_dists[:, 0].argmin()
        min_dist_ind = [ind, 0]
    elif keep == -1 and not seedfirst:
        ind = win_dists[qrep_nres[0]-1, :].argmin()
        min_dist_ind = [qrep_nres[0]-1, ind] 
    elif keep == 1 and seedfirst:   
        ind = win_dists[:, loop_query_win -1].argmin()
        min_dist_ind = [ind, loop_query_win -1]    
    elif keep == 1 and not seedfirst:
        ind = win_dists[qrep_nres[0] - loop_query_win, :].argmin()
        min_dist_ind = [qrep_nres[0] - loop_query_win, ind]
        
    min_dist = win_dists[min_dist_ind]
    
    return min_dist, min_dist_ind, qrep_nres

#Calculate minimum distance between [one aa, seq]
# return index of the aa and the aa positions of the qeq.

def cal_loop_mid_dist(sse_list):
    '''
    If two loops are at the same directions. Their distance should be within a limitation.
    For example, two loops at the same direction of a 4 helix bundles. The distance should be smaller than the 2*radius (2*R0) of the 4 helix bundles.
    '''
    n_reps = len(sse_list)
    if n_reps != 2:
        raise AssertionError('sse_list does not contain 2 sses.')

    qrep_natoms, qrep_nres, dists = cal_cdist(sse_list)

    dist = 0
    for z in range(4):
        dist += dists[qrep_natoms[0] + int(qrep_nres[0]/2)*4 + z, qrep_natoms[1] + int(qrep_nres[1]/2)*4 + z]
    return dist   

def cal_main_chian_dist(sse1_path, sse2_path):
    '''
    Use prody to get the median of min distance of each atom in sse1 to an atom in sse2.
    For a 'good' helix bundles, each helix have two neibors that are of similar distances. 
    Since there is no guaranttee of the atom selected of the second sse. This calculation is rough.
    TO THINK: There may have a better way to do this.
    '''
    sse1 = pr.parsePDB(sse1_path)
    atoms1 = sse1.select('name N CA C O')
    sse2 = pr.parsePDB(sse2_path)
    atoms2 = sse2.select('name N CA C O')
    
    #dists = pr.calcDistance(atoms1, atoms2) 
    distmat = pr.buildDistMatrix(atoms1, atoms2[:-1])
    median_ave_dist = np.median(np.amin(distmat, axis = 1))
    return median_ave_dist

def get_neighbor_dists(full_sse_list, direction):
    dists = []
    for i in range(len(direction)): 
        dist1 = cal_main_chian_dist(full_sse_list[i], full_sse_list[i-1])
        ind2 = i+1 if i < len(direction)-1 else 0
        dist2 = cal_main_chian_dist(full_sse_list[i], full_sse_list[ind2])
        dists.append(dist1)
        dists.append(dist2)
    return dists



    





