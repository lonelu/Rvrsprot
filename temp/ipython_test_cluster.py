import os
import sys
sys.path.append(r'/mnt/e/GitHub_Design/Qbits')
sys.path.append(r'/mnt/e/GitHub_Design/smallprot')

import prody as pr
import numpy as np
from scipy.sparse import csr_matrix
from smallprot import cluster_loops

#path = '/mnt/e/GitHub_Design/smallprot/huong/output1/loops_C_D/'
path = '/mnt/e/DesignData/smallprot_loops/rocker/test_cluster/'

dirs = [f for f in os.listdir(path) if f[0] != '.']

loopdir = dirs[0]
loopsize = int(loopdir)
pdbs = list()

for f in [f for f in os.listdir(path + loopdir) 
            if f[-3:] == 'pdb']:
    pdb = pr.parsePDB(path + loopdir + '/' + f)
    ca_sel = pdb.select('name CA')
    if len(ca_sel) == 14 + loopsize:
        pdbs.append(pdb)

kwargs = dict(rmsd_cutoff=1.0, 
                cluster_tag='_loop_' + str(loopsize),
                clusters_outdir=path + loopdir + '/clusters/',
                selection='name CA')

clu = cluster_loops.Cluster(**kwargs)

clu.set_pdbs(pdbs)

clu.set_coords()
'''
pdb = pdbs[0]
coords = pdb.select(clu.selection).getCoords()
type(coords) #numpy.ndarray
clu.pdb_coords.append(coords)
type(clu.pdb_coords) #list

clu.pdb_coords = list()
for pdb in clu.pdbs:
    clu._set_coords(pdb)
type(clu.pdb_coords) #list   

clu.pdb_coords = np.array(clu.pdb_coords, dtype='float32')
type(clu.pdb_coords) #numpy.ndarray
'''

#-----------------------------------------
X= clu.pdb_coords
M = X.shape[0]
N = X.shape[1]
O = X.shape[2]
D = np.zeros((M, M), dtype=np.float32)
m_com = np.zeros(O, dtype=np.float32)
t_com = np.zeros(O, dtype=np.float32)
m = np.zeros((N, O), dtype=np.float32)
mtrans = np.zeros((O, N), dtype=np.float32)
mtr = np.zeros((N, O), dtype=np.float32)
t = np.zeros((N, O), dtype=np.float32)
c = np.zeros((O, O), dtype=np.float32)
U = np.zeros((O, O), dtype=np.float32)
S = np.zeros(O, dtype=np.float32)
Wt = np.zeros((O, O), dtype=np.float32)
R = np.zeros((O, O), dtype=np.float32)
mtr_re = np.zeros(N * O, dtype=np.float32)
t_re = np.zeros(N * O, dtype=np.float32)
sub = np.zeros(N * O, dtype=np.float32)
for i in range(M):
    for j in range(i + 1, M):
        for k in range(O):
            m_com[k] = np.mean(X[i, :, k])
            t_com[k] = np.mean(X[j, :, k])
        m = np.subtract(X[i, :, :], m_com)
        for a in range(N):
            for b in range(O):
                mtrans[b, a] = m[a, b]
        t = np.subtract(X[j, :, :], t_com)
        c = np.dot(mtrans, t)
        U, S, Wt = np.linalg.svd(c)
        R = np.dot(U, Wt)
        if np.linalg.det(R) < 0.0:
            Wt[-1, :] *= -1.0
            R = np.dot(U, Wt)
        mtr = np.add(np.dot(m, R), t_com)
        q = 0
        for a in range(N):
            for b in range(O):
                mtr_re[q] = mtr[a, b]
                t_re[q] = X[j, :, :][a, b]
                q += 1
        sub = np.subtract(mtr_re, t_re)
        D[i, j] = np.sqrt(1.0 / N * np.dot(sub, sub))
#-----------------------------------------
clu.make_pairwise_rmsd_mat()

clu.rmsd_mat[0]

if not clu._square:
    print("make_square")
    clu.make_square()

if not clu._adj_mat:
    print("make_adj_mat")
    clu.make_adj_mat()

clu.rmsd_mat[0]

clu.adj_mat[0]

adj_mat = clu.adj_mat


#test greedy
all_mems = []
cents = []
indices = np.arange(adj_mat.shape[0])

while adj_mat.shape[0] > 0:
    cent = adj_mat.sum(axis=1).argmax()
    cents.append(indices[cent])
    row = adj_mat.getrow(cent)
    tf = ~row.toarray().astype(bool)[0]
    mems = indices[~tf]
    all_mems.append(mems)
    indices = indices[tf]
    adj_mat = adj_mat[tf][:, tf]

#test print_cluster
clu.mems, clu.cents = clu.greedy(clu.adj_mat)

clusters=range(1, len(clu.mems) + 1)
cluster_number = list(clusters)[0]

outpath=None

if not outpath:
    outpath = clu.clusters_outdir + '/' + str(cluster_number) + '/'

cluster_index = cluster_number - 1
cent = clu.cents[cluster_index]
mems = clu.mems[cluster_index]

R, m_com, t_com = cluster_loops.get_rot_trans(clu.pdb_coords[cent],
                                clu.pdb_coords[clu.cents[0]])
cent_coords = np.dot((clu.pdb_coords[cent] - m_com), R) + t_com

#----------------------------test get_rot_trans
mob_coords = clu.pdb_coords[0]
targ_coords = clu.pdb_coords[2]

mob_coords_com = mob_coords.mean(0)
targ_coords_com = targ_coords.mean(0)
mob_coords_cen = mob_coords - mob_coords_com
targ_coords_cen = targ_coords - targ_coords_com
cov_matrix = np.dot(mob_coords_cen.T, targ_coords_cen)
U, S, Wt = np.linalg.svd(cov_matrix)
R = np.dot(U, Wt)
if np.linalg.det(R) < 0.:
    Wt[-1] *= -1
    R = np.dot(U, Wt)

#----------------------------------------------


for i, mem in enumerate(mems):
    print(i, '-', mem)

i = 0
mem = 0
R, m_com, t_com = cluster_loops.get_rot_trans(clu.pdb_coords[mem], cent_coords)
pdb = clu.pdbs_errorfree[mem].copy()
pdb_coords = pdb.getCoords()
coords_transformed = np.dot((pdb_coords - m_com), R) + t_com
pdb.setCoords(coords_transformed)
is_cent = '_centroid' if mem == cent else ''

pr.writePDB(outpath + 'cluster_' + str(cluster_number) + '_mem_' + str(i)
            + clu.cluster_tag + is_cent + '_' + str(pdb).split()[-1] + '_test.pdb.gz', pdb)
