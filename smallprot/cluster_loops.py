import os
import sys
import prody as pr
import numpy as np
from scipy.sparse import csr_matrix

#The method is kabasch_algorithm
#https://en.wikipedia.org/wiki/Kabsch_algorithm
def get_rot_trans(mob_coords, targ_coords):
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
    return R, mob_coords_com, targ_coords_com

#The method first apply kabasch_algorithm, then calculate rmsd
#https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
def _make_pairwise_rmsd_mat(X):
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
    return D


class Cluster:
    def __init__(self, **kwargs):
        self.rmsd_cutoff = kwargs.get('rmsd_cutoff', 0.5)
        self.cluster_tag = kwargs.get('cluster_tag', '')
        self.selection = kwargs.get('selection', 'all')
        od = kwargs.get('rmsd_mat_outdir', '.')
        self.rmsd_mat_outdir = od if od[-1] == '/' else od + '/'
        od = kwargs.get('clusters_outdir', '.')
        self.clusters_outdir = od if od[-1] == '/' else od + '/'
        self.pdbs = None
        self.pdb_coords = list()
        self.pdbs_errorfree = list()
        self.rmsd_mat = None
        self.adj_mat = None
        self.mems = None
        self.cents = None
        self._square = False
        self._adj_mat = False
        self._data_set = False
        self._ifg_count = list()
        self._vdm_count = list()
        self._query_name = list()
        self._centroid = list()
        self._cluster_num = list()
        self._cluster_size = list()
        self._rmsd_from_centroid = list()
        self._cluster_type = list()
        self.df = None

    def set_pdbs(self, pdbs):
        """Load pdbs into the Cluster object to be clustered.
        VdMs must be all one residue type."""
        self.pdbs = pdbs

    def _set_coords(self, pdb):
        """Grabs coordinates of atoms in bb selection and iFG dict."""
        try:
            coords = pdb.select(self.selection).getCoords()
            self.pdb_coords.append(coords)
            self.pdbs_errorfree.append(pdb)
        except AttributeError:
            pass

    def set_coords(self):
        """Grabs coords for every vdM in pdbs."""
        for pdb in self.pdbs:
            self._set_coords(pdb)
        self.pdb_coords = np.array(self.pdb_coords, dtype='float32')

    def make_pairwise_rmsd_mat(self):
        """Uses C-compiled numba code for fast pairwise superposition
        and RMSD calculation of coords (all against all)."""
        assert isinstance(self.pdb_coords, np.ndarray), 'PDB coords must be ' \
                                                        'numpy array'
        assert self.pdb_coords.dtype == 'float32', 'PDB coords must ' \
                                                   'be dtype of float32'
        self.rmsd_mat = _make_pairwise_rmsd_mat(self.pdb_coords)

    @staticmethod
    def greedy(adj_mat):
        """Takes an adjacency matrix as input.
            All values of adj_mat are 1 or 0:  1 if <= to cutoff, 0 if > cutoff.
            Can generate adj_mat from data in column format with:
            sklearn.neighbors.NearestNeighbors(metric='euclidean',
            radius=cutoff).fit(data).radius_neighbors_graph(data)"""

        if not isinstance(adj_mat, csr_matrix):
            try:
                adj_mat = csr_matrix(adj_mat)
            except:
                print('adj_mat distance matrix must be scipy csr_matrix '
                      '(or able to convert to one)')
                return

        assert adj_mat.shape[0] == adj_mat.shape[1], 'Distance matrix is not square.'

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

        return all_mems, cents

    def make_square(self):
        self.rmsd_mat = self.rmsd_mat.T + self.rmsd_mat
        self._square = True

    def make_adj_mat(self):
        """Makes an adjacency matrix from the RMSD matrix"""
        self.adj_mat = np.zeros(self.rmsd_mat.shape)
        self.adj_mat[self.rmsd_mat <= self.rmsd_cutoff] = 1
        self.adj_mat = csr_matrix(self.adj_mat)
        self._adj_mat = True

    def cluster(self):
        """Performs greedy clustering of the RMSD matrix with a given
        RMSD cutoff (Default cutoff is 0.5 A)."""
        assert self.rmsd_mat is not None, 'Must create rmsd matrix first with ' \
                                          'make_pairwise_rmsd_mat()'
        if not self._square:
            self.make_square()

        if not self._adj_mat:
            self.make_adj_mat()

        self.mems, self.cents = self.greedy(self.adj_mat)

    def print_cluster(self, cluster_number, outpath=None):
        """Prints PDBs of a cluster after superposition of the backbone (bb_sel)
        onto that of the cluster centroid.  The backbone of the cluster centroid
        is itself superposed onto that of the largest cluster's centroid."""

        if not outpath:
            outpath = self.clusters_outdir + '/' + str(cluster_number) + '/'

        try:
            os.makedirs(outpath)
        except FileExistsError:
            pass

        cluster_index = cluster_number - 1
        cent = self.cents[cluster_index]
        mems = self.mems[cluster_index]

        # Align backbone of cluster centroid to backbone of centroid of largest cluster.
        R, m_com, t_com = get_rot_trans(self.pdb_coords[cent],
                                        self.pdb_coords[self.cents[0]])
        cent_coords = np.dot((self.pdb_coords[cent] - m_com), R) + t_com

        for i, mem in enumerate(mems):
            R, m_com, t_com = get_rot_trans(self.pdb_coords[mem], cent_coords)
            pdb = self.pdbs_errorfree[mem].copy()
            pdb_coords = pdb.getCoords()
            coords_transformed = np.dot((pdb_coords - m_com), R) + t_com
            pdb.setCoords(coords_transformed)
            is_cent = '_centroid' if mem == cent else ''
            pr.writePDB(outpath + 'cluster_' + str(cluster_number) + '_mem_' + str(i)
                        + self.cluster_tag + is_cent + '_' + str(pdb).split()[-1] + '.pdb.gz', pdb)

    def print_clusters(self, clusters):
        """Prints all clusters in the list *clusters*"""
        for cluster_number in clusters:
            self.print_cluster(cluster_number)

    def run_protocol(self, pdbs):
        """Runs a series of methods that will cluster the
        PDBs and output a pandas dataframe, RMSD matrix, and
        PDBs of top 20 clusters."""
        print('begin clustering...')
        self.set_pdbs(pdbs)
        print('setting PDB coords...')
        self.set_coords()
        print('constructing RMSD matrix...')
        self.make_pairwise_rmsd_mat()
        print('greedy clustering...')
        self.cluster()
        print('making dataframe...')
        self.print_clusters(clusters=range(1, len(self.mems) + 1))


def listdir_mac(path):
    return [f for f in os.listdir(path) if f[0] != '.']

##########
##########

def run_cluster(path, log, outfile=None):
    if outfile:
        orig_out = sys.stdout
        sys.stdout = open(outfile, 'a')
    for loopdir in listdir_mac(path):
        # if other files exist in the directory, pass over them
        if loopdir not in [str(i) for i in range(1, 100)]:
            continue
        loopsize = int(loopdir)
        print('gapLen =', loopdir)
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
        clu = Cluster(**kwargs)
        try:
            clu.run_protocol(pdbs)
        except:
            # make empty clusters directory if clustering fails
            if not os.path.exists(path + loopdir + '/clusters/'):
                os.mkdir(path + loopdir + '/clusters/')
            else:
                log.info('loop run_cluster error.')
                log.info('The path: (' + path + loopdir + '/clusters/' +  ') already exist.')
    if outfile:
        sys.stdout = orig_out
