import os
import numpy as np
import prody as pr
from ..basic import cluster, transformation
import matplotlib.pyplot as plt
from dataclasses import dataclass


def superimpose_aa_core(pdbs, rmsd, align_sel, min_cluster_size = 0, len_sel = 5):
    '''
    '''
    clu = cluster.Cluster()
    clu.rmsd_cutoff = rmsd
    clu.pdbs = []

    for pdb in pdbs:
        if align_sel is not None:
            c = pdb.select(align_sel).getCoords() 
        else:
            c=pdb.getCoords()       
        # if len(c)!= len_sel: 
        #     continue
        clu.pdb_coords.append(c)
        clu.pdbs.append(pdb)
    if len(clu.pdb_coords) <= 0:
        return 
    clu.pdb_coords = np.array(clu.pdb_coords, dtype = 'float32')

    clu.make_pairwise_rmsd_mat()  
    if not clu._square:
        clu.make_square()
    if not clu._adj_mat:
        #clu.make_adj_mat()
        clu.make_adj_mat_no_superpose()
    clu.cluster(min_cluster_size = min_cluster_size)

    return clu

def print_cluster_pdbs(clu, outdir, superimpose = True, tag = ''):

    for i in range(len(clu.mems)):
        cluster_out_dir = outdir + str(i) + '/'
        # if not os.path.exists(cluster_out_dir):
        #     os.mkdir(cluster_out_dir)
        if superimpose:     
            _print_cluster_rank_pdbs(clu, i, cluster_out_dir, tag)
        else:
            _print_cluster_rank_pdbs_nosupper(clu, i, cluster_out_dir, tag)
    return

def _print_cluster_rank_pdbs(clu, rank, outdir='./', tag=''):
    try: os.makedirs(outdir)
    except: pass
    try:
        cent = clu.cents[rank]
        mems = clu.mems[rank]

        # Align backbone of cluster centroid to backbone of centroid of largest cluster.
        R, m_com, t_com = transformation.get_rot_trans(clu.pdb_coords[cent],
                                        clu.pdb_coords[clu.cents[0]])
        cent_coords = np.dot((clu.pdb_coords[cent] - m_com), R) + t_com

        for i, mem in enumerate(mems):
            R, m_com, t_com = transformation.get_rot_trans(clu.pdb_coords[mem], cent_coords)
            pdb = clu.pdbs[mem].copy()
            pdb_coords = pdb.getCoords()
            coords_transformed = np.dot((pdb_coords - m_com), R) + t_com
            pdb.setCoords(coords_transformed)
            is_cent = '_centroid' if mem == cent else ''
            pr.writePDB(outdir + tag + 'cluster_' + str(rank) + '_mem_' + str(i)
                         + is_cent + '_' + pdb.getTitle().split('.')[0] + '.pdb', pdb)

    except IndexError:
        print('Cluster', rank, 'does not exist.')

def _print_cluster_rank_pdbs_nosupper(clu, rank, outdir='./', tag=''):
    try: os.makedirs(outdir)
    except: pass
    try:
        cent = clu.cents[rank]
        mems = clu.mems[rank]
        for i, mem in enumerate(mems):
            pdb = clu.pdbs[mem]
            is_cent = '_centroid' if mem == cent else ''
            pr.writePDB(outdir + tag + 'cluster_' + str(rank) + '_mem_' + str(i)
                         + is_cent + '_' + pdb.getTitle().split('.')[0] + '.pdb', pdb)
    except IndexError:
        print('Cluster', rank, 'does not exist.')
    return

@dataclass
class clu_info:
    clu_type: str
    clu_rmsd: float
    total_num: int
    clu_rank: int
    clu_num: int
    score: float

    def to_tab_string(self):
        clu_info = self.clu_type + '\t' + str(self.clu_rmsd) + '\t' + str(self.total_num) + '\t' + str(self.clu_rank) + '\t' + str(self.clu_num) + '\t'+ str(self.score) 
        return clu_info
    
def get_clu_info_write(outfile, pdbs, clu, rmsd, align_sel = 'heavy'):
    clu_infos = []
    n_avg = sum([len(m) for m in clu.mems])/len(clu.mems)
    for i in range(len(clu.mems)):
        c = clu_info(clu_type = align_sel, 
            clu_rmsd=rmsd, total_num = len(pdbs), clu_rank = i, 
            clu_num=len(clu.mems[i]), score = np.log(len(clu.mems[i])/n_avg) )
        clu_infos.append(c)
    write_clu_info(outfile, clu_infos)
    return clu_infos

def write_clu_info(filename, clu_infos):
    '''
    Write information of cluters.
    @ clu_infos: [clu_info]
    '''
    with open(filename, 'w') as f:
        f.write('clu_type\tclu_rmsd\ttotal_num\tclust_rank\tclust_num\tscore\n')
        for r in clu_infos:
            f.write(r.to_tab_string() + '\n')  
    return

def plot_clu_info(clu_infos, outplot):
    fig, (ax, ax1) =plt.subplots(2, 1, figsize=(12,8))
    
    x = list(range(1, len(clu_infos) + 1))
    y = [c.score for c in clu_infos]
    ax.plot(x, y)
    ax.hlines(y=0, xmin = 0, xmax = x[-1], color='r')
    ax.legend()
    #ax.set_xticks(x)
    ax.set_xlabel("Rank", fontsize = 12)
    ax.set_ylabel("vdM score", fontsize = 12)

    counts = [c.clu_num for c in clu_infos]
    ax1.bar(x, counts)

    plt.tight_layout()
    plt.savefig(outplot)
    plt.close()
    return