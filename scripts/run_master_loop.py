## activate the conda environment. conda activate env_smallprot
## please check/change the parameters and file paths according to your purpose.
## you can use ipython or just 'python run_loop.py'
import sys
sys.path.append(r'/mnt/e/GitHub_Design/Qbits')
sys.path.append(r'/mnt/e/GitHub_Design/smallprot')
import os
import prody as pr
import matplotlib.pyplot as plt
import numpy as np

from smallprot.basic import cluster_loops
from smallprot.loop import loop_sse
from smallprot.external import query, extract_master
from smallprot.basic import plot

print('Thanks for using smallprot!')

target_list='/mnt/e/GitHub_Design/master_db/list'
master_query_loop_top = 200
rmsdCut = 1.0

workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta/output_sel/_rosetta_3rdRound/output_55F_sel/53_69/'
seed_pdb = workdir + '53_69.pdb'
seed_pdb_len = len(pr.parsePDB(seed_pdb).select('protein and name CA'))

outdir = workdir + 'output/'
os.makedirs(outdir, exist_ok=True)

loop_outfile = None
loop_range = [3, 10]

clusters_exist = loop_sse.loop_search_query_search(outdir, seed_pdb, loop_outfile, loop_range, target_list, rmsdCut, master_query_loop_top, seed_pdb_len = seed_pdb_len)

def write_top_aa(workdir, df, top = 5, filename = '_top_seq.txt'):
    '''
    For the logo frequency, output the top aa.
    '''
    with open(workdir + filename, 'w') as f:
        for i in range(df.shape[0]):
            idx = np.argsort(df.iloc[i])
            aas = np.flip(df.columns[idx][-top:])
            f.write(str(i) + '\t' + ''.join([a for a in aas]) + '\n')
    return

def main(workdir, seqfile = 'seq.txt', out_name = '_info.png'):
    all_seqs, all_rmsds = extract_master._get_seq_rmsd(workdir + seqfile)
    all_seqs_filter = [s for s in all_seqs if len(s) == 16]
    fig, (ax1) =plt.subplots(1, 1, figsize=(15, 6))
    df = plot._plot_log(fig, ax1, all_seqs_filter)
    plt.tight_layout()
    plt.savefig(workdir + out_name )
    plt.close()
    write_top_aa(workdir, df, top = 5, filename = '_top_seq.txt')

