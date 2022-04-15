import os
from rvrsprot.external import query as qr


# file path for clustered chains dictionary used for handling redundancy
targetList = '/mnt/e/DesignData/Database/Qbits/pds_list_2p5.txt'
path_to_chains_dict = '/mnt/e/DesignData/Database/Qbits/db_2p5A_0p3rfree_chains_dictionary.pkl'

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/helix_noclash_0_0_0_search/A_36_B_47_C_81/qbits/'

query = workdir + 'query.pdb'
query_full = workdir + 'query_full.pdb'


outdir_master = workdir + 'output/'
match = outdir_master + 'match.txt'
seq = outdir_master + 'seq.txt'
#os.makedirs(outdir_master, exist_ok=True)
#qr.master_query(outdir_master, query, targetList, rmsdCut=0.5, topN=None, outfile= outdir_master + 'stdout.txt')

outdir_qbit = workdir + 'output_qbit/'
os.makedirs(outdir_qbit, exist_ok=True)

qr.qbits_search(query, query_full, outdir_master, 
                path_to_chains_dict, outdir_qbit, 
                window_size = 10, rmsd = 0.75, 
                top = 1000, sec_struct="H",
                antiparallel=False, min_nbrs = 1, contiguous = True)