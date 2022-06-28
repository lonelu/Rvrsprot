import os
import prody as pr
from rvrsprot.external import query as qr

'''
python /mnt/e/GitHub_Design/Rvrsprot/scripts/qbit_search/qbit_search_multiple.py
'''

# file path for clustered chains dictionary used for handling redundancy
targetList = '/mnt/e/DesignData/Database/Qbits/pds_list_2p5.txt'
path_to_chains_dict = '/mnt/e/DesignData/Database/Qbits/db_2p5A_0p3rfree_chains_dictionary.pkl'

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c3/c3_vc0_sel/'

def qbits(file):

    outdir = workdir + file.split('.')[0] + '/'
    os.makedirs(outdir, exist_ok=True)
    pdb = pr.parsePDB(workdir + file)

    pr.writePDB(outdir + 'query.pdb', pdb.select('protein and chid A B and resnum 5 6 7 8 9 10 11').toAtomGroup())
    pr.writePDB(outdir + 'query_full.pdb', pdb.select('protein and chid A B C and resnum 5 6 7 8 9 10 11').toAtomGroup())

    query = outdir + 'query.pdb'
    query_full = outdir + 'query_full.pdb'

    outdir_master = outdir + 'output_master/'
    match = outdir_master + 'match.txt'
    seq = outdir_master + 'seq.txt'
    os.makedirs(outdir_master, exist_ok=True)
    qr.master_query(outdir_master, query, targetList, rmsdCut=0.5, topN=None, outfile= outdir_master + 'stdout.txt')

    with open(seq, 'r') as f:
        if len(f.readlines()) <= 0:
            return
    outdir_qbit = outdir + 'output_qbit/'
    os.makedirs(outdir_qbit, exist_ok=True)
    
    qr.qbits_search(query, query_full, outdir_master, 
                    path_to_chains_dict, outdir_qbit, 
                    window_size = 10, rmsd = 0.75, 
                    top = 200, sec_struct="H",
                    antiparallel=True, min_nbrs = 1, contiguous = True)

    return

for file in os.listdir(workdir):
    if not '.pdb' in file:
        continue
    qbits(file)

#file = 'A_22_C3.pdb'
#qbits(file)

# query_pdb_path = query
# query_full_pdb_path = query_full
# query_dir = outdir_master
# chains_dict_path = path_to_chains_dict
# outdir = outdir_qbit
# window_size=10
# rmsd=1.5
# top=5
# sec_struct=None
# antiparallel=False
# min_nbrs=1
# contiguous=False