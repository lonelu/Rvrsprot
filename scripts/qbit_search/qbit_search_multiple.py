import os
import prody as pr
from rvrsprot.external import query as qr


# file path for clustered chains dictionary used for handling redundancy
targetList = '/mnt/e/DesignData/Database/Qbits/pds_list_2p5.txt'
path_to_chains_dict = '/mnt/e/DesignData/Database/Qbits/db_2p5A_0p3rfree_chains_dictionary.pkl'

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c3_vc0_sel/'

for file in os.listdir(workdir):
    if not '.pdb' in file:
        continue
    outdir = workdir + file.split('.')[0] + '/'
    os.makedirs(outdir, exist_ok=True)
    pdb = pr.parsePDB(workdir + file)

    pr.writePDB(outdir + pdb.getTitle() + '_AB.pdb', pdb.select('protein and chid A B and resnum 5 6 7 8 9 10 11'))
    pr.writePDB(outdir + pdb.getTitle() + '_all.pdb', pdb)

    query = outdir + pdb.getTitle() + '_AB.pdb'
    query_full = outdir + pdb.getTitle() + '_all.pdb'


    outdir_master = outdir + 'output_master/'
    match = outdir_master + 'match.txt'
    seq = outdir_master + 'seq.txt'
    os.makedirs(outdir_master, exist_ok=True)
    qr.master_query(outdir_master, query, targetList, rmsdCut=0.5, topN=None, outfile= outdir_master + 'stdout.txt')

    with open(seq, 'r') as f:
        if len(f.readlines()) <= 0:
            continue
    outdir_qbit = outdir + 'output_qbit/'
    os.makedirs(outdir_qbit, exist_ok=True)
    try:
        qr.qbits_search(query, query_full, outdir_master, 
                        path_to_chains_dict, outdir_qbit, 
                        window_size = 10, rmsd = 0.75, 
                        top = 1000, sec_struct="H",
                        antiparallel=False, min_nbrs = 1, contiguous = True)
    except:
        print('____________________Qbit Wrong________________________')