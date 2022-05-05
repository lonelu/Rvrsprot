import os
import sys
import multiprocessing as mp
import shutil
import prody as pr
import numpy as np
import itertools
import statistics
from rvrsprot.external import query


def master_pdb(workdir, helix_noclash_outdir, outdir_mid_name, pdb, targetList):
    print(pdb)
    _outdir = workdir + outdir_mid_name + pdb.split('.')[0] + '_all/'
    os.makedirs(_outdir, exist_ok=True)
    query.master_query(_outdir, helix_noclash_outdir + pdb, targetList = targetList, rmsdCut=1.0, topN=None, outfile=_outdir + 'stdout.txt')
    
    exist = False
    with open(_outdir + 'seq.txt') as f:
        lines = f.readlines()
        if len(lines) > 0:
            exist = True
    if not exist:
        shutil.rmtree(_outdir)
    return


def master_pdb_bypair(workdir, helix_noclash_outdir, outdir_mid_name, _pdb, targetList):
    '''
    If the input is three helix bundles, mater search each pair.
    '''
    pdb = pr.parsePDB(helix_noclash_outdir + _pdb)
    chids = np.unique(pdb.select('protein').getChids())
    for i, j in itertools.combinations(range(len(chids)), 2):
        a = chids[i]
        b = chids[j]
        
        _outdir = workdir + outdir_mid_name + pdb.getTitle() + '/' + a + '_' + b + '/'
        os.makedirs(_outdir, exist_ok= True)
        pdb_sel = pdb.select('chain ' + a + ' ' + b).copy()
        pr.writePDB(_outdir + a + '_' + b + '.pdb', pdb_sel)

        query.master_query(_outdir, _outdir + a + '_' + b + '.pdb', targetList = targetList, rmsdCut=0.75, topN=None, outfile=_outdir + 'stdout.txt')

    return

def find_best_all():
    '''
    After master search, we get results for all start structures (from master_pdb_bypair()).
    Then we need to summary the information to find the one with the most master hit.
    '''
    workdir = '/wynton/home/degradolab/lonelu/DesignData/_reverse_design/'
    _entrys = []
    for _d1 in os.listdir(workdir):
        if (not os.path.isdir(_d1)) or ('master_all' not in _d1):
            continue   
        print(_d1)

        for _d2 in os.listdir(workdir + _d1):
            if (not os.path.isdir(workdir + _d1 + '/' + _d2)) or ('cluster' not in _d2):
                continue   
            print(_d2) 

            rmsds = []
            with open(workdir + _d1 + '/' + _d2  + '/seq.txt', 'r') as f:         
                lines = f.readlines()
                if len(lines) <= 0:
                    continue
                for line in lines:
                    rmsds.append(float(line.strip(' ').split(' ')[0]))
            entry=(_d1, _d2, statistics.mean(rmsds), min(rmsds), len(lines))
            _entrys.append(entry)

    with open(workdir + '_master_all_summary.tsv', 'w') as f:
        f.write('topology\tentry\tChid-chid\tmean_rmsd\tmin_rmsd\tcount\n')
        for entry in _entrys:
            f.write('\t'.join(str(e) for e in entry) + '\n') 

    return 

def find_best():
    '''
    After master search, we get results for all start structures (from master_pdb_bypair()).
    Then we need to summary the information to find the one with the most master hit.
    '''
    workdir = '/wynton/home/degradolab/lonelu/DesignData/_reverse_design/'
    _entrys = []
    for _d1 in os.listdir(workdir):
        if (not os.path.isdir(_d1)) or ('master' not in _d1):
            continue   
        print(_d1)
        for _d2 in os.listdir(workdir + _d1):
            if (not os.path.isdir(workdir + _d1 + '/' + _d2)) or ('cluster' not in _d2):
                continue   
            print(_d2) 
            entry = []
            for _d3 in os.listdir(workdir + _d1 + '/' + _d2):
                rmsds = []
                with open(workdir + _d1 + '/' + _d2 + '/' + _d3 + '/seq.txt', 'r') as f:
                    lines = f.readlines()
                    if len(lines) <= 0:
                        continue
                    for line in lines:
                        rmsds.append(float(line.strip(' ').split(' ')[0]))
                entry.append((_d1, _d2, _d3, statistics.mean(rmsds), min(rmsds), len(lines)))
            _entrys.append(entry)

    with open(workdir + '_master_summary.tsv', 'w') as f:
        f.write('topology\tentry\tChid-chid\tmean_rmsd\tmin_rmsd\tcount\tChid-chid\tmean_rmsd\tmin_rmsd\tcount\tChid-chid\tmean_rmsd\tmin_rmsd\tcount\tTotal\n')
        for entry in _entrys:
            f.write(entry[0][0] + '\t' + entry[0][1] + '\t')
            f.write('\t'.join(str(e) for e in entry[0][2:]) + '\t') 
            f.write('\t'.join(str(e) for e in entry[1][2:]) + '\t') 
            f.write('\t'.join(str(e) for e in entry[2][2:]) + '\t' + str(entry[0][-1] + entry[1][-1] + entry[2][-1]) + '\n')
    return 

    
def main():
    workdir = '/home/gpu/Lei/DesignData/_reverse_design/'
    helix_noclash_outdir = workdir + 'helix_noclash_0-1_vdm/'
    outdir_mid_name = 'helix_noclash_0-1_vdm_master/'

    #pdbs = ['A_0_B_20_C_33.pdb' , 'A_109_B_13_C_51.pdb']
    pdbs = []
    for pdb in os.listdir(helix_noclash_outdir):
        if '.pdb' not in pdb:
            continue
        pdbs.append(pdb)

    pool = mp.Pool(mp.cpu_count())
    [pool.apply_async(master_pdb, args=(workdir, helix_noclash_outdir, outdir_mid_name, pdb, '/home/gpu/Nick/master_db/list')) for pdb in pdbs]
    pool.close()
    pool.join()
    return 


def main_wynton():
    workdir = '/wynton/home/degradolab/lonelu/DesignData/_reverse_design/'
    helix_noclash_outdir = workdir + 'helix_noclash_1-1-1_cent/'
    outdir_mid_name = 'helix_noclash_1-1-1_cent_master/'

    pdbs = []
    for pdb in os.listdir(helix_noclash_outdir):
        if '.pdb' not in pdb:
            continue
        pdbs.append(pdb)

    ind = int(sys.argv[1]) -1
    _pdb = pdbs[ind]

    master_pdb_bypair(workdir, helix_noclash_outdir, outdir_mid_name, _pdb, targetList='/wynton/home/degradolab/lonelu/DesignData/Database/master_db/list')

    return 


if __name__=='__main__':
    main_wynton()
