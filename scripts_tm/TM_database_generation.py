'''
The script here is used to remove duplicate transmembrane proteins in OPM(https://opm.phar.umich.edu/).

'''

import os
import prody as pr
import pandas as pd
import shutil
import numpy as np
from metalprot.basic import utils



### Copy selected pdb (resolution <= 3.5) files into a new directory. 
def sel_pdb_by_resolution():
    workdir = '/mnt/e/DesignData/tm/database/'
    proteins = set()

    for prot in os.listdir(workdir + 'tms/'):
        if '.pdb' not in prot:
            continue
        proteins.add(prot)


    protein_info = pd.read_csv(workdir + 'proteins-2022-01-31_filter3.5.csv')

    outdir = workdir + 'tms_filter3.5/'
    os.makedirs(outdir, exist_ok=True)

    not_in_db = []
    for prot in protein_info['pdbid']:
        if prot + '.pdb' in proteins:
            shutil.copy(workdir + 'tms/' + prot + '.pdb', outdir + prot + '.pdb')
        else:
            not_in_db.append(prot)
            print('Protein pdb: ' + prot + ' not in database')

'''
### correctly copy the failed entries.
not_in_db_cr = ['1e12', '2e74', '3e86', '5e94', '6e59', '7e33', '7e32', '7e26', '7e27', '7e14']
for prot in not_in_db_cr:
    if prot + '.pdb' in proteins:
        shutil.copy(workdir + 'tms/' + prot + '.pdb', outdir + prot + '.pdb')
    else:
        not_in_db.append(prot)
        print('Protein pdb: ' + prot + ' not in database')
'''

### Extracting transmembrane domain

def ext_tm_domain(pdb, outdir, zdist = 19):
    '''
    The function is used to extract transmembrane domain of pdb from OPM database.
    zdist is the z axis for the tm domain cutoff.
    '''
    chid_resns = {}  
    for chid in np.unique(pdb.select('protein').getChids()):
        ind_sel = np.unique(pdb.select('protein and chid ' + chid + ' and name CA and -' + str(zdist) + '<=z<=' + str(zdist)).getResindices())
        ind_sel = np.unique(pdb.select('protein and name CA and -' + str(zdist) + '<=z<=' + str(zdist)).getResnums())
        
        #TO DO: rm  len < 7
        ind_list = []
        inds = []
        for i in range(len(ind_sel)):
            if len(inds) ==0:
                inds.append(ind_sel[i])
            else:
                if ind_sel[i] == inds[-1] + 1:
                    inds.append(ind_sel[i])
                    if i == len(ind_sel)-1:
                        ind_list.append([x for x in inds])
                else:
                    if len(inds) >= 7:
                        ind_list.append([x for x in inds])
                    inds = []

        #TO DO: rm no helix. The count of connected A is more than 7 aa.
          
        for inds in ind_list:
            p_sel = pdb.select('protein and chid ' + chid + ' and resnum ' + ' '.join([str(x) for x in inds])).toAtomGroup()
            abples, phipsi = utils.seq_get_ABPLE(p_sel)
            if  'AAAAAAA' in ''.join(abples):
                if chid in chid_resns.keys():
                    chid_resns[chid].extend(inds)
                else:
                    chid_resns[chid] = inds

    if len(chid_resns) <= 0:
        return False

    for chid in chid_resns.keys():
        pdb_sel = pdb.select('protein and chid ' + chid + ' and resnum ' + ' '.join([str(x) for x in chid_resns[chid]]))       
        if not pdb_sel.select('name CA and -5<=z<=5'):
            return False
        pr.writePDB(outdir + pdb.getTitle() + '_' + chid, pdb_sel)
    return True


def run():
    workdir = '/mnt/e/DesignData/tm/database/'
    outdir = workdir + 'tms_filter3.5_tmdomain2/'
    os.makedirs(outdir, exist_ok=True)

    pdbs_no_tmdomain = []
    failed = []
    for file in os.listdir(workdir + 'tms_filter3.5/'):
        if '.pdb' not in file:
            continue
        try:
            pdb = pr.parsePDB(workdir + 'tms_filter3.5/' + file)
            x = ext_tm_domain(pdb, outdir, zdist=19)
            if not x:
                pdbs_no_tmdomain.append(file)
        except:
            failed.append(file)

    with open(outdir + '_failed.txt' , 'w') as f:
        f.write('\n'.join(failed))

    with open(outdir + '_notmdomain.txt' , 'w') as f:
        f.write('\n'.join(pdbs_no_tmdomain))
    

def prepare_fasta():
    '''
    prepare fasta file of each chain each tm protein.
    '''
    workdir = '/mnt/e/DesignData/tm/database/'
    pdbs = []
    failed = []
    for file in os.listdir(workdir + 'tms_filter3.5/'):
        if '.pdb' not in file:
            continue
        try:
            pdb = pr.parsePDB(workdir + 'tms_filter3.5/' + file)
            pdbs.append(pdb)

        except:
            failed.append(file)

    with open(workdir + 'allseq.fasta', 'w') as f:
        for pdb in pdbs:
            for chid in pdb.select('protein').getChids():
                c = pdb.select('protein and chid ' + chid)
                f.write('>' + pdb.getTitle() + '_' + chid + '\n')
                f.write(c.select('name CA').getSequence() + '\n')

from pymol import cmd
import glob
import os


def prepare_fasta_by_pymol():

    workdir = '/mnt/e/DesignData/tm/database/tms_filter3.5/'  

    outdir = workdir + 'seq/'
    os.makedirs(outdir, exist_ok=True)
    cmd.cd(outdir)

    for _file in glob.glob(os.path.join(workdir, '*.pdb')):
        cmd.load(_file)
        obj = cmd.get_names()[0]
        print(obj)

        for ch in cmd.get_chains(obj):
            #if len(ch) >=1:
            if len(ch) < 1:
                name = obj + '_' + ch
                cmd.create(name, 'model %s and chain %s and polymer.protein' % (obj, ch))
                cmd.save(name + '.fasta', name)
                cmd.delete(name)
        cmd.delete(obj)
        

