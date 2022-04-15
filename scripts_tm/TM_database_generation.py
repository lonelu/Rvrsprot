'''
The script here is used to remove duplicate transmembrane proteins in OPM(https://opm.phar.umich.edu/).

'''

import os
import prody as pr
import pandas as pd
import shutil
import numpy as np
from metalprot.basic import utils

workdir = '/mnt/e/DesignData/tm/database/'

### Copy selected pdb (resolution <= 3.5) files into a new directory. 
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
    ind_sel = np.unique(pdb.select('protein and name CA and -' + str(zdist) + '<=z<=' + str(zdist)).getResindices())
    
    #TO DO: rm  len < 7
    ind_list = []
    inds = []
    for i in range(len(ind_sel)-1):
        if len(inds) ==0:
            inds.append(ind_sel[i])
        else:
            if ind_sel[i] == inds[-1] + 1:
                inds.append(ind_sel[i])
            else:
                if len(inds) >= 7:
                    ind_list.append([x for x in inds])
                inds = []

    #TO DO: rm no helix. The count of connected A is more than 7 aa.
    ind_sel2 = []    
    for inds in ind_list:
        p_sel = pdb.select('resindex ' + ' '.join([str(x) for x in inds])).toAtomGroup()
        abples, phipsi = utils.seq_get_ABPLE(p_sel)
        if  'AAAAAAA' in ''.join(abples):
            ind_sel2.extend(inds)

    if len(ind_sel2) <= 0:
        return False

    pdb_sel = pdb.select('resindex ' + ' '.join([str(x) for x in ind_sel2]))
    
    if not pdb_sel.select('name CA and -5<=z<=5'):
        return False
    #pr.writePDB(outdir + pdb.getTitle(), pdb_sel)
    return True

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
    
