'''
The script here is used to remove duplicate transmembrane proteins in OPM(https://opm.phar.umich.edu/).

'''

import os
import prody as pr
import pandas as pd
import shutil
import numpy as np
from metalprot.basic import utils

from pymol import cmd
import glob

### Copy selected pdb (resolution <= 3.5) files into a new directory. 
def sel_pdb_by_resolution():
    workdir = '/mnt/e/DesignData/tm/database/'
    proteins = set()

    for prot in os.listdir(workdir + 'tms/'):
        if '.pdb' not in prot:
            continue
        proteins.add(prot)


    protein_info = pd.read_csv(workdir + 'proteins-2022-01-31_filter3.5.csv')

    outdir = workdir + 'tms_filter35/'
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
not_in_db_cr = ['1e12', '7e14']
for prot in not_in_db_cr:
    if prot + '.pdb' in proteins:
        shutil.copy(workdir + 'tms/' + prot + '.pdb', outdir + prot + '.pdb')
    else:
        not_in_db.append(prot)
        print('Protein pdb: ' + prot + ' not in database')
'''

def ext_protein_by_pymol(workdir, outdir):
    '''
    There are a few proteins that prody could not read. 
    '''
    cmd.cd(outdir)
    all_files = []

    for _file in glob.glob(os.path.join(workdir, '*.pdb')):
        all_files.append(_file)

    for _file in all_files:
        cmd.load(_file)
        objs = cmd.get_names()
        
        obj = objs[0]   ### only the first entry contains protein.

        obj_prot = obj + '_prot'
        count = cmd.count_atoms('model %s and polymer.protein' % (obj))
        if count <= 10:
            print(_file)
            continue

        cmd.create(obj_prot ,'model %s and polymer.protein' % (obj))

        cmd.save(_file.split('.')[0] + '.pdb', obj_prot)
        cmd.delete("all")
    return

#workdir = '/mnt/e/DesignData/tm/database/tms_filter35_prots/'
#outdir = '/mnt/e/DesignData/tm/database/tms_filter35_prots/'
#ext_protein_by_pymol(workdir, outdir)

### Extracting transmembrane domain

def merge_prody_objs(pdbs, title):
    '''
    combine vdms into target and mutate all other aa 2 ala or gly for ligand position.
    '''

    ag = pr.AtomGroup(title)
    coords = []
    chids = []
    names = []
    resnames = []
    resnums = []

    for c in pdbs:
        coords.extend(c.getCoords())
        chids.extend(c.getChids())
        names.extend(c.getNames())
        resnames.extend(c.getResnames())
        resnums.extend(c.getResnums())

    ag.setCoords(np.array(coords))
    ag.setChids(chids)
    ag.setNames(names)
    ag.setResnames(resnames)
    ag.setResnums(resnums)
    return ag

def ext_tm_domain(pdb, outdir, zdist = 18):
    '''
    The function is used to extract transmembrane domain of pdb from OPM database.
    zdist is the z axis for the tm domain cutoff.
    '''
    print('Extract tm domain ' + pdb.getTitle())
    chid_resns = {}  
    for chid in np.unique(pdb.select('protein').getChids()):
        #ind_sel = np.unique(pdb.select('protein and chid ' + chid + ' and name CA and -' + str(zdist) + '<=z<=' + str(zdist)).getResindices())
        _sel = pdb.select('protein and chid ' + chid + ' and name CA and -' + str(zdist) + '<=z<=' + str(zdist))
        if _sel is None or len(_sel) <= 0:
            print('The chid is empty ' + chid)
            continue
        ind_sel = np.unique(_sel.getResnums())
        
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
    _pdbs = []
    for chid in chid_resns.keys():
        pdb_sel = pdb.select('protein and chid ' + chid + ' and resnum ' + ' '.join([str(x) for x in chid_resns[chid]])).copy().toAtomGroup()       
        if not pdb_sel.select('name CA and -5<=z<=5'):
            continue
        _pdbs.append(pdb_sel)
        
    if len(_pdbs) <= 0:
        return False
    _pdb = merge_prody_objs(_pdbs, title=pdb.getTitle())

    pr.writePDB(outdir + _pdb.getTitle(), _pdb)
    return True


def run():
    workdir = '/mnt/e/DesignData/tm/database/'
    outdir = workdir + 'tms_filter35_tmdomain/'
    os.makedirs(outdir, exist_ok=True)
    files = []
    for file in os.listdir(workdir + 'tms_filter35_prots/'):
        if '.pdb' not in file:
            continue
        files.append(file)

    pdbs_no_tmdomain = []
    failed = []

    for file in files:
        try:
            pdb = pr.parsePDB(workdir + 'tms_filter35_prots/' + file)
            x = ext_tm_domain(pdb, outdir, zdist=18)
            if not x:
                pdbs_no_tmdomain.append(file)
        except:
            failed.append(file)

    with open(outdir + '_failed.txt' , 'w') as f:
        f.write('\n'.join(failed))

    with open(outdir + '_notmdomain.txt' , 'w') as f:
        f.write('\n'.join(pdbs_no_tmdomain))
    return

def ext_tm_domain_by_pymol():
    '''
    There are a few proteins that prody could not read. 
    Not implemented.
    '''
    return

'''
python /mnt/e/GitHub_Design/Rvrsprot/scripts_tm/TM_database_generation.py
'''
if __name__=='__main__':
    print('run main function.')
    run()
    