'''
Extract his-h2o-asp and superimpose on COMBS ideal ala.
'''
import os
import prody as pr
import itertools
from rvrsprot.basic import constant

import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool

'''
python /Users/lonelu/GitHub_Design/Rvrsprot/scripts/construct_h2o_vdm/h2o_vdm_extract.py
'''

def get_pairs(resnums, chids, distance=4):
    pairs = []
    
    for cb in itertools.combinations(list(range(len(resnums))), 2):
        chid_a = chids[cb[0]]
        chid_b = chids[cb[1]]
        rn_a = resnums[cb[0]]
        rn_b = resnums[cb[1]]

        if chid_a == chid_b and abs(rn_a -rn_b) < distance:
            continue
        pairs.append((chid_a, str(rn_a), chid_b, str(rn_b)))
    return pairs


def _extract_h2o_vdm(pdb):
    failed_extraction = []
    h2o_chids = pdb.select('water').getChids()
    h2o_resnums = pdb.select('water').getResnums()

    for i_h2o in range(len(h2o_resnums)):
        h2o_rn = h2o_resnums[i_h2o]
        h2o_chid = h2o_chids[i_h2o]
        all_near = pdb.select('protein and within 3.6 of ( chid ' + h2o_chid + ' and resnum ' + str(h2o_rn) + ' )')
        if not all_near: 
            continue

        all_near = all_near.select('nitrogen or oxygen')
        if not all_near: 
            continue
        chids = all_near.getChids()
        inds = all_near.getResnums()
        if len(inds) < 2:
            continue

        for pair in get_pairs(inds, chids):
            sel_str = '( chid ' + pair[0] + ' resnum ' + pair[1] + ' ) or ( chid ' + pair[2] + ' resnum ' + pair[3] + ' ) or ( chid ' + h2o_chid + ' and resnum ' + str(h2o_rn) + ' )'
            try:
                resnames = pdb.select(sel_str).getResnames()
            except:
                print(pdb_name + ' | ' + sel_str )
                failed_extraction.append(pdb_name + ' \t ' + sel_str + '\n')
                continue
            if not 'ASP' in resnames or not 'TRP' in resnames:
                continue
            ext = pdb.select(sel_str)

            coords = ext.select('resname ASP and name N CA C').getCoords()
            if len(coords) != 3:
                continue
            tf = pr.calcTransformation(coords, constant.ideal_ala_coords)
            tf.apply(ext)
            name = pdb_name.split('.')[0] + '_' + pair[0] + pair[1] + '_' + pair[2] + pair[3] + '.pdb'
            pr.writePDB(workdir_vdM + name, ext)
    return failed_extraction
    

def extract_h2o_vdm(pdb_name):
    failed_extraction = []

    pdb = pr.parsePDB(workdir_ext + pdb_name)
    print(pdb_name)
    if not pdb.select('water'):
        failed_extraction.append(pdb_name + '\t no water\n')
        return failed_extraction
    
    failed_extraction.extend(_extract_h2o_vdm(pdb))

    return failed_extraction


workdir = '/Users/lonelu/DesignData/Database/vdMh2o/'
workdir_ext = workdir + 'parent_pdb_exts/'
#workdir_vdM = workdir + 'vdm_pdb_TYR-pro/'
workdir_vdM = workdir + 'vdm_pdb_ASP_indole/'
os.makedirs(workdir_vdM, exist_ok=True)


extracted = set()
for file in os.listdir(workdir_vdM):
    if file.endswith(".pdb"):
        extracted.add(file.split('_')[0].upper())


pdbs = []
for pdb_name in os.listdir(workdir_ext):
    if '.pdb' not in pdb_name:
        continue

    #if pdb_name.split('.')[0].upper() in extracted:
    #    continue

    pdbs.append(pdb_name)

'''
### Somehow using pool run into an error.
num_cores = int(mp.cpu_count() - 1)
print('pool: {}'.format(num_cores))
pool = ThreadPool(num_cores)
results = pool.map(extract_h2o_vdm, pdbs)
pool.close()
pool.join()
'''

results = []
for pdb_name in pdbs:
    result = extract_h2o_vdm(pdb_name)
    results.append(result)


with open(workdir + 'failed_extraction_ASP-indole.txt', 'w') as f:
    for fes in results:
        for fe in fes:
            f.write(fe)




