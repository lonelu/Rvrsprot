'''
From 200 aa long 6 helix bundles extract portions for tm protein design.
'''

import os
import prody as pr
import numpy as np
from scipy.spatial.distance import cdist

workdir = '/mnt/e/DesignData/tm/Kehan/'

pdb = pr.parsePDB(workdir + 'hp1_85.pdb')

chain_len = len(pdb.select('chid A and name CA'))

skip = 3
conn = 18
chids = ['A', 'B', 'C', 'D', 'E', 'F']

sel_chidresnum_set = []
for i in range(0, chain_len - conn, skip):
    sel_chidresnums = ''
    resnums_odd = [i+ c for c in range(conn)]
    resnums_evn = [chain_len-i-c for c in range(conn)]

    for chid in chids:
        if chid in ['A', 'C', 'E']:
            sel_chidresnums += 'chid ' + chid + ' and resnum ' + ' '.join([str(x) for x in resnums_odd])
        else:
            sel_chidresnums += 'chid ' + chid + ' and resnum ' + ' '.join([str(x) for x in resnums_evn])
        sel_chidresnums += ' or '
    sel_chidresnum_set.append(sel_chidresnums[:-4])

outdir = workdir + 'all/'
os.makedirs(outdir, exist_ok=True)
for i in range(len(sel_chidresnum_set)):
    pr.writePDB(outdir + 'test_' + str(i)  + '.pdb', pdb.select(sel_chidresnum_set[i]))

pdb = pr.parsePDB(outdir + 'test_0pdb.pdb')

def calc_mindist_between_chid(pdb, chid1 = 'A', chid2 = 'B'):
    '''
    As the name suggested.
    '''
    chid1_coord = pdb.select('name CA and chid ' + chid1).getCoords()
    chid2_coord = pdb.select('name CA and chid ' + chid2).getCoords()

    y = cdist(chid1_coord, chid2_coord)
    return y.min()


def filter_by_chain_dist(pdb):
    '''
    Special usage for 6 helix bundles 'A B C D E F'
    '''

    dists_A = [calc_mindist_between_chid(pdb, 'A', x) for x in ['B', 'C', 'D', 'E', 'F']]

    dists_B = [calc_mindist_between_chid(pdb, 'B', x) for x in ['A', 'C', 'D', 'E', 'F']]

    if sum(np.array(dists_A) < 6) > 2 or sum(np.array(dists_B) < 6) > 2:
        return False
    return True


outdir_filter = workdir + 'filter/'
os.makedirs(outdir_filter, exist_ok= True)

for file in os.listdir(outdir):
    if '.pdb' not in file:
        continue
    pdb = pr.parsePDB(outdir + file)

    if filter_by_chain_dist(pdb):
        pr.writePDB(outdir_filter + pdb.getTitle() + '.pdb', pdb)
    

### Extract 3 Helixs from the 6 helix bundles. 
outdir_3helix = workdir + 'filter_3helix/'
os.makedirs(outdir_3helix, exist_ok= True)

for file in os.listdir(outdir_filter):
    if '.pdb' not in file:
        continue
    pdb = pr.parsePDB(outdir_filter + file)
    dists_B = [calc_mindist_between_chid(pdb, 'B', x) for x in ['C', 'D', 'E', 'F']]

    min_dist2B_chid = 'B'
    min_dist2B = np.Infinity
    for i in range(len(dists_B)):
        if  dists_B[i] < min_dist2B:
            min_dist2B = dists_B[i]
            min_dist2B_chid = ['C', 'D', 'E', 'F'][i]

    pdb_sel = pdb.select('chid A B ' + min_dist2B_chid)
    pr.writePDB(outdir_3helix + pdb.getTitle() + '.pdb', pdb_sel)

