'''
After searching the loops for target proteins, previously loops are manually inserted.
The script here is to help quickly inserted the loops into the target.
'''

import os
import prody as pr
import numpy as np
from scipy.spatial.distance import cdist
from metalprot.basic import prody_ext

def calc_ca_dist(target, loop):
    target_coords = target.select('name CA').getCoords()
    loop_coords = loop.select('name CA').getCoords()

    dists = cdist(target_coords, loop_coords)
    
    chidreses = []
    values = []
    for i in range(loop_coords.shape[0]):
        j = np.argmin(dists[:,i])
        chidres = (target.select('resindex ' + str(j)).getResnums()[0], loop.select('resindex ' + str(i)).getResnums()[0])
        v = dists[j,i]
        
        chidreses.append(chidres)
        values.append(round(v, 2))
    return chidreses, values

workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/search_unknow_result/output__helix6a_10-14-139_HDH.pdb_20220624-222024/loop_ss_20220630-095815/merge/'

target = pr.parsePDB(workdir + 'helix6a_preloop.pdb')

loopfiles = ['A-93-99-A-124-130_cent_rg_13_clu_81.pdb', 'A-138-144-A-6-12_cent_rg_20_clu_14.pdb', 'A-22-28-A-48-54_cent_rg_8_clu_158.pdb']

loops = [pr.parsePDB(workdir + f) for f in loopfiles]

for loop in loops:
    chidreses, values = calc_ca_dist(target, loop)
    print('--------------------------------------')
    print(chidreses)
    print(values)

(96, 157), (128, 178)
(141, 94),  (10, 122)
(27, 70), (51, 83)

pdb_pre = pr.parsePDB(workdir + 'A-142-148-A-76-81_cent_rg_15_clu_39.pdb')
pdb_ext = pr.parsePDB(workdir + 'A-66-71-A-4-10_cent_rg_17_clu_22.pdb')

ag_pre = pdb_pre.select('resnum 51to58')
ag_ext = pdb_ext.select('resnum 352to359').toAtomGroup()


ag1 = target.select('protein').select('resnum ' + ' '.join(str(x) for x in range(77, 96))).toAtomGroup()

ag2 = loops[0].select('resnum ' + ' '.join(str(x) for x in range(157, 178 + 1))).toAtomGroup()

ag3 = target.select('resnum ' + ' '.join(str(x) for x in range(128 + 1, 141))).toAtomGroup()

ag4 = loops[1].select('resnum ' + ' '.join(str(x) for x in range(94, 122 + 1))).toAtomGroup()

ag5 = target.select('resnum ' + ' '.join(str(x) for x in range(10 + 1, 27))).toAtomGroup()

ag6 = loops[2].select('resnum ' + ' '.join(str(x) for x in range(70, 83 + 1))).toAtomGroup()

ag7 = target.select('resnum ' + ' '.join(str(x) for x in range(51 + 1, 70))).toAtomGroup()

ag8 = target.select('not protein').toAtomGroup()
ags = [ag_pre, ag1, ag2, ag3, ag4, ag5, ag6, ag7, ag_ext, ag8]
ag_all = prody_ext.combine_ags_into_one_chain(ags, 'merge2')

[pr.writePDB(workdir + 'ag_' + str(i) + '.pdb', ags[i])  for i in range(len(ags))]

pr.writePDB(workdir + 'merge2.pdb', ag_all)

ag_all_ala = prody_ext.target_to_all_gly_ala(ag_all, 'merge_allala', keep_no_protein=True)


