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

workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/search_unknow_result/output__helix6a_10-14-139_HDH.pdb_20220624-222024/loop_ss_20220625-103345/merge/'

target = pr.parsePDB(workdir + 'merge_before.pdb')

loopfiles = ['A-16-22-A-48-54_cent_rg_14_clu_19.pdb', 'match4.pdb', 'A-91-97-A-118-124_cent_rg_6_clu_18.pdb']

loops = [pr.parsePDB(workdir + f) for f in loopfiles]

for loop in loops:
    chidreses, values = calc_ca_dist(target, loop)
    print('--------------------------------------')
    print(chidreses)
    print(values)

(19, 100) (49, 119)
(68, 48)  (81, 63)
(91, 74) (121, 90)


ag1 = target.select('protein').select('resnum ' + ' '.join(str(x) for x in range(0, 19))).toAtomGroup()

ag2 = loops[0].select('resnum ' + ' '.join(str(x) for x in range(100, 120))).toAtomGroup()

ag3 = target.select('resnum ' + ' '.join(str(x) for x in range(50, 68))).toAtomGroup()

ag4 = loops[1].select('resnum ' + ' '.join(str(x) for x in range(48, 64))).toAtomGroup()

ag5 = target.select('resnum ' + ' '.join(str(x) for x in range(82, 91))).toAtomGroup()

ag6 = loops[2].select('resnum ' + ' '.join(str(x) for x in range(74, 91))).toAtomGroup()

ag7 = target.select('resnum ' + ' '.join(str(x) for x in range(122, 149))).toAtomGroup()

ag8 = target.select('not protein').toAtomGroup()
ags = [ag1, ag2, ag3, ag4, ag5, ag6, ag7, ag8]
ag_all = prody_ext.combine_ags_into_one_chain(ags, 'merge')

[pr.writePDB(workdir + 'ag_' + str(i) + '.pdb', ags[i])  for i in range(len(ags))]

pr.writePDB(workdir + 'merge.pdb', ag_all)

ag_all_ala = prody_ext.target_to_all_gly_ala(ag_all, 'merge_allala', keep_no_protein=True)
