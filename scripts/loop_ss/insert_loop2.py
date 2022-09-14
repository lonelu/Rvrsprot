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

workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/search_unknow_result/output__helix6a_10-14-139_HDH.pdb_20220624-222024/loop_ss_20220630-211324/merge/'

target = pr.parsePDB(workdir + 'helix6a_preloop.pdb')

loopfiles = ['A-66-71-A-4-10_cent_rg_9_clu_105.pdb', 'A-16-22-A-128-134_cent_rg_17_clu_11.pdb', 'A-140-146-A-78-83_cent_rg_8_clu_188.pdb']

loops = [pr.parsePDB(workdir + f) for f in loopfiles]

for loop in loops:
    chidreses, values = calc_ca_dist(target, loop)
    print('--------------------------------------')
    print(chidreses)
    print(values)



ag1 = target.select('protein').select('resnum ' + ' '.join(str(x) for x in range(47, 71))).toAtomGroup()

ag2 = loops[0].select('resnum ' + ' '.join(str(x) for x in range(35, 47 + 1))).toAtomGroup()

ag3 = target.select('resnum ' + ' '.join(str(x) for x in range(6 + 1, 22))).toAtomGroup()

ag4 = loops[1].select('resnum ' + ' '.join(str(x) for x in range(92, 112 + 1))).toAtomGroup()

ag5 = target.select('resnum ' + ' '.join(str(x) for x in range(130 + 1, 142))).toAtomGroup()

ag6 = loops[2].select('resnum ' + ' '.join(str(x) for x in range(201, 217 + 1))).toAtomGroup()

ag7 = target.select('resnum ' + ' '.join(str(x) for x in range(81 + 1, 102))).toAtomGroup()

ag8 = target.select('not protein').toAtomGroup()
ags = [ag1, ag2, ag3, ag4, ag5, ag6, ag7,  ag8]
ag_all = prody_ext.combine_ags_into_one_chain(ags, 'merge')

#[pr.writePDB(workdir + 'ag_' + str(i) + '.pdb', ags[i])  for i in range(len(ags))]

pr.writePDB(workdir + 'merge.pdb', ag_all)

#ag_all_ala = prody_ext.target_to_all_gly_ala(ag_all, 'merge_allala', keep_no_protein=True)