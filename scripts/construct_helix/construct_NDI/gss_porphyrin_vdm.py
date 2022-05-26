import os
import prody as pr
import numpy as np
from metalprot.basic import prody_ext
from sklearn.neighbors import NearestNeighbors

def generate_porphyrin_vdm(vdm_dir, outdir, std, porphyrin):
    '''
    Attach the porphyrin on the his vdm, which could be used for distance check later.
    '''
    os.makedirs(outdir, exist_ok=True)
    for file in os.listdir(vdm_dir):
        if not '.pdb' in file:
            continue
        v = pr.parsePDB(vdm_dir + file)
        std_coords = [std.select('serial ' + str(i)).getCoords()[0] for i in [4, 5, 2, 3, 1]]
        pr.calcTransformation(v.select('name CG ND1 CD2 CE1 NE2').getCoords(), np.array(std_coords)).apply(v)
        print(v.getTitle())
        ag = prody_ext.combine_ags([v.select('not name ZN'), v.select('name ZN'), porphyrin.copy()], title=v.getTitle(), ABCchids=['A', 'M', 'P'])
        pr.writePDB(outdir + file, ag)
    return

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/SamM/'

std = pr.parsePDB(workdir + 'A_0.pdb')
NDI = pr.parsePDB(workdir + 'NDI.pdb')
porphyrin = pr.parsePDB(workdir + 'Ph_porphyrin.pdb')

dist = pr.calcDistance(porphyrin.select('name MN1')[0], porphyrin.select('name C5')[0])
pr.calcTransformation(porphyrin.select('name C5 C10 MN1'), np.array([[0, 0, dist], [0, dist, 0], [0, 0, 0]])).apply(porphyrin)


vdm_dir = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211015_AAext3/AAA_H_vdms_name/'
outdir = workdir + 'porphyrin_his_vdm/'
generate_porphyrin_vdm(vdm_dir, outdir, std, porphyrin)

