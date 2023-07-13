'''
Once we anchor a ligand with some defined interactions.
We can generate ligand ensembels for vdm sampling.
The ligands ensembels can be generated with MOE or Schrodinger-Maestro.
'''

import os
import prody as pr

'''
python /mnt/e/GitHub_Design/Rvrsprot/scripts/sample_vdm/search_ligand_space_vel.py
'''

########################################################################################
def load_ligs(lig_path):
    ligs = []
    for lig_name in os.listdir(lig_path):
        if not '.pdb' in lig_name:
            continue
        lig = pr.parsePDB(lig_path + lig_name)
        ligs.append(lig)
    return ligs

########################################################################################
#Get middle maestro lig by rmsd. The middle maestro lig should be the one have minium sum rmsd with all others.
workdir = '/Users/lonelu/DesignData/Chemodrugs/VEL/'

ligands = load_ligs(workdir + 'maestro/')

sel = 'name C5 C7 C8 C9 C10 C11 C12 C13 N2 N3'
rmsd_mean = []
for i in range(len(ligands)):
    lig = ligands[i]
    rmsds = []
    for j in range(len(ligands)):
        lig2 = ligands[j]
        rmsd = pr.calcRMSD(lig.select(sel), lig2.select(sel))
        rmsds.append(rmsd)
    rmsd_mean.append(sum(rmsds)/len(ligands))

min_value = min(rmsd_mean)
min_index = rmsd_mean.index(min_value)
print(min_index)
ligands[min_index].getTitle()

# middle lig is '11960529_12.pdb'

########################################################################################
#superimpost maestro ligs to designed lig on middle lig.
import numpy as np
workdir = '/Users/lonelu/DesignData/Chemodrugs/VEL/'

outdir = workdir + 'ligs_2nd_tf2middle/'
os.makedirs(outdir, exist_ok= True)

design_lig = pr.parsePDB(workdir + 'lig_2nd_round.pdb')
middle_lig = pr.parsePDB(workdir + 'maestro/11960529_12.pdb')
#y = np.array([middle_lig.select('name ' + x).getCoords()[0] for x in ['C5', 'C10', 'C9']])
#y = middle_lig.select('name C5 C9 C10').getCoords()
#x = design_lig.select('name C4 C6 C10').getCoords()
#tf = pr.calcTransformation(x, y)
tf = pr.calcTransformation(middle_lig.select('name C5 C9 C10'), design_lig.select('name  C4 C6 C10'))

ligands = load_ligs(workdir + 'Maestro/')

for lig in ligands:
    tf.apply(lig)
    pr.writePDB(outdir + lig.getTitle(), lig)

########################################################################################