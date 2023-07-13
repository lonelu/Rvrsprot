'''
Once we anchor a ligand with some defined interactions.
We can generate ligand ensembels for vdm sampling.
The ligands ensembels can be generated with MOE or Schrodinger-Maestro.
'''

import os
import prody as pr

'''
python /mnt/e/GitHub_Design/Rvrsprot/scripts/sample_vdm/search_ligand_space_mef.py
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
workdir = '/Users/lonelu/DesignData/Chemodrugs/MEF/'

ligands = load_ligs(workdir + 'maestro/')

sel = 'name C1 C2 C3 C5 C6 C10 C14 C15 C16 O1'
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

# middle lig is '71488522_38.pdb'

########################################################################################
#superimpost maestro ligs to designed lig on middle lig.
workdir = '/Users/lonelu/DesignData/Chemodrugs/MEF/'

outdir = workdir + 'ligs_2nd_tf2middle/'
os.makedirs(outdir, exist_ok= True)

design_lig = pr.parsePDB(workdir + 'lig_2nd_round.pdb')
middle_lig = pr.parsePDB(workdir + 'maestro/71488522_38.pdb')
tf = pr.calcTransformation(middle_lig.select('name C2 C6 C10'), design_lig.select('name  C4 C6 C10'))

ligands = load_ligs(workdir + 'Maestro/')

for lig in ligands:
    tf.apply(lig)
    pr.writePDB(outdir + lig.getTitle(), lig)

########################################################################################