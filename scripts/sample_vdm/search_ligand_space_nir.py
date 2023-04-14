'''
Once we anchor a ligand with some defined interactions.
We can generate ligand ensembels for vdm sampling.
The ligands ensembels can be generated with MOE or Schrodinger-Maestro.
'''

import os
import prody as pr

'''
python /mnt/e/GitHub_Design/Rvrsprot/scripts/sample_vdm/search_ligand_space_nir.py
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
workdir = '/mnt/e/DesignData/Chemodrugs/NIR/'

ligands = load_ligs(workdir + 'maestro/')

sel = 'name C12 C13 C14 C15 C16 C17 C18 C19 N2 N3'
rmsd_mean = []
for i in range(30):
    lig = ligands[i]
    rmsds = []
    for j in range(30):
        lig2 = ligands[j]
        rmsd = pr.calcRMSD(lig.select(sel), lig2.select(sel))
        rmsds.append(rmsd)
    rmsd_mean.append(sum(rmsds)/30)

min_value = min(rmsd_mean)
min_index = rmsd_mean.index(min_value)
print(min_index)
ligands[min_index].getTitle()

# middle lig is '24958200_3.pdb'

########################################################################################
#superimpost maestro ligs to designed lig on middle lig.
workdir = '/mnt/e/DesignData/Chemodrugs/NIR/'

outdir = workdir + 'ligs_1_tf2middle/'
os.makedirs(outdir, exist_ok= True)

design_lig = pr.parsePDB(workdir + 'lig_1st_round.pdb')
middle_lig = pr.parsePDB(workdir + 'maestro/24958200_3.pdb')
tf = pr.calcTransformation(middle_lig.select('name N2 C15 C16'), design_lig.select('name  C4 C6 C10'))

ligands = load_ligs(workdir + 'Maestro/')

for lig in ligands:
    tf.apply(lig)
    pr.writePDB(outdir + lig.getTitle(), lig)

########################################################################################