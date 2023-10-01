'''
Once we anchor a ligand with some defined interactions.
We can generate ligand ensembels for vdm sampling.
The ligands ensembels can be generated with MOE or Schrodinger-Maestro.
'''

import os
import prody as pr

'''
python /mnt/e/GitHub_Design/Rvrsprot/scripts/sample_vdm/search_ligand_space.py
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
#superimpost maestro ligs to designed lig on indole group.
workdir = '/mnt/e/DesignData/Chemodrugs/HB_RUC__ligs/'

outdir = workdir + 'ligs_1st_tf2indole/'
os.makedirs(outdir, exist_ok= True)

design_lig = pr.parsePDB(workdir + 'lig_1st_round.pdb')

ligands = load_ligs(workdir + 'Maestro/')

for lig in ligands:
    tf = pr.calcTransformation(lig.select('name  N1 C1 C2 C4 C5 C6 C10 C11 C12'), design_lig.select('name  N1 C1 C2 C4 C5 C6 C10 C11 C12')).apply(lig)
    #rmsd = pr.calcRMSD(lig.select('name N2 C7'), design_lig.select('name N2 C7'))
    #if rmsd > 1.0:
    #    continue
    pr.writePDB(outdir + lig.getTitle(), lig)


########################################################################################
#Get middle maestro lig by rmsd. The middle maestro lig should be the one have minium sum rmsd with all others.
workdir = '/mnt/e/DesignData/Chemodrugs/HB_RUC__ligs/'

ligands = load_ligs(workdir + 'Maestro/')

sel = 'name C8 C4 N1 C5 C2 C1 C3 C6 C9 C11 C12 C10 F1'
rmsd_mean = []
for i in range(50):
    lig = ligands[i]
    rmsds = []
    for j in range(50):
        lig2 = ligands[j]
        rmsd = pr.calcRMSD(lig.select(sel), lig2.select(sel))
        rmsds.append(rmsd)
    rmsd_mean.append(sum(rmsds)/50)

min_value = min(rmsd_mean)
min_index = rmsd_mean.index(min_value)
print(min_index)
ligands[min_index].getTitle()

# middle lig is '9931954_24.pdb'


########################################################################################
#superimpost maestro ligs to designed lig on middle lig.
workdir = '/mnt/e/DesignData/Chemodrugs/HB_RUC__ligs/'

outdir = workdir + 'ligs_1st_tf2middle/'
os.makedirs(outdir, exist_ok= True)

design_lig = pr.parsePDB(workdir + 'lig_2nd_round.pdb')
middle_lig = pr.parsePDB(workdir + 'Maestro/9931954_24.pdb')
tf = pr.calcTransformation(middle_lig.select('name  N1 C1 C2 C4 C5 C6 C10 C11 C12'), design_lig.select('name  N1 C1 C2 C4 C5 C6 C10 C11 C12'))

ligands = load_ligs(workdir + 'Maestro/')

for lig in ligands:
    tf.apply(lig)
    pr.writePDB(outdir + lig.getTitle(), lig)

########################################################################################
#superimpost maestro ligs to designed lig on middle lig.
workdir = '/Users/lonelu/DesignData/Chemodrugs/HB_RUC__ligs/'

outdir = workdir + 'ligs_2nd_tf2middle/'
os.makedirs(outdir, exist_ok= True)

design_lig = pr.parsePDB(workdir + 'lig_2nd_round.pdb')
middle_lig = pr.parsePDB(workdir + 'Maestro/9931954_24.pdb')
tf = pr.calcTransformation(middle_lig.select('name C4 C6 C10'), design_lig.select('name C4 C6 C10'))

ligands = load_ligs(workdir + 'Maestro/')

for lig in ligands:
    tf.apply(lig)
    pr.writePDB(outdir + lig.getTitle(), lig)

#########################################################################################superimpost maestro ligs to designed lig on middle lig.
workdir = '/Users/lonelu/DesignData/Chemodrugs/HB_RUC__ligs/'

outdir = workdir + 'ligs_rucA_tf2middle/'
os.makedirs(outdir, exist_ok= True)

design_lig = pr.parsePDB(workdir + 'ruc_A_lig.pdb')
middle_lig = pr.parsePDB(workdir + 'Maestro/9931954_24.pdb')
tf = pr.calcTransformation(middle_lig.select('name C4 C6 C10'), design_lig.select('name C8 C21 C24'))

ligands = load_ligs(workdir + 'Maestro/')

for lig in ligands:
    tf.apply(lig)
    pr.writePDB(outdir + lig.getTitle(), lig)

########################################################################################