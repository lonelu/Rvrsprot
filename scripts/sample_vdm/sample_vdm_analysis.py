'''

'''
import os
import prody as pr


'''
python /Users/lonelu/GitHub_Design/Rvrsprot/scripts/sample_vdm/sample_vdm_analysis.py
'''


workdir = '/Users/lonelu/DesignData/Chemodrugs/RUC_vdm/HB_RUC_rucA_vdm/vdm_sampling_res29_ASP_hidhie2/'

outdir = workdir + 'duplicate_remove/'
os.makedirs(outdir, exist_ok= True)

def load_ligs_no_rep(lig_path):
    keys = []
    ligs = []
    for lig_name in os.listdir(lig_path):
        if not '.pdb' in lig_name:
            continue
        if not 'key' in lig_name:
            continue
        if lig_name.split('__')[-1] in keys:
            continue
        keys.append(lig_name.split('__')[-1])
        lig = pr.parsePDB(lig_path + lig_name)
        ligs.append(lig)
    print(keys)
    return ligs

ligands = load_ligs_no_rep(workdir)

for lig in ligands:
    pr.writePDB(outdir + lig.getTitle(), lig)



