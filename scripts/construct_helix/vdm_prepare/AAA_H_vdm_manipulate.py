'''
Change the name of vdms.
'''

import os
import shutil

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/SamM/'

indir = workdir + 'porphyrin_his_vdm0/'

indir2 = '/mnt/e/DesignData/ligands/ZN_rcsb_datesplit/20211013/20211015_AAext3/AAA_H_vdms_original/'

outdir = workdir + 'porphyrin_his_vdm/'
os.makedirs(outdir, exist_ok=True)

clu_name_dict = {}
for file in os.listdir(indir2):
    if not '.pdb' in file:
        continue
    clu_name = '_'.join(file.split('_')[0:3])
    clu_name_dict[clu_name] = '_'.join([file.split('_')[7], file.split('_')[-2], file.split('_')[-1]])


for file in os.listdir(indir):
    if not '.pdb' in file:
        continue
    clu_name = '_'.join(file.split('_')[0:3])
    file_new = file[0:-4] + '_' + clu_name_dict[clu_name]
    shutil.copy(indir + file, outdir + file_new)
