import os
import pandas as pd
import shutil

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2/'

indir = workdir + '_z_fix_c2_enriched/'
outdir = workdir + '_z_fix_c2_enriched_cccp/'
os.makedirs(outdir, exist_ok=True)

parinfo = pd.read_csv(indir + '_summary_par_info.txt', sep='\t')

parinfo_x = parinfo[(parinfo['error'] < 1.0)]

for i in range(parinfo_x.shape[1]):
    file = parinfo_x.iloc[i][0][0:-7] + 'pdb'
    shutil.copy(indir + file, outdir + file)