import os

workdir = '/mnt/e/DesignData/Database/vdMh2o/'
workdir_vdM = workdir + 'vdm_pdb/'
outdir_vdM = workdir + 'vdm_pdb_h/'
os.makedirs(outdir_vdM, exist_ok=True)

for pdb in os.listdir(workdir_vdM):
    pdbpath = workdir_vdM + pdb
    os.system('reduce -BUILD -FLIP -Quiet {} > {}'.format( pdbpath, outdir_vdM + pdb))