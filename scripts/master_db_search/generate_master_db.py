### create file list

import os

#Create PDB list file
workdir = '/mnt/e/DesignData/Metalloprotein/DPP/6HelixLarge/'

with open(workdir + 'list.txt', 'w') as f:
    for file in os.listdir(workdir):
        if not '.pdb' in file:
            continue
        f.write(workdir + file + '\n')

#Master createPDS
#./createPDS --type target --pdbList /mnt/e/DesignData/Metalloprotein/DPP/6HelixLarge/list.txt

#Create PDS list file
workdir = '/mnt/e/DesignData/Metalloprotein/DPP/6HelixLargeDB/'

with open(workdir + 'list', 'w') as f:
    for file in os.listdir(workdir):
        if not '.pds' in file:
            continue
        f.write(workdir + file + '\n')


# Create query pds
#./createPDS --type query --pdb /mnt/e/DesignData/Metalloprotein/DPP/fragment.pdb

#Master searchPDS
#./master --query /mnt/e/DesignData/Metalloprotein/DPP/fragment.pds --targetList /mnt/e/DesignData/Metalloprotein/DPP/6HelixLargeDB/list --rmsdCut 0.6 --matchOut /mnt/e/DesignData/Metalloprotein/DPP/out_6hL/matchOut.txt

from rvrsprot.external import query 
outdir = '/mnt/e/DesignData/Metalloprotein/DPP/out_6hL/'
os.makedirs(outdir, exist_ok=True)
pdb_path = '/mnt/e/DesignData/Metalloprotein/DPP/fragment.pdb'
targetList = '/mnt/e/DesignData/Metalloprotein/DPP/6HelixLargeDB/list'
rmsdCut = 0.5
query.master_query(outdir, pdb_path, targetList, rmsdCut)