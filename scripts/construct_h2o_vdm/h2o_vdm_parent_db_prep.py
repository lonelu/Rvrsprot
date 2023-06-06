import os
import prody as pr
 

# Get pdb chid list. The file is from vdM paper S1.
pdb_name_file = '/mnt/e/GitHub_Design/Rvrsprot/scripts/construct_h2o_vdm/abb8330_data_s1.txt'

pdb_chid = {} 
with open(pdb_name_file, 'r') as f:
    lines = f.readlines()
    for line in lines:
        if '#' in line:
            continue
        pdb = line[:4]
        if pdb in pdb_chid.keys():
            pdb_chid[pdb].append(line[4])
        else:
            pdb_chid[pdb] = [line[4]]

# Download pdb database.

workdir = '/mnt/e/DesignData/Database/vdMh2o/'
workdir_ori = workdir + 'parent_pdbs/'
os.makedirs(workdir_ori, exist_ok=True)

exist_pdb = set()
for file in os.listdir(workdir_ori):
    if file.endswith(".pdb.gz"):
        exist_pdb.add(file.split('.')[0].upper())

pr.pathPDBFolder(workdir_ori)
for p in pdb_chid.keys():
    if p in exist_pdb: continue
    if p.upper() in exist_pdb: continue
    pr.fetchPDBviaFTP(p.upper(), compressed = False) 


# Extract chid and h2o.
workdir_ext = workdir + 'parent_pdb_exts/'
os.makedirs(workdir_ext, exist_ok=True)

for pdb_name in pdb_chid.keys():
    if pdb_name.upper() not in exist_pdb: continue
    pdb = pr.parsePDB(workdir_ori + pdb_name.upper() + '.pdb.gz')
    chid = pdb_chid[pdb_name]
    pdb_ext = pdb.select('water or protein and chid ' + ' '.join(c for c in chid))
    pr.writePDB(workdir_ext + pdb_name + '.pdb', pdb_ext)



