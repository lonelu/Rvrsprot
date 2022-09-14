'''
Kehan created ~100 6 Helix bundles.
I previously use reverse vdm method generate C3 helix metal motifs. 
It is possible to use master search to search each the motif with the 6 helix bundles.
To do so, I will first createPDS of the 6 Helix bundles.
'''
import os
import rvrsprot.external.query as query
'''
python /mnt/e/GitHub_Design/Rvrsprot/scripts_tm/TM_6HelixBundle_MaterSearch.py
'''

workdir = '/mnt/e/DesignData/tm/Kehan/20220804/'

'''
#Generate pdblist.txt
with open(workdir + 'queries/pdblist.txt', 'w') as f:
    for file in os.listdir(workdir + 'queries/'):
        if not '.pdb' in file:
            continue

        f.write(workdir + 'queries/' + file + '\n')

#Create PDS
outdir = workdir  + 'queries_pds/'
os.makedirs(outdir, exist_ok=True)
for file in os.listdir(workdir + 'queries/'):
    if not '.pdb' in file:
        continue
    query.create_target_pds(workdir + 'queries/' + file, outdir)

#Generate pdslist.txt
with open(workdir + 'queries_pds/pdslist.txt', 'w') as f:
    for file in os.listdir(workdir + 'queries_pds/'):
        if not '.pds' in file:
            continue

        f.write(workdir + 'queries_pds/' + file + '\n')

'''
#MasterSearch
targetlist = workdir + 'queries_pds/pdslist.txt'
c3indir = '/mnt/e/DesignData/ReverseDesign/Phosphate/rvrs/'
c3outdir = workdir + 'phosph_output/'
os.makedirs(c3outdir, exist_ok=True)
for file in os.listdir(c3indir):
    if not 'pdb' in file:
        continue
    query.master_query(c3outdir, c3indir + file, targetlist, rmsdCut=1.2)


