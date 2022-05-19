'''
After master searching the loops. 
Summarize the count of loops for each structure.
'''

import os

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c3/c3_vc0/loop2/'

info = {}
for folder in os.listdir(workdir):

    with open(workdir + folder + '/ABC_frt/A-9-2-B-8-2/seq.txt') as f:
        count = len(f.readlines())
    with open(workdir + folder + '/ABC_frt/A-9-2-C-8-2/seq.txt') as f:
        count2 = len(f.readlines())
    info[folder] = (count, count2)
