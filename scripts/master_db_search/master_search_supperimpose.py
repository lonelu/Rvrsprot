import os
import prody as pr
import re

workdir = '/mnt/e/DesignData/Metalloprotein/DPP/out_6h/'

file = 'matchOut.txt'

match_infos = []

with open(workdir + file, 'r') as f:
    for line in f.readlines():
        xs = line.split('/')[-1]
        #print(xs)
        name = xs.split(' ')[0][:-4]
        number = re.findall(r'\d+', xs.split('[')[1])
        match_infos.append((name, number))

query = pr.parsePDB('/mnt/e/DesignData/Metalloprotein/DPP/fragement_lig.pdb')

targets = {}
db_dir = '/mnt/e/DesignData/Metalloprotein/DPP/6Helix/'
for file in os.listdir(db_dir):
    if not '.pdb' in file:
        continue
    targets[file[:-4]] = pr.parsePDB(db_dir + file)

for i in range(100, 150):
    name = match_infos[i][0]
    number = match_infos[i][1]
    target_sel = 'name CA and resindex ' + str(number[0]) + 'to' + str(number[1]) + ' '+ str(number[2]) + 'to' + str(number[3])
    query_cp = query.copy()
    pr.calcTransformation(query.select('name CA'), targets[name].select(target_sel)).apply(query_cp)
    pr.writePDB(workdir + name + '_' +str(i) + '_query.pdb', query_cp)
