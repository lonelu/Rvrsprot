'''
For all the designed tm proteins, we would summary the count to see which topology&position is most designable. 
'''
import os
import statistics

workdir = '/mnt/e/DesignData/tm/Kehan/output/output_sd2/'

topo_dict = {}

for _d1 in os.listdir(workdir):
    if not os.path.isdir(workdir + _d1):
        continue
    for _f1 in os.listdir(workdir + _d1):
        if not 'summary' in _f1:
            continue
        with open(workdir + _d1 + '/' + _f1, 'r') as f:
            lines = f.readlines()
            ks = set() # In the same file, we only count each topo once.
            for l in lines[1:]:
                sps = l.split('\t')
                info = (float(sps[5]), float(sps[6]), float(sps[11]))
                key = (sps[0], sps[1])
                if key in ks:
                    continue
                ks.add(key)
                if key in topo_dict:
                    topo_dict[key].append(info)
                else:
                    topo_dict[key] = [info]

with open(workdir + 'metalprot_summary.tsv', 'w') as f:
    f.write('topo\tcomp\tcount\tmean_cluscore\tmean_overscore\tmean_overnum\n')
    for k in topo_dict.keys():
        num = len(topo_dict[k])
        mean_cluscore = round( statistics.mean([topo_dict[k][i][0] for i in range(num)]), 2)
        mean_overscore = round( statistics.mean([topo_dict[k][i][1] for i in range(num)]), 2)
        mean_overnum = round( statistics.mean([topo_dict[k][i][2] for i in range(num)]), 2)
        f.write(k[0] + '\t' + k[1] + '\t' + str(num) + '\t' + str(mean_cluscore)  + '\t' + str(mean_overscore) + '\t' + str(mean_overnum) + '\n')