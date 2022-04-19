import os

'''
python /mnt/e/GitHub_Design/Rvrsprot/scripts_tm/CalcPropensity.py
'''

def read_zs(workdir, file): 

    with open(workdir + file, 'r') as f:
        xxx = f.readline()

    xs = [float(x) for x in xxx.split(', ') if x!='']

    return xs

def bin_num(start, end, step, input):
    bins = {}

    x = start
    key = 0
    while x < end:
        if x + step <= end:
            bins[key]= [(x, x + step), 0]
        else:
            bins[key] = [(x, end), 0]     
        key += 1 
        x += step

    for k in bins.keys():
        r = bins[k][0]
        for inp in input:
            if r[0] <= inp and inp < r[1] :
                bins[k][1] += 1
    return bins


def write_bin_zs(workdir, bins, outfile):
    with open(workdir + outfile, 'w') as f:
        f.write('k\tbin\tcount\n')
        for b in bins.keys():
            f.write(str(b) + '\t(' + str(bins[b][0][0]) + ',' + str(bins[b][0][1]) +  ')\t'  + str(bins[b][1]) + '\n')


if __name__=='__main__':
    workdir = '/mnt/e/DesignData/tm/Kehan/_temp/'
    file = 'ALA.txt'
    start = -15
    end = 15
    step = 2
    zs = read_zs(workdir, file)
    bins = bin_num(start, end, step, zs)
    write_bin_zs(workdir, bins, outfile = 'zs_bin.tsv')