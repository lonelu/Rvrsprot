'''
Kehan is preparing TM protein propensity along y-axis. 
He need to use ahull to calculate the position of aas. He need to mutate all 'G' to 'A' to use ahull with Rosetta. 
There are a few proteins Rosetta could not read.
'''
import os
import prody as pr
from metalprot.basic import prody_ext

in_dir  = '/mnt/e/DesignData/tm/database/tms_filter35_prots/'
out_dir = in_dir + 'mutate/'
os.makedirs(out_dir, exist_ok=True)

proteins = ['6orv.pdb','6p9x.pdb','5ogl.pdb','7dz7.pdb','4pd4.pdb','3rce.pdb','1nek.pdb','1dxr.pdb',
    '2cfq.pdb','2cfp.pdb','1h2s.pdb','5nmi.pdb','3nog.pdb','1fft.pdb','1ezv.pdb','7e0k.pdb']

for  p in proteins:
    prot = pr.parsePDB(in_dir + p).select('protein')

    win_to_mutate = prot.select('name CA and resname GLY').getResindices()
    title = p.split('.')[0] + '_mut'
    try:
        prot_mut = prody_ext.target_mutation(prot, title, win_to_mutate, aa = 'ALA')
        pr.writePDB(out_dir +  title, prot_mut)
    except:
        print(title)