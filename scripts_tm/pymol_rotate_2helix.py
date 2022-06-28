'''
Using pymol to rotate 2 helix bundles and build 6 helix bundles.
'''

from pymol import cmd
import glob
import os
import sys
# sys.path.append(r'/mnt/e/DesignData/tm/Kehan/20220602_2helix/')
# from flatten_obj import flatten_obj

indir = '/mnt/e/DesignData/tm/Kehan/20220602_2helix/'
outdir = indir + 'out/'
os.makedirs(outdir, exist_ok=True)

for _file in glob.glob(os.path.join(indir, '*.pdb')):
    cmd.load(_file)
seleobjs = cmd.get_names()

obj0 = seleobjs[0]
cmd.alter('model ' + obj0, 'seg = ""')

cmd.translate([0, -7, 0], object = obj0)
cmd.save(outdir + 'obj0.pdb', obj0)

cmd.copy('obj1', obj0)
cmd.copy('obj2', obj0)
cmd.alter('model obj1 and chain A', 'chain = "C"')
cmd.alter('model obj1 and chain B', 'chain = "D"')
cmd.alter('model obj2 and chain A', 'chain = "E"')
cmd.alter('model obj2 and chain B', 'chain = "F"')
cmd.rotate(axis = 'z', angle = 120, object = 'obj1')
cmd.rotate(axis = 'z', angle = 240, object = 'obj2')
cmd.save(outdir + 'obj1.pdb', 'obj1')
cmd.save(outdir + 'obj2.pdb', 'obj2')


cmd.save(outdir + 'out6helix.pdb', obj0 + ' or obj1 or obj2')
#cmd.delete('all')

