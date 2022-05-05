'''

'''

import prody as pr
import os
from scipy.spatial.transform import Rotation
from numpy.linalg import norm
import math
from itertools import permutations
import numpy as np
from sklearn.neighbors import NearestNeighbors
from metalprot.database import database_cluster
from metalprot.basic import prody_ext
from metalprot.basic import vdmer
import shutil
import sys
sys.path.append(r'/mnt/e/GitHub_Design/Rvrsprot/scripts/construct_helix/construct_linear')
from gss_rvdmH_basics import clash, clash_rvs, cal_angle, calc_z_direction, rot_rvdmH_C2

'''
python /mnt/e/GitHub_Design/Rvrsprot/scripts/construct_helix/construct_linear/gss_zfix_c2.py
'''

def angle_filter(zpdb, z_sel, rpdb, r_sel, angle_min = 0, angle_max = 150):
        a_coords = zpdb.select(z_sel).getCoords()
        b_coords = rpdb.select(r_sel).getCoords()
        ae = cal_angle(a_coords[0] - a_coords[1], b_coords[0] - b_coords[1])
        if ae < angle_min or ae > angle_max:
            return False
        return True


def generate_rvdm2C2(z_fix_dir, rvdmH_dir, lig, outdir):
    '''
    Each vdmH, use the std_rot key (such as 'A_0') in the title to find the candidate rvdmH.
    Then combine them if not clash. 
    Then rot to generate C2.
    '''
    os.makedirs(outdir, exist_ok=True)
    for z in os.listdir(z_fix_dir):
        if not '.pdb' in z:
            continue
        zpdb = pr.parsePDB(z_fix_dir + z)
        folder = z.split('_')[5] + '_' + z.split('_')[6]
        for r in os.listdir(rvdmH_dir + folder + '/'):
            if not '.pdb' in r:
                continue
            rpdb = pr.parsePDB(rvdmH_dir + folder + '/' + r)
            #The distance between bb atoms from two parallel or antiparallel helixs is an important factor to consider.
            if clash(zpdb.select('bb').getCoords(), rpdb.select('bb').getCoords(), 5.0):
                continue
            
            if not angle_filter(zpdb, 'resnum 1 15 and name CA', rpdb,'resnum 1 15 and name CA', 150, 180):
                continue

            title = zpdb.getTitle() + '_' + rpdb.getTitle()
            ag = prody_ext.combine_ags([zpdb, rpdb], title)

            ag = rot_rvdmH_C2(ag, ABCchids=['A', 'M', 'B', 'C', 'D'])

            if clash(ag.select('bb and chid B').getCoords(), ag.select('bb and chid C').getCoords(), 5.0):
                continue

            if clash_rvs(ag.select('bb and chid B').getCoords(), ag.select('bb and chid C').getCoords(), 10.0):
                continue
            
            if not angle_filter(ag, 'chid A and resnum 1 15 and name CA', ag,'chid C and resnum 1 15 and name CA', 0, 30):
                continue
            
            #Clashing ligand. The rotation of the ligand is not considered here. 
            if clash(lig.select('name C2 C3 C4 C9 C10 C11').getCoords(), ag.select('bb and chid B D').getCoords(), 1.8):
                continue
            pr.writePDB(outdir + title + '.pdb', ag)
    return 


workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2/'
z_fix_dir = workdir + '_z_fix_H/'
rvdmH_dir = workdir + 'std_rots_rvdm/'
outdir = workdir + '_z_fix_c2/'
lig = pr.parsePDB(workdir + 'Ligand.pdb')
generate_rvdm2C2(z_fix_dir, rvdmH_dir, lig, outdir)