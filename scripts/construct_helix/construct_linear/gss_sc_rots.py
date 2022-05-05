from matplotlib.pyplot import contour
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

### Generate vdm rotate with different angles.

def _generate_sc_rots(outdir, std, std_cp_sel, tf_sel, tag_pre):
    dist = pr.calcDistance(std.select('serial 1'), std.select('serial 2'))[0]
    xcoord = np.array([[0, 0, 0], [dist , 0, 0]])
    tf = pr.calcTransformation(std.select(tf_sel).getCoords(), xcoord)
    tf_rv = pr.calcTransformation(xcoord, std.select(tf_sel).getCoords())
    
    for i in range(120):
        std_cp = std.select(std_cp_sel).copy()
        tf.apply(std_cp)
        v = std_cp.getCoords()
        theta = 2*math.pi/120 * i
        axis = xcoord[1]/norm(xcoord[1])
        rot = Rotation.from_rotvec(theta * axis)
        new_v = rot.apply(v) 
        new_v = tf_rv.apply(new_v) 
        std_cp.setCoords(new_v)
        pr.writePDB(outdir + tag_pre + str(i) +'.pdb', std_cp)
    return

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2/'
outdir = workdir + 'std_rots/'
os.makedirs(outdir, exist_ok = True)

std = pr.parsePDB(workdir + 'LinearHistidinesforLei.pdb')
#std.getSerials()

_generate_sc_rots(outdir, std.copy(), 'serial 2 3 4 5 7', 'serial 1 2', 'A_')
_generate_sc_rots(outdir, std.copy(), 'serial 13 14 15 16 18', 'serial 1 2', 'B_')

