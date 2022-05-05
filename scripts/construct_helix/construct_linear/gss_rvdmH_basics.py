'''
Basic funcitons for using rvdm to generate helix bundles.
'''

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


def clash(Acoords, Bcoords, vdm_vdm_clash_dist = 2.6):
    neigh_y = NearestNeighbors(radius= vdm_vdm_clash_dist)
    neigh_y.fit(Acoords)
    x_in_y = neigh_y.radius_neighbors(Bcoords)
    x_has_y = any([True if len(a) >0 else False for a in x_in_y[1]])
    if x_has_y:
        return True
    return False 


def clash_rvs(Acoords, Bcoords, vdm_vdm_clash_dist = 10.0):
    neigh_y = NearestNeighbors(radius= vdm_vdm_clash_dist)
    neigh_y.fit(Acoords)
    x_in_y = neigh_y.radius_neighbors(Bcoords)
    x_has_y = any([True if len(a) >0 else False for a in x_in_y[1]])
    if x_has_y:
        return False
    return True 


def cal_angle(a, b):
    unit_vector1 = a / np.linalg.norm(a)
    unit_vector2 = b / np.linalg.norm(b)

    dot_product = np.dot(unit_vector1, unit_vector2)
    angle = np.arccos(dot_product)*180/math.pi
    return angle    


def calc_z_direction(helix_std, sel = 'resindex 0 14 and name CA'):
    '''
    helix_std: std helix which is already superimposed to rotated imido.
    '''
    a_coords = helix_std.select(sel).getCoords()
    a = a_coords[0] - a_coords[1]
    b = [0, 0, 1] 
    b_zy = [0, a[1], a[2]] #Important here. Which angle does we want, the exact relative to z-axis or relative to z_y plate.
    b_zx = [a[0], 0, a[2]]

    angle = cal_angle(a, b) #angle in radian to degree
    angle_zy = cal_angle(a, b_zy) 
    angle_zx = cal_angle(a, b_zx) 

    return angle, angle_zy, angle_zx


def rot_rvdmH_C2(rvdmH, ABCchids = None):
    dist = pr.calcDistance(rvdmH.select('name NE2'), rvdmH.select('name ZN'))[0]
    xcoord = np.array([[dist , 0, 0], [0, 0, 0]])
    tf = pr.calcTransformation(rvdmH.select('name NE2 ZN').getCoords(), xcoord)

    tf.apply(rvdmH)
    rvdmH_cp = rvdmH.select('protein').copy()
    theta = math.pi
    axis = [0 , 0, dist]/norm([0 , 0, dist])
    rot = Rotation.from_rotvec(theta * axis)
    new_v = rot.apply(rvdmH_cp.getCoords()) 
    rvdmH_cp.setCoords(new_v)
    ag = prody_ext.combine_ags([rvdmH, rvdmH_cp], rvdmH.getTitle(), ABCchids)

    return ag

