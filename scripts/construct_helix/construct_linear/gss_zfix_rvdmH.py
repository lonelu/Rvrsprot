'''
Add helix to rvdms for cg hid of GLU, ASP, SER, THR. 
Superimpose on std motif rotations.
Keep the rvdmH is z-fixed. 
z-fix is the helix is in the same direction with z-axis.
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

from metalprot.combs import search_lig_indep
from metalprot.basic import transformation
from metalprot.basic import prody_ext
import numpy as np
import prody as pr
import pandas as pd
from rvrsprot.combs import reverse_vdm
import sys
sys.path.append(r'/mnt/e/GitHub_Design/Rvrsprot/scripts/construct_helix/construct_linear')
from gss_rvdmH_basics import clash, clash_rvs, cal_angle, calc_z_direction, rot_rvdmH_C2


def filter_rvdm_for_special_his(path_to_database, cg, aa, helix, helix_resind, vdm_cg_sel, std_cp, outdir):
    '''
    Here we only want to keep rvdms with direction agree z-axis.
    '''
    df_vdm = search_lig_indep.load_old_vdm(path_to_database, cg, aa)
    df_vdm_f = search_lig_indep.filter_db(df_vdm, use_enriched = True, use_abple = True, abple = 'A')
    df_vdm_rep = df_vdm_f[['CG', 'rota', 'probe_name']].drop_duplicates()

    vdm_pd = df_vdm[(df_vdm['CG'] == df_vdm_rep.iloc[0]['CG']) 
                    & (df_vdm['rota'] == df_vdm_rep.iloc[0]['rota']) 
                    & (df_vdm['probe_name'] == df_vdm_rep.iloc[0]['probe_name'])]
    vdm0 = search_lig_indep.df2ag(vdm_pd)

    std_coords = np.array([std_cp.select('serial ' + str(i)).getCoords()[0] for i in [4, 5, 2, 3, 1]])


    a_coods = vdm0.select(vdm_cg_sel).getCoords()
    b_coods = std_coords
    coords4change = df_vdm_f[['c_x', 'c_y', 'c_z']].values
    R, m_com, t_com = transformation.get_rot_trans(a_coods, b_coods)
    df_vdm_f.loc[:, ['c_x', 'c_y', 'c_z']] = np.dot(coords4change - m_com, R) + t_com

    for i in range(0, df_vdm_rep.shape[0]):
        _cg = df_vdm_rep.iloc[i]['CG']
        _rota = df_vdm_rep.iloc[i]['rota']
        _probe_name = df_vdm_rep.iloc[i]['probe_name']
        v = df_vdm_f[(df_vdm_f['CG'] == _cg) & (df_vdm_f['rota'] == _rota) & (df_vdm_f['probe_name'] == _probe_name)]
        #Filter hbond. Filter ABPLE_3mer.
        if not v[['contact_hb']].any()[0]: # or v[['ABPLE_3mer']][0] != 'AAA':
            continue
        ag = reverse_vdm.vdm_add_helix(df_vdm_f, _cg, _rota, _probe_name, helix, helix_resind)
        if ag:
            if clash(np.zeros((1, 3), dtype=float), ag.getCoords(), 3.0):
                continue
            angle, angle_zy, angle_zx = calc_z_direction(ag, sel = 'resindex 0 14 and name CA')
            #print('angle: {}, angle_zx: {}'.format(angle, angle_zx))
            if ( angle < 30 and (angle_zx <20 or angle_zx > 160)) or (angle > 150 and (angle_zx <20 or angle_zx > 160)):
                pr.writePDB(outdir + aa + '_' + ag.getTitle(), ag)
    return


def generate_std_rots_rvdm(outdir, std_rots_dir, path_to_database, cg, aa, helix_path, helix_resind, vdm_cg_sel):
    helix = pr.parsePDB(helix_path)
    os.makedirs(outdir, exist_ok=True)
    for x in os.listdir(std_rots_dir):
        if not '.pdb' in x or not 'A' in x:
            continue
        std_cp = pr.parsePDB(std_rots_dir + x)

        _outdir = outdir + x.split('.')[0] + '/'
        os.makedirs(_outdir, exist_ok=True)

        filter_rvdm_for_special_his(path_to_database, cg, aa, helix, helix_resind, vdm_cg_sel, std_cp, _outdir)
    return

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2_coil/'

path_to_database='/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/vdMs/'
cg = 'hid'

vdm_cg_sel = 'chid Y and name CG ND1 CD2 CE1 NE2'
helix_path =  '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2_coil/CC_15-Mer.pdb'
helix_resind = 7

outdir = workdir + 'std_rots_rvdm/'
std_rots_dir = workdir + 'std_rots/'

aa = 'GLU'
generate_std_rots_rvdm(outdir, std_rots_dir, path_to_database, cg, aa, helix_path, helix_resind, vdm_cg_sel)

aa = 'ASP'
generate_std_rots_rvdm(outdir, std_rots_dir, path_to_database, cg, aa, helix_path, helix_resind, vdm_cg_sel)

aa = 'SER'
generate_std_rots_rvdm(outdir, std_rots_dir, path_to_database, cg, aa, helix_path, helix_resind, vdm_cg_sel)

aa = 'THR'
generate_std_rots_rvdm(outdir, std_rots_dir, path_to_database, cg, aa, helix_path, helix_resind, vdm_cg_sel)

'''
#Test a single case.
workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2/'

path_to_database='/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/vdMs/'
cg = 'hid'
aa = 'GLU'

helix_path =  '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2/15mer_ALA.pdb'
helix_resind = 7

vdm_cg_sel = 'chid Y and name CG ND1 CD2 CE1 NE2'
std_cp = pr.parsePDB(workdir + 'std_rots/A_0.pdb')

outdir = workdir + 'std_rots_rvdm/A_0/'
os.makedirs(outdir, exist_ok = True)

helix = pr.parsePDB(helix_path)
filter_rvdm_for_special_his(path_to_database, cg, aa, helix, helix_resind, vdm_cg_sel, std_cp, outdir)

'''