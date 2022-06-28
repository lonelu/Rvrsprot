'''
Phosphate is tetrahydral, which is similar to Metal binding tetrahydral geometry.
Similarly, we can use the same strategy to reversely generate C3 for 6 helix bundles.
Treat each O-P-O bond as cg 'COO'.
'''

import os
import sys
import prody as pr
import numpy as np

from metalprot.basic import utils, transformation, prody_ext
from metalprot.combs import gvdm_helper, search_lig_indep, search_lig_indep_wrap, search_lig_indep_inpair

from rvrsprot.vdm_rev import reverse_vdm

sys.path.append(r'/mnt/e/GitHub_Design/Rvrsprot/scripts/construct_helix/construct_linear')
from gss_rvdmH_basics import clash, clash_rvs, cal_angle, calc_z_direction, rot_rvdmH_C2

path_to_database='/mnt/e/DesignData/Database/rvdMs/'

workdir = '/mnt/e/DesignData/ReverseDesign/Phosphate/'
outdir = workdir + 'rvrs/'

helix_path =  '/mnt/e/DesignData/ReverseDesign/c2_coil/CC_15-Mer.pdb'
helix = pr.parsePDB(helix_path)
helix_resnum = ('A', 7)

#>>> Manipulate phosphate.
# lig = pr.parsePDB(workdir + 'phosphate.pdb')
# bond_dist = pr.calcDistance(lig.select('name P')[0], lig.select('name O3')[0])
# z_axis = np.array([[0, 0, 0], [0, 0, bond_dist]])

# pr.calcTransformation(lig.select('name P O3').getCoords(), z_axis).apply(lig)
# pr.writePDB(workdir + 'phosphate_std.pdb', lig)

lig = pr.parsePDB(workdir + 'phosphate_std.pdb')


def generate_phosphate_c3(outdir, cg, aa, correspond_resname, represent_name, vdm_cg_sel, lig_sel):
    os.makedirs(outdir, exist_ok=True)
    df_vdm = gvdm_helper.load_old_vdm(path_to_database, cg, aa)
    df_vdm_f = gvdm_helper.filter_db(df_vdm, use_enriched = True, use_abple=True, abple='A')

    df_vdm_rep = df_vdm_f[(df_vdm_f['chain'] == 'Y') & (df_vdm_f['resname'] == correspond_resname) & (df_vdm_f['name'] == represent_name)][['CG', 'rota', 'probe_name']]

    vdm_pd = df_vdm[(df_vdm['CG'] == df_vdm_rep.iloc[0]['CG']) 
                    & (df_vdm['rota'] == df_vdm_rep.iloc[0]['rota']) 
                    & (df_vdm['probe_name'] == df_vdm_rep.iloc[0]['probe_name'])]
    vdm0 = gvdm_helper.df2ag(vdm_pd)

    std_coords = np.array([lig.select('name ' + x).getCoords()[0] for x in lig_sel])

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

        if v.empty:
            print('vdM is None')
            continue
        if not v[['contact_hb']].any()[0]: # or v[['ABPLE_3mer']][0] != 'AAA':
            continue    

        ag = reverse_vdm.vdm_add_helix(df_vdm_f, _cg, _rota, _probe_name, helix, helix_resnum, check_clash = True)
        if not ag:
            continue
        if clash(np.zeros((1, 3), dtype=float), ag.select('bb').getCoords(), 5.0):
            continue
        angle, angle_zy, angle_zx = calc_z_direction(ag, sel = 'resindex 0 14 and name CA')
        #print('angle: {}, angle_zx: {}'.format(angle, angle_zx))
        if ( angle < 30 and (angle_zx <20 or angle_zx > 160)) or (angle > 150 and (angle_zx <20 or angle_zx > 160)):
            ag1 = ag.copy()
            pr.calcTransformation(lig.select('name P O0 O3'), lig.select('name P O1 O3')).apply(ag1)
            if clash(ag.select('bb').getCoords(), ag1.select('bb').getCoords(), 5.0):
                continue
            ag2 = ag.copy()
            pr.calcTransformation(lig.select('name P O0 O3'), lig.select('name P O2 O3')).apply(ag2)
            ag_all = prody_ext.combine_ags([ag, ag1, ag2], ag.getTitle(), ['A', 'B', 'C'])
            pr.writePDB(outdir + aa + '_' + ag.getTitle(), ag_all)
    return
    

cg = 'coo'
aa = 'SER'
correspond_resname = 'GLU'
represent_name = 'CG'
vdm_cg_sel = 'chid Y and name CG OE1 OE2'
lig_sel = ['P', 'O0', 'O3']
outdir = workdir + 'rvrs0/'
generate_phosphate_c3(cg, aa, correspond_resname, represent_name, vdm_cg_sel, lig_sel)
        

cg = 'coo'
aa = 'SER'
correspond_resname = 'GLU'
represent_name = 'CG'
vdm_cg_sel = 'chid Y and name CG OE1 OE2'
lig_sel = ['P', 'O3', 'O0']
outdir = workdir + 'rvrs1/'
generate_phosphate_c3(cg, aa, correspond_resname, represent_name, vdm_cg_sel, lig_sel)
     

