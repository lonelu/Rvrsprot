'''
Generating 2 helix bundles for 2 cgs of belinostat.
The pair selected here is ['S', 'O1', 'O2'] and ['O3', 'C12', 'N2']

The current purpose is to create a seudo target for the 'run_pos_lig_by_pair_vdm.py' 
'''
import os
import prody as pr
import numpy as np

from metalprot.basic import utils, transformation
from metalprot.combs import gvdm_helper, search_lig_indep, search_lig_indep_wrap, search_lig_indep_inpair

from rvrsprot.vdm_rev import reverse_vdm

'''
python /mnt/e/GitHub_Design/Metalprot/scrips/position_ligand/belinostat/pose_lig_by_pair_vdm_dev.py
'''

path_to_database='/mnt/e/DesignData/Database/rvdMs/'

workdir = '/mnt/e/DesignData/Metalloenzyme/belinostat/'
outdir = workdir + 'rvrs/'
os.makedirs(outdir, exist_ok=True)

helix_path =  '/mnt/e/DesignData/ReverseDesign/c2_coil/CC_15-Mer.pdb'
helix = pr.parsePDB(helix_path)
helix_resnum = ('A', 7)
lig = pr.parsePDB(workdir + 'ligs/meo_50g_amber14eht_md_out/50g_md_0.pdb')


def add_one_helix(cg, aa, correspond_resname, represent_name, vdm_cg_sel, lig_sel):
    df_vdm = gvdm_helper.load_old_vdm(path_to_database, cg, aa)
    df_vdm_f = gvdm_helper.filter_db(df_vdm, use_enriched = True, use_abple=True, abple='A')

    #df_vdm_rep = df_vdm_f[['CG', 'rota', 'probe_name']].drop_duplicates()
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

    #for i in range(0, df_vdm_rep.shape[0]):
    for i in range(0, 10):
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
        pr.writePDB(outdir + aa + '_' + ag.getTitle(), ag)
        #break

    return

cg = 'coo'
aa = 'TYR'
correspond_resname = 'GLU'
represent_name = 'CD'
vdm_cg_sel = 'chid Y and name CD OE1 OE2'
lig_sel = ['S', 'O1', 'O2']
add_one_helix(cg, aa, correspond_resname, represent_name, vdm_cg_sel, lig_sel)

cg = 'conh2'
aa = 'GLU'
correspond_resname = 'ASN'
represent_name = 'CG'
vdm_cg_sel = 'chid Y and name CG OD1 ND2'
lig_sel = ['C12', 'O3','N2']
add_one_helix(cg, aa, correspond_resname, represent_name, vdm_cg_sel, lig_sel)