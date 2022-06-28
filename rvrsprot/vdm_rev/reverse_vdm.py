'''
In Combs2 database, the vdMs are all superimposed on aa backbone. 
We assume that reverse the vdMs on cg could be useful for function based design. 
The functions here are trying to prepare the rev vdms based on Combs2 database.
'''

from metalprot.combs import search_lig_indep, gvdm_helper
from metalprot.basic import transformation
from metalprot.basic import prody_ext
import numpy as np
import prody as pr
import pandas as pd

def generate_rev_vdm_lib(path_to_database, outdir, cg, aa):
    df_vdm = gvdm_helper.load_old_vdm(path_to_database, cg, aa)

    df_vdm_rep = df_vdm[['CG', 'rota', 'probe_name']].drop_duplicates()

    First_cg_coords_y =  df_vdm[(df_vdm['chain'] == 'Y') 
                        & (df_vdm['CG'] == df_vdm_rep.iloc[0]['CG']) 
                        & (df_vdm['rota'] == df_vdm_rep.iloc[0]['rota']) 
                        & (df_vdm['probe_name'] == df_vdm_rep.iloc[0]['probe_name'])][['c_x', 'c_y', 'c_z']].values

    for i in range(1, df_vdm_rep.shape[0]):
        i_cg_coords_y =  df_vdm[(df_vdm['chain'] == 'Y') 
                        & (df_vdm['CG'] == df_vdm_rep.iloc[i]['CG']) 
                        & (df_vdm['rota'] == df_vdm_rep.iloc[i]['rota']) 
                        & (df_vdm['probe_name'] == df_vdm_rep.iloc[i]['probe_name'])][['c_x', 'c_y', 'c_z']].values
        
        i_cg =  df_vdm[(df_vdm['CG'] == df_vdm_rep.iloc[i]['CG']) 
                        & (df_vdm['rota'] == df_vdm_rep.iloc[i]['rota']) 
                        & (df_vdm['probe_name'] == df_vdm_rep.iloc[i]['probe_name'])]
        i_cg_coords_all = i_cg[['c_x', 'c_y', 'c_z']].values

        R, m_com, t_com = transformation.get_rot_trans(i_cg_coords_y, First_cg_coords_y)
        
        # The following is not workding. https://www.dataquest.io/blog/settingwithcopywarning/
        #i_cg.loc[i_cg.index][['c_x', 'c_y', 'c_z']] = np.dot(i_cg_coords_all - m_com, R) + t_com
        df_vdm.loc[i_cg.index, ['c_x', 'c_y', 'c_z']] = np.dot(i_cg_coords_all - m_com, R) + t_com

    df_vdm.to_parquet(outdir + aa + '.parquet.gzip', compression = 'gzip')
    return


def write_vdm(path_to_database, cg, aa, outdir, total = None):
    df_vdm = gvdm_helper.load_old_vdm(path_to_database, cg, aa)
    df_vdm_rep = df_vdm[['CG', 'rota', 'probe_name']].drop_duplicates()

    if total:
        if total > df_vdm.shape[0]:
            total = df_vdm.shape[0]
    else:
        total = df_vdm.shape[0]

    for i in range(0, total):
        i_cg =  df_vdm[(df_vdm['CG'] == df_vdm_rep.iloc[i]['CG']) 
                        & (df_vdm['rota'] == df_vdm_rep.iloc[i]['rota']) 
                        & (df_vdm['probe_name'] == df_vdm_rep.iloc[i]['probe_name'])]
        agi = gvdm_helper.df2ag(i_cg)
        pr.writePDB(outdir + cg + '_' + aa  + '_' + str(i) + '.pdb', agi)     
    return


def vdm_add_helix(df_vdm, CG, rota, probe_name, helix, helix_chidres, helix_cg = None, check_clash = True):
    '''
    Attach vdM in a std helix.
    '''
    vdm_pd = df_vdm[(df_vdm['CG'] == CG) 
                    & (df_vdm['rota'] == rota) 
                    & (df_vdm['probe_name'] == probe_name)]
    vdm = gvdm_helper.df2ag(vdm_pd)
    if not vdm:
        print('Wrong sel: h_' + str(CG) + '_' + str(rota) + '_' + str(probe_name))
        return
    #print(vdm.select('chid X and resnum 10 and name N CA C'))
    _helix = helix.copy()
    tf = pr.calcTransformation(_helix.select('chid ' + helix_chidres[0] + ' and resnum ' + str(helix_chidres[1]) + ' and name N CA C'), vdm.select('chid X and resnum 10 and name N CA C'))
    tf.apply(_helix)
    helix_resind = _helix.select('chid ' + helix_chidres[0] + ' and resnum ' + str(helix_chidres[1])).getResindices()[0]
    ind_vdm_dict = {helix_resind: vdm}
    title = str(CG) + '_' + str(rota) + '_' + str(probe_name) + '_sc_' + str(round(vdm_pd['C_score_ABPLE_A'].iloc[0], 1))
    #print(ind_vdm_dict)
    if check_clash and search_lig_indep.clash_filter_protein_single(_helix, (helix_chidres[0], helix_chidres[1]), vdm):
        #print('Clash '+ title)
        return
    ag = prody_ext.mutate_vdm_target_into_ag2(_helix, ind_vdm_dict, title, vdm_sel = 'chid X and resnum 10')
    if helix_cg:
        if search_lig_indep.clash_filter_proteins([helix_cg, ag], 2.5):
            #print('HH Clash '+ title)
            return
        ag = prody_ext.combine_ags([helix_cg, ag], ag.getTitle())
    return ag

def generate_vdmH(outdir, path_to_database, cg, aa, helix_path, helix_resind):
    
    df_vdm = gvdm_helper.load_old_vdm(path_to_database, cg, aa)
    df_vdm_f = gvdm_helper.filter_db(df_vdm, use_enriched = True, use_abple = True, abple = 'A')
    df_vdm_rep = df_vdm_f[['CG', 'rota', 'probe_name']].drop_duplicates()
    helix = pr.parsePDB(helix_path)    

    for i in range(0, df_vdm_rep.shape[0]):
        _cg = df_vdm_rep.iloc[i]['CG']
        _rota = df_vdm_rep.iloc[i]['rota']
        _probe_name = df_vdm_rep.iloc[i]['probe_name']
        ag = vdm_add_helix(df_vdm_f, _cg, _rota, _probe_name, helix, helix_resind)
        if ag:
            pr.writePDB(outdir + ag.getTitle(), ag)
    return 


def vdm_cg_add_helix(vdm, vdm_cg_sel, helix_cg, helix_cg_sel):
    '''
    The helix here has an aa which will be used to superimpose the cg of vdm.
    '''
    _helix_cg = helix_cg.copy()
    tf = pr.calcTransformation(_helix_cg.select(helix_cg_sel), vdm.select(vdm_cg_sel))
    tf.apply(_helix_cg)
    return _helix_cg

def generate_vdm2H(path_to_database, cg, aa, helix_path, helix_resind, helix_cg_path, helix_cg_sel, vdm_cg_sel):
    ags = []
    df_vdm = gvdm_helper.load_old_vdm(path_to_database, cg, aa)
    df_vdm_f = gvdm_helper.filter_db(df_vdm, use_enriched = True, use_abple = True, abple = 'A')
    df_vdm_rep = df_vdm_f[['CG', 'rota', 'probe_name']].drop_duplicates()
    helix = pr.parsePDB(helix_path)    

    helix_cg = pr.parsePDB(helix_cg_path)  
    vdm_pd = df_vdm[(df_vdm['CG'] == df_vdm_rep.iloc[0]['CG']) 
                    & (df_vdm['rota'] == df_vdm_rep.iloc[0]['rota']) 
                    & (df_vdm['probe_name'] == df_vdm_rep.iloc[0]['probe_name'])]
    vdm0 = gvdm_helper.df2ag(vdm_pd)
    _helix_cg = vdm_cg_add_helix(vdm0, vdm_cg_sel, helix_cg, helix_cg_sel)
    
    for i in range(0, df_vdm_rep.shape[0]):
        _cg = df_vdm_rep.iloc[i]['CG']
        _rota = df_vdm_rep.iloc[i]['rota']
        _probe_name = df_vdm_rep.iloc[i]['probe_name']
        ag = vdm_add_helix(df_vdm_f, _cg, _rota, _probe_name, helix, helix_resind, _helix_cg)
        if ag:
            ags.append(ag)
    return ags

