import prody as pr
import pickle
import pandas as pd
import numpy as np


def load_new_vdm(path_to_database, cg, aa):
    '''
    The function is used to load vdM database after 2022.02.
    Pelase check combs2_da.ipynb to learn the loaded database.
    '''
    cg_aa = cg + '/' + aa 
    df_vdm = pd.read_parquet(path_to_database + 'vdMs/' + cg_aa + + '.parquet.gzip')

    with open(path_to_database + 'vdMs_gr_indices/' + cg_aa + '.pkl', 'rb') as f:
        df_gr = pickle.load(f)

    df_score = pd.read_parquet(path_to_database + 'nbrs/vdMs_cg_nbrs_scores/' + cg_aa + '.parquet.gzip')

    return df_vdm, df_gr, df_score


def load_old_vdm(path_to_database, cg, aa):
    '''
    The function is used to load vdM database before 2022.02.
    Pelase check combs2_da.ipynb to learn the loaded database.
    '''
    cg_aa = cg + '/' + aa 
    df_vdm = pd.read_parquet(path_to_database + cg_aa + '.parquet.gzip')

    return df_vdm


def ligscoords_2_ideal_ala(pos, ligs, ideal_ala_coords):
    '''
    For all the pre generated ligs, combine them to each position in bb and then transform the pos_lig to 'ideal alanine'.
    The 'ideal alanine' is the alanine used for the vdM database.
    TO DO: this function can be more efficient to calc the matrix instead of calc each ligs.
    '''
    pos_coords = pos.select('name N CA C').getCoords()
    tf = pr.calcTransformation(pos_coords, ideal_ala_coords)
    _ligs = []
    for lig in ligs:
        _lig = lig.copy()
        tf.apply(_lig)
        _ligs.append(_lig)
    
    tf_rev = pr.calcTransformation(ideal_ala_coords, pos_coords)
    return _ligs, tf_rev, tf


def get_vdm_labels_coords_4old_vdm_db(dfa, correspond_resname, represent_name, correspond_names):
    '''
    Example:
    correspond_resname = 'ASP'
    represent_name = 'OD2'
    correspond_names = ['CG', 'OD1', 'OD2']
    '''
    labels = dfa[(dfa['chain'] == 'Y') & (dfa['resname'] == correspond_resname) & (dfa['name'] == represent_name)][['CG', 'rota', 'probe_name']]
    #labels = dfa[['CG', 'rota', 'probe_name']].drop_duplicates()

    vdm_coords =[]
    for k in correspond_names:
        df_contacts = dfa[
            (dfa['resname'] == correspond_resname) 
            & (dfa['chain'] == 'Y') 
            & (dfa['name'] == k)
        ]
        vdm_coords.append(df_contacts[['c_x', 'c_y', 'c_z']].values.T)
    try:
        vdm_coords = np.concatenate(vdm_coords).T
    except:
        vdm_coords = np.array([])
        print(dfa.head(50))
        print(correspond_resname)
        print(correspond_names)
        print('Error: get_vdm_labels_coords_4old_vdm_db. labels length is ' + str(len(labels)))

    return labels, vdm_coords


def filter_db(df_vdm, use_enriched = True, use_abple = True, abple = 'A'):
    '''
    Filter database based on ABPLE, enriched. Check combs2._sample._load_res().
    '''
    if use_enriched and use_abple:
        return df_vdm[df_vdm['C_score_' + 'ABPLE_' + abple] > 0]
    if use_enriched and not use_abple:
        return df_vdm[df_vdm['C_score_bb_ind'] > 0]
    return df_vdm


def df2ag(df, title = None, b_factor_column=None):
    '''
    The df here only contains one vdM.
    '''
    df = df.copy()
    if 'chain' not in df.columns:
        df['chain'] = 'A'
    ag = pr.AtomGroup()
    ag.setCoords(df[['c_x','c_y','c_z']].values)
    ag.setResnums(df['resnum'].values)
    ag.setResnames(df['resname'].values)
    ag.setNames(df['name'].values)
    ag.setChids(df['chain'].values)
    ag.setSegnames(df['chain'].values)
    heteroflags = ag.getSegnames() == 'L'
    ag.setFlags('hetatm', heteroflags)
    if title:
        ag.setTitle(title)
    if 'beta' in df.columns and b_factor_column is None:
        ag.setBetas(df['beta'].values)
    elif b_factor_column is not None:
        ag.setBetas(df[b_factor_column].values)
    if 'occ' not in df.columns:   
        df['occ'] = 1
    #pr.writePDB(outpath + prefix + filename + tag + '.pdb.gz', ag, occupancy=df['occ'].values)
    return ag


