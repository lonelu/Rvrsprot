import numpy as np
import prody as pr
import pandas as pd


def rvdm_transform(ligand, lgd_sel_atomnames, cg_sel_atomnames, df_vdm):
    lig_coord = np.array([ligand.select('name ' + x).getCoords()[0] for x in lgd_sel_atomnames])
    df_vdm_first = df_vdm[['CG', 'rota', 'probe_name']].drop_duplicates().iloc[0]

    First_cg =  df_vdm[(df_vdm['CG'] == df_vdm_first['CG']) 
                        & (df_vdm['rota'] == df_vdm_first['rota']) 
                        & (df_vdm['probe_name'] == df_vdm_first['probe_name'])]
    
    First_cg_sel_coords = np.array([First_cg[(df_vdm['chain'] == 'Y') & (df_vdm['name']== x)][['c_x', 'c_y', 'c_z']].values[0] for x in cg_sel_atomnames])
    
    tf = pr.calcTransformation(First_cg_sel_coords, lig_coord)
    df_vdm.loc[:, ['c_x', 'c_y', 'c_z']] = tf.apply(df_vdm[['c_x', 'c_y', 'c_z']].to_numpy())
    return 