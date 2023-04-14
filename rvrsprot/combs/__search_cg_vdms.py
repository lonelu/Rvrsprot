import prody as pr
import numpy as np
from sklearn.neighbors import NearestNeighbors
from . import gvdm_helper

class CgCombInfo:
    def __init__(self):
        self.ligand_id = None
        self.vdm_cgs = {}


def get_ligand_coords(filtered_ligands, lgd_sel):
    # Get all ligand coords.
    ligand_coords = []
    for i in range(len(filtered_ligands)):
        lgd = filtered_ligands[i]
        coords = []
        for lsa in lgd_sel:
            coords.extend(lgd.select('name ' + lsa).getCoords().flatten())
        ligand_coords.append(coords)
    ligand_coords = np.array(ligand_coords)

    return ligand_coords


def get_vdm_labels_coords(dfa, represent_name, correspond_resname, correspond_names):
    '''
    Example:
    correspond_resname = 'ASP'
    represent_name = 'OD2'
    correspond_names = ['CG', 'OD1', 'OD2']
    '''
    labels = dfa[(dfa['chain'] == 'Y') & (dfa['resname'] == correspond_resname) & (dfa['name'] == represent_name)][['CG', 'rota', 'probe_name']]

    vdm_coords =[]
    for k in correspond_names:
        df_contacts = dfa[
            (dfa['resname'] == correspond_resname) 
            & (dfa['chain'] == 'Y') 
            & (dfa['name'] == k)
        ]
        vdm_coords.append(df_contacts[['c_x', 'c_y', 'c_z']].values.T)
    vdm_coords = np.concatenate(vdm_coords).T

    return labels, vdm_coords


def get_nearest_vdms(vdm_coords, ligand_coords, radius = 2):
    '''
    #to be predecated. 
    '''
    # Nearest Neighbor
    nbr = NearestNeighbors(radius=radius).fit(vdm_coords)
    adj_matrix = nbr.radius_neighbors_graph(ligand_coords).astype(bool)
    print(adj_matrix.sum())
    print(adj_matrix.shape)

    m_adj_matrix = adj_matrix.tolil()

    all_inds = {}
    for r in range(m_adj_matrix.shape[0]):
        inds = m_adj_matrix.rows[r]
        if len(inds) <= 0:
            continue
        for ind in inds:
            all_inds.add(ind)

    return all_inds


def write_vdms(outdir, all_inds, labels, dfa, prefix):

    for ind in all_inds:
        x = labels.iloc[ind]
        # if x['seg_chain_resnum'][2] != 45:
        #     #print(x['seg_chain_resnum'][2])
        #     continue
        v = dfa[(dfa['CG'] == x['CG']) & (dfa['rota'] == x['rota']) & (dfa['probe_name'] == x['probe_name']) & (dfa['seg_chain_resnum'] == x['seg_chain_resnum'])]
        # if v['resname'].iloc[-1] != 'ARG':
        #     continue
        print(ind)
        ag = gvdm_helper.df2ag(v)
        pr.writePDB(outdir + prefix + '_' + str(ind),  ag)
        #combs2.design.functions.print_dataframe(v, outpath=outdir, tag = '_' + str(ind), prefix = prefix)
    return


def get_nearest_vdms_rmsd(vdm_coords, ligand_coords, radius = 1):
    # Nearest Neighbor
    nbr = NearestNeighbors(radius=radius).fit(vdm_coords)
    dists, inds = nbr.radius_neighbors(ligand_coords)

    return dists, inds


def vdm_ligand_clash(vdm, ligand, clash_radius = 2.7):
    '''
    The vdm X sidechain may clash with the ligand.
    return True if clash.
    '''
    vdm_coords =[]
    for i in range(vdm[(vdm['chain'] == 'X')].shape[0]):
        k = vdm[(vdm['chain'] == 'X')].iloc[i]
        if k['name'][0] == 'H':
            continue
        vdm_coords.append(k[['c_x', 'c_y', 'c_z']].values)

    ligand_coords = ligand.select('heavy').getCoords()


    if len(vdm_coords) == 0:
        print(vdm[['CG', 'rota', 'probe_name', 'chain', 'name']])
        return False
    nbr = NearestNeighbors(radius=clash_radius).fit(np.array(vdm_coords))
    adj_matrix = nbr.radius_neighbors_graph(ligand_coords).astype(bool)

    if adj_matrix.sum() > 0:
        return True
    return False


def search_vdm(dfa, ligands, cg_id, input_dict, labels_cgs, df_cgs, dist_ind_cgs, rmsd = 0.5):
    '''
    s: combs2.design._sample.Sample

    cg_id: (0, 0)  #The first is to record the cg in ligand, the second is to record the cg in vdms. 

    input_dict[cg_id] = {

        cg : 'coo'  #dfa = s.cg_dict['coo']

        lgd_sel : ['C9', 'O3', 'O4']  # correpond to coo: ASP ['CG', 'OD1', 'OD2'] 

        represent_name : 'OD2'

        correspond_resname : 'ASP'

        correspond_names : ['CG', 'OD1', 'OD2']
    }
    '''

    #dfa = cg_dict[input_dict[cg_id]['cg']]

    ligand_coords = get_ligand_coords(ligands, input_dict[cg_id]['lgd_sel'])

    labels, vdm_coords = get_vdm_labels_coords(dfa, input_dict[cg_id]['represent_name'], input_dict[cg_id]['correspond_resname'], input_dict[cg_id]['correspond_names'])

    if vdm_coords.shape[0] <=0:
        print('No vdM coords are generated.')
        return
    
    num_cg_atoms = len(input_dict[cg_id]['lgd_sel'])
    radius = np.sqrt(num_cg_atoms) * rmsd

    dists, inds = get_nearest_vdms_rmsd(vdm_coords, ligand_coords, radius = radius)

    labels_cgs[cg_id] = labels
    df_cgs[cg_id] = dfa
    dist_ind_cgs[cg_id] = (dists, inds)
    print('{} Found {} vdm'.format(cg_id, sum([len(x) for x in inds])))
    return 


def construct_vdm_write(outdir, ligands, labels_cgs, input_dict, df_cgs, dist_ind_cgs, clash_radius = 2.7, benchmark_filters = []):
    '''
    dist_ind_cgs: dict. {cg: (dists, inds)}, where dists in shape: (len(ligands), )
    df_cgs: {cg: df}
    '''
    CgCombInfoDict = {}

    for i in range(len(ligands)):
        cgCombInfo = CgCombInfo()
        cgCombInfo.ligand_id = i

        for cg in dist_ind_cgs.keys():
            info = []

            num_cg_atoms = len(input_dict[cg]['lgd_sel'])
            labels = labels_cgs[cg]
            dfa = df_cgs[cg]

            dists = dist_ind_cgs[cg][0][i]
            inds = dist_ind_cgs[cg][1][i]

            for j in range(len(inds)):
                ind = inds[j]
                rmsd = dists[j]/np.sqrt(num_cg_atoms)
                if ind > labels.shape[0]:
                    print(cg)
                    print(labels.shape)
                    print(inds)
                x = labels.iloc[ind]
                v = dfa[(dfa['CG'] == x['CG']) & (dfa['rota'] == x['rota']) & (dfa['probe_name'] == x['probe_name'])]
                if vdm_ligand_clash(v, ligands[i], clash_radius):
                    print('clash')
                    continue

                info.append((round(rmsd, 2), round(v['C_score_bb_ind'].iloc[0], 2), v)) 
                prefix = 'Lig-'+ str(i) + '_key_' + '-'.join(str(z) for z in cg) + '_rmsd_' + str(round(rmsd, 2))  + '_' + str(round(v['C_score_bb_ind'].iloc[0], 2)) + '_'                              
                #combs2.design.functions.print_dataframe(v, outpath=outdir, tag = '_' + str(ind), prefix = prefix)
                ag = gvdm_helper.df2ag(v)
                pr.writePDB(outdir + prefix + '_' + str(ind),  ag)
            cgCombInfo.vdm_cgs[cg] = info
        CgCombInfoDict[i] = cgCombInfo
        
        if any([len(cgCombInfo.vdm_cgs[k]) > 0 for k in cgCombInfo.vdm_cgs.keys()]):
            pr.writePDB(outdir + 'Lig-' + str(i) + '_' + ligands[i].getTitle(), ligands[i])

    return CgCombInfoDict

            
def write_summary(outdir, CgCombInfoDict, name = '_summary.tsv'):
    with open(outdir + name, 'w') as f:
        f.write('ligand_id\tcg_id\trmsd\n')
        for k in CgCombInfoDict.keys():          
            cgInfo = CgCombInfoDict[k]
            for cg in cgInfo.vdm_cgs.keys():
                info = cgInfo.vdm_cgs[cg]
                for o in info:
                    f.write(str(cgInfo.ligand_id) + '\t')
                    f.write(str(cg) + '\t')
                    f.write(str(o[0]) + '\t')
                    f.write(str(o[1]) +'\n')
    return