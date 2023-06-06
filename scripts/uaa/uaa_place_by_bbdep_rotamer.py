'''
Use the bbdep rotamer library to place aa on selected positions.
Then align the shared cg of uaa to the aa. 
Remove clashes.
'''

import os
import prody as pr
import pandas as pd
import numpy as np

from rdkit import Chem
import rdkit.Chem.rdMolTransforms as rmt
from rdkit.Chem.rdmolfiles import MolToPDBFile 
from itertools import combinations

from sklearn.neighbors import NearestNeighbors
from rvrsprot.basic import convex_hull


bbdep_database_dir = '/mnt/e/GitHub_Design/Rvrsprot/scripts/uaa/MEDFORD-rotamer-library-main/'

def load_bbdep_database(bbdep_database_dir, AA):
    _db = bbdep_database_dir + AA + '.bbdep.rotamers.lib'
    db = pd.read_table(_db, skiprows=1, header=None, delim_whitespace=True)
    with open(_db, 'r') as f:
        title = f.readline()
    names = title[1:].split()
    db.columns = names

    phipsi_dict = {}
    for x in range(-18, 18):
        for y in range(-5, 5):
            z = x*10 + y
            if z < -180:
                z = 360+z
            phipsi_dict[z]=x*10
    return db, phipsi_dict


def extract_phipsi(target, chid, resnum):
    res = target.getHierView().getResidue(chid, resnum)
    phi = pr.calcPhi(res)
    psi = pr.calcPsi(res)
    return phi, psi



def extract_bbdepd(phi, psi, bbdep_rotamer_db, phipsi_dict):
    _phi = phipsi_dict[int(phi)]
    _psi = phipsi_dict[int(psi)]
    sel_db = bbdep_rotamer_db[(bbdep_rotamer_db['Phi'] == _phi) & (bbdep_rotamer_db['Psi'] == _psi)]
    return sel_db


def generate_rotamers(uaa_rd, chi1, chi2, angle_step = np.linspace(-16,16,8)):
    '''
    
    '''
    uaas = []
    for ags1, ags2 in combinations(angle_step, 2):
        _chi1 = chi1 + ags1
        _chi2 = chi2 + ags2
        uaa_cp = Chem.Mol(uaa_rd)
        ucc = uaa_cp.GetConformer()
        rmt.SetDihedralDeg(ucc, 0, 1, 4, 5, _chi1)
        rmt.SetDihedralDeg(ucc, 1, 4, 5, 6, _chi2)
        uaas.append(uaa_cp)
        #Rotate 180
        uaa_cp_cp = Chem.Mol(uaa_cp)
        ucc2 = uaa_cp_cp.GetConformer()
        x = rmt.GetDihedralDeg(ucc, 9, 10, 11, 12)
        rmt.SetDihedralDeg(ucc2, 9, 10, 11, 12, x+180)
        uaas.append(uaa_cp_cp)
    return uaas


def align_uaa_target_by_pose(outdir, uaa, target,  resnum, chid= 'A'):
    sel = 'chid ' + chid + ' resnum ' + str(resnum) + ' name N CA C'
    pr.calcTransformation(uaa.select('name N CA C'), target.select(sel)).apply(uaa)
    pr.writePDB(outdir + uaa.getTitle() + '_' + chid + '_' + str(resnum), uaa )
    return


def extract_uaard_coords(uaas, ignore_atom_ids = [0, 1, 2, 3, 4, 19]):
    '''
    
    '''
    coords = []
    for uaa in uaas:
        for i, atom in enumerate(uaa.GetAtoms()):
            if i in ignore_atom_ids:
                continue
            positions = uaa.GetConformer().GetAtomPosition(i)
            coords.append([positions.x, positions.y, positions.z])
    return np.array(coords), (len(uaas), len(uaa.GetAtoms())-len(ignore_atom_ids))


def extract_single_uaard_coords(uaa, ignore_atom_ids = [0, 1, 2, 3, 4, 19]):
    '''
    
    '''
    coords = []

    for i, atom in enumerate(uaa.GetAtoms()):
        if i in ignore_atom_ids:
            continue
        positions = uaa.GetConformer().GetAtomPosition(i)
        coords.append([positions.x, positions.y, positions.z])
    return np.array(coords)


def get_clashing_uaa_id(target_coords, uaa_coords, uaa_num, uaa_atom_num, rmsd = 1.5):
    '''
    #to be predecated. 
    '''
    # Nearest Neighbor
    nbr = NearestNeighbors(radius=rmsd).fit(uaa_coords)
    adj_matrix = nbr.radius_neighbors_graph(target_coords).astype(bool)
    #print(adj_matrix.sum())
    #print(adj_matrix.shape)

    m_adj_matrix = adj_matrix.tolil()

    #get clashed inds.
    all_inds = set()
    for r in range(m_adj_matrix.shape[0]):
        inds = m_adj_matrix.rows[r]
        if len(inds) <= 0:
            continue
        for ind in inds:
            all_inds.add(int(ind/uaa_atom_num))

    return all_inds


def load_uaa_by_pos_filter_by_clash(dir_uaa_res, outdir, target, chid, resnum, rmsd):
    '''
    
    '''
    #target_coords = target.select('protein and name N CA C O').getCoords()
    target_coords = target.select('( protein and name N CA C O ) or resname RPW').getCoords()
    phi, psi = extract_phipsi(target, chid, resnum)
    sel_db = extract_bbdepd(phi, psi, bbdep_rotamer_db, phipsi_dict)

    file_name = 'AzoPhe_' + chid + '_' + str(resnum) + '.pdb'
    uaa_rd = Chem.MolFromPDBFile(dir_uaa_res + file_name)

    for ikk in range(sel_db.shape[0]):
        chi1 = sel_db.iloc[ikk]['chi1Val']
        chi2 = sel_db.iloc[ikk]['chi2Val']
        prob = sel_db.iloc[ikk]['Probabil']

        uaas = generate_rotamers(uaa_rd, chi1, chi2)
        #[MolToPDBFile(uaas[i], workdir + '_test_dir/' + str(i) + 'test.pdb') for i in range(len(uaas))]
        uaa_coords, (uaa_num, uaa_atom_num) = extract_uaard_coords(uaas)
        all_inds = get_clashing_uaa_id(target_coords, uaa_coords, uaa_num, uaa_atom_num, rmsd)
        potential_ids = [_ix for _ix in range(len(uaas)) if _ix not in all_inds]

        
        candidate_ids = []
        for pid in potential_ids:
            uaa_coord = extract_single_uaard_coords(uaas[pid])
            if not filterout_by_ahull(target, uaa_coord):
                candidate_ids.append(pid)
        print((len(potential_ids), len(candidate_ids)))
        out_file = outdir + chid + '_'+ str(resnum) + '_' + str(int(phi)) + '_' + str(int(psi)) 
        out_file = out_file + '_' + str(ikk) + '_' + str(sel_db.shape[0]) + '_' + str(round(prob, 2))
        [MolToPDBFile(uaas[izz], out_file + '_' + str(izz) + '.pdb') for izz in candidate_ids]

    return


def filterout_by_ahull(target, uaa_coord):
    ahull = convex_hull.AlphaHull(alpha=9) 
    ahull.set_coords(target)
    ahull.calc_hull()
    in_hull = ahull.pnts_in_hull(uaa_coord)
    # dist_to_hull = ahull.get_pnts_distance(uaa_coord) 
    if sum(in_hull) > len(uaa_coord)/3:
        return False
    return True

###-------------------------------------------------------------------------

bbdep_rotamer_db, phipsi_dict = load_bbdep_database(bbdep_database_dir, 'phe')


workdir = '/mnt/e/DesignData/Chemodrugs/HB_RUC_uaa/'
target = pr.parsePDB(workdir + 'PiB.pdb')
uaa = pr.parsePDB(workdir + 'AzoPhe.pdb')

dir_uaa_res = workdir + 'uaa_res/'
os.makedirs(dir_uaa_res, exist_ok=True)
for resnum in target.select('name CA').getResnums():
    align_uaa_target_by_pose(dir_uaa_res, uaa, target, resnum, chid= 'A')


#phi, psi = extract_phipsi(target, 'A', 15)
#sel_db = extract_bbdepd(phi, psi, bbdep_rotamer_db, phipsi_dict)
#chi1 = sel_db.iloc[0]['chi1Val']
#chi2 = sel_db.iloc[0]['chi2Val']


workdir = '/mnt/e/DesignData/Chemodrugs/HB_RUC_uaa/'
#uaa_rd = Chem.MolFromPDBFile(workdir + 'AzoPhe.pdb')
#uc = uaa_rd.GetConformer()

#rmt.GetAngleDeg(uc, 1,2, 3)
#rmt.SetAngleDeg(uc, 1, 2, 3, 60)
#MolToPDBFile(uaa_rd, workdir + 'test.pdb')

#rmt.GetDihedralDeg(uc, 0, 1, 4, 5)
#rmt.GetDihedralDeg(uc, 1, 4, 5, 6)

#rmt.SetDihedralDeg(uc, 1, 4, 5, 6, 60)
#MolToPDBFile(uaa_rd, workdir + 'test2.pdb')

outdir = workdir + 'outdir2/'
os.makedirs(outdir, exist_ok=True)

for resnum in target.select('protein and name CA').getResnums()[1:146]:
    print(resnum)
    load_uaa_by_pos_filter_by_clash(dir_uaa_res, outdir, target, 'A', resnum, 2.3)

