import prody as pr
import numpy as np
from sklearn.neighbors import NearestNeighbors

from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.rdmolfiles import MolToPDBFile

from . import fragment_handle 


def getSele(residues):
    sele = 'name N CA C O or resnum '
    xs = ' '.join([str(r) for r in residues])
    sele = sele + xs
    return sele

def calcPDBcoords(prot, residues):
    sele = getSele(residues)
    #print(sele)
    coords = prot.select(sele).getCoords()
    return coords


def run_one_mol(mol, m_start, prot_coords, radius):
    ccc = fragment_handle.Get_mol_coords(mol, m_start)
    if ccc is None:
        return None
    mol_coords, conf_num, heavy_atom_num = ccc
    nbr = NearestNeighbors(radius=radius).fit(prot_coords)
    adj_matrix = nbr.radius_neighbors_graph(mol_coords).astype(bool)

    return adj_matrix, conf_num, heavy_atom_num


def extract_excluded(adj_matrix, conf_num, heavy_atom_num):

    m_adj_matrix = adj_matrix.tolil()

    all_inds = np.zeros(conf_num * heavy_atom_num)

    for r in range(conf_num * heavy_atom_num):
        inds = m_adj_matrix.rows[r]
        if len(inds) <= 0:
            continue
        all_inds[r] = 1

    #print(sum(all_inds))

    all_inds = all_inds.reshape((conf_num, heavy_atom_num))

    results = np.sum(all_inds, axis = 1)

    excluded = results >= 1
    
    return excluded


def get_all_items(workdir, pdbname, m_start_name = 'CMR.sdf', mol_name = 'molport_2.sdf'):
    
    #workdir = '/Users/lonelu/DesignData/fragment_design/FragmentScreen_yuda/'
    #pdbname = '7OCoumarin fragment ABLE.pdb'
    prot = pr.parsePDB(workdir + pdbname)

    suppl = Chem.SDMolSupplier(workdir + m_start_name)
    m_start = [x for x in suppl][0]

    molpath = workdir +  mol_name

    suppl = Chem.SDMolSupplier(molpath)

    mols = [x for x in suppl]

    return prot, m_start, mols


def write_summary(outdir, infos, key = 'Molport ID'):
    with open(outdir + '_summary.txt', 'w') as f:
        f.write('mol_id\tcid\n')
        for info in infos:
            mol, cids, excluded = info
            if sum(1-excluded) <=0:
                continue
            f.write(mol.GetProp(key) + '\t' + str(sum(1-excluded)) + '\n')





