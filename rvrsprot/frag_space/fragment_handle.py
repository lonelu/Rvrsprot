from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.rdmolfiles import MolToPDBFile

import prody as pr
import numpy as np

def gen_atomMap(mol_mobile, m_start):
    #need to be used before addH and embed.
    if not mol_mobile.HasSubstructMatch(m_start):
        return None
    match_atom_ids = mol_mobile.GetSubstructMatches(m_start)
    atomMap = [(match_atom_ids[0][i], i) for i in range(len(match_atom_ids[0]))]
    return atomMap

# def gen_atomMaps(mols, m_start):
#     atomMaps = []
#     for mol in mols:
#         atomMap = gen_atomMap(mol, m_start)
#     atomMaps.append(atomMap)


def embedMolecule(mol, numConfs = 1):
    Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=numConfs, randomSeed=0xf00d)
    return list(cids)

def embedMolecule2(mol, numConfs=100, maxAttempts=1000, pruneRmsThresh=0.05, 
                   useExpTorsionAnglePrefs=True, useBasicKnowledge=True, enforceChirality=True):
    Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=numConfs, maxAttempts=maxAttempts, pruneRmsThresh=pruneRmsThresh, useExpTorsionAnglePrefs=useExpTorsionAnglePrefs, useBasicKnowledge=useBasicKnowledge, enforceChirality=enforceChirality, numThreads=0)
    return list(cids)


# def embedAllMolecules(mols, numConfs = 100, maxAttempts =1000, pruneRmsThresh = 0.2):
#     #In the future, implement multi-thread.
#     '''
#     We may do not want the numConfs to be too large. 
#     From a practical consideration, we may only want to design ligands with less flexibility so far.
#     '''
#     cidsLists = []
#     for mol in mols:
#         cids = embedMolecule2(mol, numConfs,maxAttempts, pruneRmsThresh)
#         cidsLists.append(cids)
#     return cids


def mol_super_mol(mol_mobile, m_start, _atomMap, cids):
    for id in cids:
        Chem.rdMolAlign.AlignMol(mol_mobile, m_start, prbCid = id, atomMap = _atomMap)
    return 


def calc_mol_coords(mol_mobile, cids):
    coords = []
    heavy_atom_num = Chem.rdMolDescriptors.CalcNumHeavyAtoms(mol_mobile)
    for id in cids:
        conf = mol_mobile.GetConformer(id)
        coord = conf.GetPositions()[:heavy_atom_num]
        coords.append(coord)
    return np.array(coords).reshape(-1, 3), heavy_atom_num


def Get_mol_coords(mol_mobile, m_start):
    _atomMap = gen_atomMap(mol_mobile, m_start)
    if _atomMap is None:
        return None
    cids = embedMolecule2(mol_mobile)
    mol_super_mol(mol_mobile, m_start, _atomMap, cids)
    coords, heavy_atom_num = calc_mol_coords(mol_mobile, cids)

    return coords, cids, heavy_atom_num


def removeConformers(mol, cids, excluded):
    for cid in cids:
        if excluded[cid]:
            mol.RemoveConformer(cid)
    return


def writeSDF(outdir, mol, cids, excluded, key = 'Molport ID'):

    title = mol.GetProp(key) 
    if sum(1-excluded) <=0:
        print('no conformer for ' + title)
        return

    writer = Chem.SDWriter(outdir + title + '.sdf')
    for cid in cids:
        if excluded[cid]:
            continue
        mol.SetProp('ID', f'_conf_{cid}')
        writer.write(mol, confId=cid)
    return

def writeMolToPDB(outdir, mol, cids, excluded, key ='Molport ID'):
    for cid in cids:
        if excluded[cid]:
            continue
        MolToPDBFile(mol, outdir + mol.GetProp(key) + '_conf' + str(cid) + '.pdb', confId = cid)


