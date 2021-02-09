import os
import sys
import string

import numpy as np
import prody as pr
import numba

from scipy.linalg import toeplitz
from scipy.spatial import Delaunay
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser, PDBIO, Select
from qbits import convex_hull, pdb, clash
from itertools import combinations

parser = PDBParser(QUIET=True)
io = PDBIO()

class NotDisorderedOrH(Select):
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == "A" or \
               not "H" in atom.get_name()

def get_struct(name, pdb_path, min_nbrs=0):
    """Using the Biopython parser, get structure from a PDB file.

    Parameters
    ----------
    name : str
        Name for the output structure object.
    pdb_path : str
        Path to PDB file containing the structure.
    min_nbrs : int, optional
        Minimum number of neighbors of residues to be included in the 
        output structure. The b-factor on the alpha-carbon of a residue 
        is assumed to the the log of the number of neighbors.

    Returns
    -------
    struct : Bio.PDB.Structure
        Structure extracted from PDB file.
    """
    struct = parser.get_structure(name, pdb_path)
    # remove residues with fewer than the minimum number of neighbors
    if min_nbrs:
        min_b = np.log(min_nbrs)
    for chain in struct.get_chains():
        ids_to_detach = []
        for res in chain.get_residues():
            if min_nbrs:
                bfac = res['CA'].get_bfactor()
                if bfac is not None and bfac < min_b:
                    ids_to_detach.append(res.id)
            # set segid to the empty string
            res.segid = ''
        int_ids = [res_id[1] for res_id in ids_to_detach]
        for i, res_id in enumerate(ids_to_detach):
            # ensure residues are not removed from the middle of the struct
            if int_ids[:i] == list(range(min(int_ids), int_ids[i])) or \
                    int_ids[i:] == list(range(int_ids[i], max(int_ids) + 1)):
                chain.detach_child(res_id)
    return struct


def merge_structs(structs, slices=[]):
    """Merge a list of BioPython structs into one.

    Parameters
    ----------
    structs : list
        List of BioPython structs.
    slices : list, optional
        List of str to be evaluated as Python-style list slices of residues 
        from each PDB structure to be included in the merged structure.  
        Useful for creating MASTER loop queries.
    """
    idx = 0
    # Set chains in structures and move to first structure
    for x, structure in enumerate(structs):
        if len(slices) == len(structs):
            res_ids = [res.id[1] for res in structure.get_residues()]
        for chain in structure.get_chains():
            chain.id = string.ascii_uppercase[idx]
            idx += 1
            # detach residues that do not meet the neighbor threshold
            if len(slices) == len(structs):
                ids_to_keep = eval('res_ids[{}]'.format(slices[x]))
                ids_to_detach = [res_id for res_id in res_ids if 
                                 res_id not in ids_to_keep]
            else:
                ids_to_detach = []
            # ensure there are no floating residues
            diffs = [ids_to_detach[i + 1] - ids_to_detach[i] for i in 
                     range(len(ids_to_detach) - 1)]
            for res_id, diff in zip(ids_to_detach, diffs):
                if diff != max(diffs) and diff > 1:
                    for j in range(1, diff):
                        ids_to_detach.append(res_id + j)
            # remove residues in the region beyond the high-neighbor core
            for res_id in set(ids_to_detach):        
                chain.detach_child((' ', res_id, ' '))
            # Don't move chains of struct[0]
            if x == 0:
                continue
            chain.detach_parent()
            structs[0][0].add(chain)
    return structs[0]

def merge_save_struct(out_path, structs, slices=[]):
    idx = 0
    # Set chains in structures and move to first structure
    
    for x, structure in enumerate(structs):
        if len(slices) == len(structs):
            res_ids = [res.id[1] for res in structure.get_residues()]
        chain = list(structure.get_chains())[0]
        chain.id = string.ascii_uppercase[idx]
        idx += 1
        # detach residues that do not meet the neighbor threshold
        if len(slices) == len(structs):
            ids_to_keep = eval('res_ids[{}]'.format(slices[x]))
            ids_to_detach = [res_id for res_id in res_ids if 
                                res_id not in ids_to_keep]
        else:
            ids_to_detach = []
        # ensure there are no floating residues
        diffs = [ids_to_detach[i + 1] - ids_to_detach[i] for i in 
                    range(len(ids_to_detach) - 1)]
        for res_id, diff in zip(ids_to_detach, diffs):
            if diff != max(diffs) and diff > 1:
                for j in range(1, diff):
                    ids_to_detach.append(res_id + j)
        # remove residues in the region beyond the high-neighbor core
        for res_id in set(ids_to_detach):        
            chain.detach_child((' ', res_id, ' '))
        # Don't move chains of struct[0]
        if x == 0:
            chain0 = list(structs[0].get_chains())[0]
            continue
        chain.detach_parent()
        init_id = max([res.id[1] for res in chain0]) + 1
        for res_id, res in enumerate(chain.get_residues()):
            res.detach_parent()
            res.id = (' ', init_id + res_id, ' ')
            chain0.add(res)
    io.set_structure(chain0)
    io.save(out_path, select=NotDisorderedOrH())

def merge_pdbs(pdb_paths, out_path, min_nbrs=0, set_bfac=None):
    """Merge a list of PDB files into one countaining all structure from each.

    Parameters
    ----------
    pdb_paths : list
        List of paths to the PDB files for the structures.
    out_path : str
        Path to the output PDB file.
    min_nbrs : int, optional
        Minimum number of neighbors of residues to be included in the 
        output structure. The b-factor on the alpha-carbon of a residue 
        is assumed to the the log of the number of neighbors.
    set_bfac : float, optional
        Float value to which to set b-factors in the output PDB.
    """
    structs = []
    counter = 0
    for i in range(len(pdb_paths)):
        title = 'struct' + str(i)
        structs.append(get_struct(title, pdb_paths[i], min_nbrs))
        # relabel chains with unique letters
        for chain in structs[-1].get_chains():
            chain.id = string.ascii_uppercase[25 - counter]
            counter += 1
    final_struct = merge_structs(structs)
    if set_bfac is not None:
        for atom in final_struct.get_atoms():
            atom.set_bfactor(set_bfac)
    io.set_structure(final_struct)
    io.save(out_path, select=NotDisorderedOrH())


def split_pdb(pdb_path, outdir, min_nbrs=0, set_bfac=None,
              n_truncations=[], c_truncations=[]):
    """Split a PDB file into one PDB file for each chain.

    Parameters
    ----------
    pdb_path : str
        Path to PDB file to be split into chains.
    outdir : str
        Path to directory in which to output chain PDB files.
    min_nbrs : int, optional
        Minimum number of neighbors of residues to be included in the 
        output structures. The b-factor on the alpha-carbon of a residue 
        is assumed to the the log of the number of neighbors.
    set_bfac : float, optional
        Float value to which to set b-factors in the output PDB.
    n_truncations : list
        List of ints corresponding to numbers of residues to 
        remove from the N-terminus of each chain.
    c_truncations : list
        List of ints corresponding to numbers of residues to 
        remove from the C-terminus of each chain.

    Returns
    -------
    outpaths : list
        Paths to output PDB files.
    """
    struct = get_struct("struct0", pdb_path, min_nbrs)
    # set b-factor of each atom if necessary
    if set_bfac is not None:
        for atom in struct.get_atoms():
            atom.set_bfactor(set_bfac)
    chains = list(struct.get_chains())
    outpaths = []
    for i, chain in enumerate(chains):
        chain.detach_parent()
        outpath = outdir + \
                  '/chain_{}.pdb'.format(string.ascii_uppercase[i])
        outpaths.append(outpath)
        # if necessary, truncate from N- and C-termini
        residues = list(chain.get_residues())
        if len(n_truncations) > i:
            for j, res in enumerate(residues):
                if j < n_truncations[i]:
                    chain.detach_child(res.id)
        if len(c_truncations) > i:
            for j, res in enumerate(residues[::-1]):
                if j < c_truncations[i]:
                    chain.detach_child(res.id)
        io.set_structure(chain)
        io.save(outpath, select=NotDisorderedOrH())
    return outpaths

def gen_loop_query(pdb_paths, out_path, min_nbrs=0):
    """Generate a query PDB for MASTER loop searches given input structures.

    Parameters
    ----------
    pdb_paths : list
        List of paths to the PDB files for the input structures.
    out_path : str
        Path to the output PDB file.
    min_nbrs : int, optional
        Minimum number of neighbors of residues to be included in the 
        output structure. The b-factor on the alpha-carbon of a residue 
        is assumed to the the log of the number of neighbors.

    Returns
    -------
    slice_lengths : list
        Numbers of residues included from the termini of the two input 
        structures in the query PDB.
    """
    assert len(pdb_paths) == 2
    structs = []
    counter = 0
    for i in range(2):
        title = 'struct' + str(i)
        structs.append(get_struct(title, pdb_paths[i], min_nbrs))
        for chain in structs[-1].get_chains():
            chain.id = string.ascii_uppercase[25 - counter]
            counter += 1
    # compute slice lengths by increasing the length until the number of 
    # residues spans an end-to-end distance greater than 10 Angstroms
    res0 = list(structs[0].get_residues())
    res1 = list(structs[1].get_residues())
    xyz_cterm = res0[-1]['CA'].get_coord()
    xyz_nterm = res1[0]['CA'].get_coord()
    slice_lengths = np.array([0, 0])
    idx = 1
    while 0 in slice_lengths:
        if idx < len(res0): 
            xyz = res0[-1 - idx]['CA'].get_coord()
            if np.linalg.norm(xyz - xyz_cterm) > 10.:
                slice_lengths[0] = idx
        else:
            slice_lengths[0] = len(res0)
        if idx < len(res1):
            xyz = res1[idx]['CA'].get_coord()
            if np.linalg.norm(xyz - xyz_nterm) > 10.:
                slice_lengths[1] = idx
        else:
            slice_lengths[1] = len(res1)
        idx += 1
    slices = ['-{}:'.format(str(slice_lengths[0])), 
              ':{}'.format(str(slice_lengths[1]))]
    # merge structures and output PDB file
    final_struct = merge_structs(structs, slices)
    io.set_structure(final_struct)
    io.save(out_path, select=NotDisorderedOrH())
    return slice_lengths


def satisfied_termini(pdb_path, max_nc_dist):
    """Calculate array of of N/C terminus pairs within a given distance.

    Parameters
    ----------
    pdb_path : str
        Path to PDB file containing the structure.
    max_nc_dist : float
        Maximum distance (in Angstroms) between an N- and C-terminus in 
        order for the pair to be considered as "satisfied."

    Returns
    -------
    sat : np.array [n_chains x n_chains]
        Array of 0s and 1s indicating which C- (rows) and N- (columns)
        termini are within max_nc_dist of one another.
    """
    struct = get_struct("struct0", pdb_path)
    chains = list(struct.get_chains())
    n_chains = len(chains)
    # extract positions of N-terminal N atoms and C-terminal C atoms
    n_term = np.zeros((n_chains, 3))
    c_term = np.zeros((n_chains, 3))
    for i, chain in enumerate(chains):
        residues = list(chain.get_residues())
        n_term[i] = residues[0]['N'].get_coord()
        c_term[i] = residues[-1]['C'].get_coord()
    # calculate distances between the N- and C-terminal atoms
    dists = cdist(c_term, n_term)
    # calculate and return the matrix of which distances exceed the threshold
    sat = (dists < max_nc_dist).astype(int)
    return sat


def stitch(pdb_paths, out_path, overlaps=7, min_nbrs=0, 
           seq_replace=None, from_closest=False):
    """Stitch together the structures in a list of 1-chain PDB files.

    Parameters
    ----------
    pdb_paths : list
        List of paths to the PDB files for the structures, organized 
        in an N- to C-terminal order.
    out_path : str
        Path to the output PDB file.
    overlaps : int or list, optional
        List of numbers of overlapping residues between each 
        subsequent pair of structures. If an int is passed, it 
        is expanded into a list of appropriate length.
    min_nbrs : int, optional
        Minimum number of neighbors of residues to be included in the 
        output structure. The b-factor on the alpha-carbon of a residue 
        is assumed to the the log of the number of neighbors.
    seq_replace : list, optional
        List of ints, one for each subsequent pair of structures, 
        denoting which of the two structures should supply the 
        sequence for the final structure (0 for the first, 1 for 
        the second).
    from_closest : bool, optional
        If True, begin linear interpolation from nearest alpha-carbon 
        within the overlap region.

    Returns
    -------
    clashing : bool
        If True, at least two alpha-carbons are within 3 Angstroms.
    """
    if type(overlaps) is int:
        overlaps = [overlaps] * (len(pdb_paths) - 1)
    if not seq_replace:
        # assume SSE's are even and loops are odd, and supply 
        # sequence from the SSE's
        seq_replace = [1] * (len(pdb_paths) - 1)
        seq_replace[::2] = [0] * (len(pdb_paths) // 2)
    structs = [get_struct("struct0", pdb_paths[0], min_nbrs)]
    # iterate over additional structures, adding them to their predecessor
    for i in range(1, len(pdb_paths)):
        title = 'struct' + str(i)
        structs.append(get_struct(title, pdb_paths[i], min_nbrs))
        overlap = overlaps[i - 1]
        res0 = list(structs[i - 1].get_residues())
        res1 = list(structs[i].get_residues())
        # get backbone atom coordinates in overlap region
        xyz0 = np.zeros((overlap * 4, 3))
        xyz1 = np.zeros((overlap * 4, 3))
        for j in range(overlap):
            xyz0[4 * j] = res0[-overlap + j]['N'].get_coord()
            xyz0[4 * j + 1] = res0[-overlap + j]['CA'].get_coord()
            xyz0[4 * j + 2] = res0[-overlap + j]['C'].get_coord()
            xyz0[4 * j + 3] = res0[-overlap + j]['O'].get_coord()
            xyz1[4 * j] = res1[j]['N'].get_coord()
            xyz1[4 * j + 1] = res1[j]['CA'].get_coord()
            xyz1[4 * j + 2] = res1[j]['C'].get_coord()
            xyz1[4 * j + 3] = res1[j]['O'].get_coord()
        if from_closest:
            dist_sq = np.sum((xyz0 - xyz1)[1::4] * 
                             (xyz0 - xyz1)[1::4], axis=1)
            close_idx = np.argmin(dist_sq)
        else:
            close_idx = 0
        # determine coefficients of linear interpolation
        coeffs = np.zeros((overlap * 4, 3))
        for j in range(4):
            for k in range(3):
                coeffs[j::4, k] = np.hstack(
                    [np.zeros(close_idx), 
                     np.linspace(0., 1., overlap - close_idx)])
        xyz = (1. - coeffs) * xyz0 + coeffs * xyz1
        # replace sequence of overlap region with appropriate source
        if seq_replace[i - 1]:
            chain = list(structs[i - 1].get_chains())[0]
            r_ids = list(chain.child_dict.keys())[-overlap:]
            [chain.detach_child(r_id) for r_id in r_ids]
            for j in range(overlap):
                res1[j]['N'].set_coord(xyz[4 * j])
                res1[j]['CA'].set_coord(xyz[4 * j + 1])
                res1[j]['C'].set_coord(xyz[4 * j + 2])
                res1[j]['O'].set_coord(xyz[4 * j + 3])
        else:
            chain = list(structs[i].get_chains())[0]
            r_ids = list(chain.child_dict.keys())[:overlap]
            [chain.detach_child(r_id) for r_id in r_ids]
            for j in range(overlap):
                res0[-overlap + j]['N'].set_coord(xyz[4 * j])
                res0[-overlap + j]['CA'].set_coord(xyz[4 * j + 1])
                res0[-overlap + j]['C'].set_coord(xyz[4 * j + 2])
                res0[-overlap + j]['O'].set_coord(xyz[4 * j + 3])
        # Move residues from structure i - 1 to structure i
        chain0 = list(structs[i - 1].get_chains())[0]
        chain1 = list(structs[i].get_chains())[0]
        init_id = max([res.id[1] for res in chain0]) + 1
        for res_id, res in enumerate(chain1.get_residues()):
            res.detach_parent()
            res.id = (' ', init_id + res_id, ' ')
            chain0.add(res)
        structs = structs[:-1]
        structs.append(structs[-1])
    # calculate minimum distance between C-alphas (adding the 
    # Toeplitz matrix excludes distances of neighboring residues)
    xyz = np.array([at.get_coord() for at in structs[-1].get_atoms()
                    if at.get_name() == 'CA'])
    toep_col = np.roll([1, 1, 1] + [0] * (xyz.shape[0] - 3), -1)
    dists = cdist(xyz, xyz) + toeplitz(10000. * toep_col)
    clashing = np.min(dists) < 3.
    # output structure
    io.set_structure(structs[-1])
    io.save(out_path, select=NotDisorderedOrH())
    return clashing


def calc_compactness(pdb_path):
    """Calculate a measure of protein compactness for a PDB file.

    Parameters
    ----------
    pdb_path : str
        Path to the PDB to be analyzed.

    Returns
    -------
    compactness : float
        The alpha hull volume of the protein divided by the alpha
        hull surface area and the radius of gyration (adjusted by a 
        fitting constant).  Should be 0.170 +/- 0.028 for a compact 
        single-domain protein.
    """
    atoms = pr.parsePDB(pdb_path).select('protein')
    ahull = convex_hull.AlphaHull()
    ahull.set_coords(atoms)
    ahull.calc_hull()
    V = ahull.get_volume()
    A = ahull.get_surface_area()
    r_gyr = pr.measure.measure.calcGyradius(atoms)
    compactness = V / A / (r_gyr + 7.35)
    return compactness


def calc_compactness_convex(pdb_path):
    """Calculate a measure of protein compactness for a PDB file.

    Parameters
    ----------
    pdb_path : str
        Path to the PDB to be analyzed.

    Returns
    -------
    compactness : float
        The convex hull volume of the protein divided by the convex 
        hull surface area and the radius of gyration (adjusted by a 
        fitting constant).  Should be 0.170 +/- 0.028 for a compact 
        single-domain protein.
    """
    atoms = pr.parsePDB(pdb_path).select('protein')
    xyz = atoms.getCoords()
    tri = Delaunay(xyz)
    V, A = tri_v_and_a(xyz, tri.simplices, tri.convex_hull)
    r_gyr = pr.measure.measure.calcGyradius(atoms)
    compactness = V / A / (r_gyr + 7.35)
    return compactness


@numba.jit(nopython=True)
def tri_v_and_a(xyz, simplices, faces):
    V = 0.
    A = 0.
    M = np.zeros((3, 3))
    n_simplices = len(simplices)
    n_faces = len(faces)
    for i in range(max(n_simplices, n_faces)):
        if i < n_simplices:
            M[0] = xyz[simplices[i, 0]] - xyz[simplices[i, 3]]
            M[1] = xyz[simplices[i, 1]] - xyz[simplices[i, 3]]
            M[2] = xyz[simplices[i, 2]] - xyz[simplices[i, 3]]
            V += np.abs(np.linalg.det(M)) / 6.
        if i < n_faces:
            M[0] = xyz[faces[i, 0]] - xyz[faces[i, 2]]
            M[1] = xyz[faces[i, 1]] - xyz[faces[i, 2]]
            A += np.linalg.norm(np.cross(M[0], M[1])) / 2.
    return V, A


def check_clashes(pdb_paths, res_ids_to_keep=[]):
    """Check whether there are steric clashes in a list of protein structures.

    Parameters
    ----------
    pdb_paths : list
        Paths to the PDB files of the structures to be analyzed.
    res_ids_to_keep : list, optional
        List of lists of residue IDs from each PDB file to be fully retained, 
        including side chains, in the structures when clashes are checked.

    Returns
    -------
    clashing : bool
        If True, there are steric clashes among the structures.
    """
    assert len(res_ids_to_keep) in [0, len(pdb_paths)]
    clashing = False
    # check clashes between structures pairwise
    for pair in combinations(pdb_paths, 2):
        atoms0 = pr.parsePDB(pair[0])
        atoms1 = pr.parsePDB(pair[1])
        sel0 = ('name CA C O N H HA HA2 HA2 '
                'CB HB3 HB2 HB1 HB 1HB 2HB 3HB')
        sel1 = sel0
        # if necessary, keep sidechains of key residues
        if len(res_ids_to_keep) == len(pdb_paths):
            idx0 = pdb_paths.index(pair[0])
            idx1 = pdb_paths.index(pair[1])
            min_rn0 = min(atoms0.getResnums())
            min_rn1 = min(atoms1.getResnums())
            ritk0 = [str(val + min_rn0) for val in res_ids_to_keep[idx0]]
            ritk1 = [str(val + min_rn1) for val in res_ids_to_keep[idx1]]
            if len(ritk0) > 0:
                sel0 += ' or resnum {}'.format(' '.join(ritk0))
            if len(ritk1) > 0:
                sel1 += ' or resnum {}'.format(' '.join(ritk1))
        atoms0 = atoms0.select('protein').select(sel0)
        atoms1 = atoms1.select('protein').select(sel1)
        pdb0 = pdb.PDB(atoms0)
        pdb1 = pdb.PDB(atoms1)
        cla = clash.FastClash(pdb0, pdb1, copy_pq=False, noH=False)
        n_res = len(cla.pq.resindices)
        cla.whittle()
        if len(cla.pq.resindices) != n_res:
            clashing = True
            break
    return clashing


def check_gaps(pdb_path):
    """Check if any gaps larger than 2 Angstroms exist in the backbone.

    Parameters
    ----------
    pdb_path : str
        Path to PDB file to be analyzed.

    Returns
    -------
    gaps_exist : bool
        If True, gaps larger than 2 Angstroms exist in the backbone.
    """
    atoms = pr.parsePDB(pdb_path).select('name CA C N')
    xyz = atoms.getCoords()
    dists = cdist(xyz, xyz)
    gaps_exist = False
    for i in range(len(dists) - 1):
        if dists[i, i+1] > 2.:
            gaps_exist = True
    return gaps_exist
