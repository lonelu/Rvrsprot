import os
import prody as pr
from . import pdbutils, convex_hull

def cal_phipsi(pdb_path):
    '''
    calculate phi psi angle for each amino acids.
    return [phi], [psi], [aa]
    '''
    protein = pr.parsePDB(pdb_path)
    seq = []
    phi_180 = []
    psi_180 = []
    for p in protein.iterResidues():
        seq.append(p.getResname())
        try:
            phi_180.append(pr.calcPhi(p))
        except:
            phi_180.append(None)
        try:
            psi_180.append(pr.calcPsi(p))
        except:
            psi_180.append(None)
    return phi_180, psi_180, seq

def cal_phipsi_depre(pdb_path):
    print(pdb_path)
    structname = pdb_path.split('/')[-1]
    struct = pdbutils.parser.get_structure(structname, pdb_path)
    pp = pdbutils.ppb.build_peptides(struct)[0] #there is a 'bug' here, it cannot read the whole chain sometimes.
    print(pp.get_sequence())
    phipsi = pp.get_phi_psi_list()
    phipsi_180 = []
    for hs in phipsi:
        if hs[0] == None:
            phipsi_180.append(0)  
        else:
            phipsi_180.append(hs[0]/3.14159265*180)

        if hs[1] == None:
            phipsi_180.append(0)
        else:  
            phipsi_180.append(hs[1]/3.14159265*180)
    return phipsi_180

def cal_ahull(_full_pdb, alpha=5):
    '''
    get the ahull of the pdb and calculate the ratio of ca that are in the ahull.
    '''
    prody_pdb = pr.parsePDB(_full_pdb)
    ahull = convex_hull.AlphaHull(alpha)
    ahull.set_coords(prody_pdb)
    ahull.calc_hull()
    volume = ahull.get_volume()
    surface_area = ahull.get_surface_area()

    ca_sel = prody_pdb.select('name CA')       
    ca_coords = ca_sel.getCoords()
    ca_in_hull = ahull.pnts_in_hull(ca_coords)
    return sum(ca_in_hull)/len(ca_in_hull), volume, surface_area

def cal_ca_dihedral(pdb_path):
    '''
    Calculate Ca dihedral from the four atoms: Ca_i+2, Ca_i+1, Ca_i, Ca_i-1.
    return [dihedral], [aa]
    '''
    protein = pr.parsePDB(pdb_path)
    all_ca = protein.select('name CA')

    res_ids = list(range(1, all_ca.numAtoms() - 2))
    dihedral = []
    seq = []
    for i in res_ids:
        d = pr.calcDihedral(all_ca[i-1], all_ca[i], all_ca[i+1], all_ca[i+2])
        dihedral.append(d)
        seq.append(all_ca[i].getResname())
    return res_ids, seq, dihedral

def compare_ca_dihedral(pdb_path1, pdb_path2):
    res_ids1, seq1, dihedral1 = cal_ca_dihedral(pdb_path1)
    res_ids2, seq2, dihedral2 = cal_ca_dihedral(pdb_path2)

    dihedral_diff = [a - b for a, b in zip(dihedral1, dihedral2)]

    filename = os.path.join(os.path.dirname(pdb_path1), 'dihedral_comparison.txt')
    with open(filename, 'w') as f:
        for r in res_ids1:
            f.write(str(r) + '\t')
        f.write('\n')
        for r in res_ids2:
            f.write(str(r) + '\t')
        f.write('\n')
        for s in seq1:
            f.write(s + '\t')
        f.write('\n')
        for s in seq2:
            f.write(s + '\t')
        f.write('\n')
        for d in dihedral1:
            f.write(str(d) + '\t')
        f.write('\n')
        for d in dihedral2:
            f.write(str(d) + '\t')
        f.write('\n')





