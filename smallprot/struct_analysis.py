import prody as pr
import qbits

def meature_phipsi(structpath):
    protein = pr.parsePDB(structpath)
    seq = []
    phi_180 = []
    psi_180 = []
    for p in protein.iterResidues():
        seq.append(p.getResname())
        try:
            phi_180.append(prody.calcPhi(p))
        except:
            phi_180.append(None)
        try:
            psi_180.append(prody.calcPsi(p))
        except:
            psi_180.append(None)
    return phi_180, psi_180, seq

    # print(structpath)
    # structname = structpath.split('/')[-1]
    # struct = parser.get_structure(structname, structpath)
    # pp = ppb.build_peptides(struct)[0] #there is a 'bug' here, it cannot read the whole chain sometimes.
    # print(pp.get_sequence())
    # phipsi = pp.get_phi_psi_list()
    # phipsi_180 = []
    # for hs in phipsi:
    #     if hs[0] == None:
    #         phipsi_180.append(0)  
    #     else:
    #         phipsi_180.append(hs[0]/3.14159265*180)

    #     if hs[1] == None:
    #         phipsi_180.append(0)
    #     else:  
    #         phipsi_180.append(hs[1]/3.14159265*180)
    # return phipsi_180


def get_in_ahull_ratio(_full_pdb, alpha=5):
    '''
    get the ahull of the pdb and calculate the ratio of ca that are in the ahull.
    '''
    prody_pdb = pr.parsePDB(_full_pdb)
    ahull = qbits.convex_hull.AlphaHull(alpha)
    ahull.set_coords(prody_pdb)
    ahull.calc_hull()
    ca_sel = prody_pdb.select('name CA')       
    ca_coords = ca_sel.getCoords()
    ca_in_hull = ahull.pnts_in_hull(ca_coords)
    return sum(ca_in_hull)/len(ca_in_hull)

