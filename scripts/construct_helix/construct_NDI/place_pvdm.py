import os
import prody as pr
import numpy as np
from metalprot.basic import prody_ext
from sklearn.neighbors import NearestNeighbors
from rvrsprot.basic import convex_hull

'''
python /mnt/e/GitHub_Design/Rvrsprot/scripts/construct_helix/construct_NDI/place_pvdm.py
'''

def clash(Acoords, Bcoords, vdm_vdm_clash_dist = 2.6):
    neigh_y = NearestNeighbors(radius= vdm_vdm_clash_dist)
    neigh_y.fit(Acoords)
    x_in_y = neigh_y.radius_neighbors(Bcoords)
    x_has_y = any([True if len(a) >0 else False for a in x_in_y[1]])
    if x_has_y:
        return True
    return False 

def calc_ahull(target):
    ahull = convex_hull.AlphaHull(alpha=9)  # alpha sphere diameter is 9 Angstroms
    ahull.set_coords(target) # Can set coordinates from any
                                # arbitrary array (arr) by ahull.coords = arr
    ahull.calc_hull()
    return ahull

def in_hull(ahull, query_coords):
    in_hull = ahull.pnts_in_hull(query_coords)
    return in_hull

def get_resind_cand(target, pvdm0):
    resind_candidates = []

    #ahull_query_coords = []
    xs = target.select('within 9 of name C19 C20').getResindices()
    for resind in np.unique(target.select('protein and chid B').getResindices()):
        if resind in xs:
            continue
        pvdm0_cp = pvdm0.copy()
        pr.calcTransformation(pvdm0_cp.select('resindex 1 and name N CA C'), target.select('name N CA C and resindex ' + str(resind))).apply(pvdm0_cp)

        if not clash(pvdm0_cp.select('resname SMU and heavy and not name C21 C22 C23 C24 C25 C26').getCoords(), target.select('resname NDI and heavy').getCoords(), 20.0):
            continue
        resind_candidates.append(resind)
        #ahull_query_coords.append(pvdm0_cp.select('name MN1').getCoords()[0])

    ### The ahull is not really working great here.
    #ahull = calc_ahull(target)
    #in_hull = calc_ahull(ahull, ahull_query_coords)
    #resind_candidates = [resind_candidates[x] for x in range(len(resind_candidates)) if in_hull[x]]

    return resind_candidates



def place_vdm(target, pvdm_dir, outdir):
    os.makedirs(outdir, exist_ok= True)

    ahull = calc_ahull(target)

    pvdms = []
    for file in os.listdir(pvdm_dir):
        if not '.pdb' in file:
            continue
        if float(file.split('_')[4][0:-4]) < -1.0:
                continue
        pvdms.append(pr.parsePDB(pvdm_dir + file))
    pvdm0 = pvdms[0]

    #use the first pvdm to find the ok position for the porphyrin.  
    resind_candidates = get_resind_cand(target, pvdm0)
    
    ag = target.select('resindex ' + ' '.join([str(x) for x in resind_candidates]))
    pr.writePDB(outdir + '_test5.pdb', ag)
    
    for resind in resind_candidates:
        for pvdm in pvdms:
            pvdm_cp = pvdm.copy()

            print(resind)
            print(pvdm_cp.getTitle())
            pr.calcTransformation(pvdm_cp.select('resindex 1 and name N CA C'), target.select('name N CA C and resindex ' + str(resind))).apply(pvdm_cp) 
            
            #pr.writePDB(outdir + 'resind_'+ str(resind) + '_' + pvdm_cp.getTitle(), pvdm_cp)

            if clash(pvdm_cp.select('name CG ND1 CD2 CE1 NE2').getCoords(), target.select('name N C CA O').getCoords(), 2.5):
                continue

            if clash(pvdm_cp.select('resname SMU and heavy and not name C21 C22 C23 C24 C25 C26').getCoords(), target.select('resname NDI and heavy').getCoords(), 3.6):
                continue

            if clash(pvdm_cp.select('resname SMU and heavy and not name C21 C22 C23 C24 C25 C26').getCoords(), target.select('bb').getCoords(), 2.5):
                continue

            if not clash(pvdm_cp.select('resname SMU and heavy').getCoords(), target.select('resname NDI and heavy').getCoords(), 6.6):
                continue

            ahull_query_coords =pvdm_cp.select('name MN1').getCoords()
            if not in_hull(ahull, ahull_query_coords)[0]:
                continue

            pr.writePDB(outdir + 'resind_'+ str(resind) + '_' + pvdm_cp.getTitle(), pvdm_cp)

    return 

def main0():
    workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/SamM/'
    target = pr.parsePDB(workdir + 'NDI_20_pose107.pdb')

    pvdm_dir = workdir + 'porphyrin_his_vdm/'
    outdir = workdir + 'vdm_candidate_placement_33/'
    place_vdm(target, pvdm_dir, outdir)
    return 

def main():
    workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/SamM/'
    indir = workdir + 'new_bb/'
    for file in os.listdir(indir):
        if not '.pdb' in file:
            continue
        target = pr.parsePDB(indir + file)
        pvdm_dir = workdir + 'porphyrin_his_vdm/'
        outdir = indir + file.split('.')[0] + '/'
        place_vdm(target, pvdm_dir, outdir) 
    return

if __name__=='__main__':
    main()