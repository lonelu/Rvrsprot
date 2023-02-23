import numpy as np
import prody as pr
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
import pandas as pd
import re

'''
The script specially work for helix bundles.
'''

def combine_ags_into_one_chain(ags, title):
    '''
    Same in Metalprot.basic.combine_ags_into_one_chain.
    '''
    ag = pr.AtomGroup(title)
    coords = []
    chids = []
    names = []
    resnames = []
    resnums = []
    resnum = 1
    for _ag_all in ags:
        for cd in np.unique(_ag_all.getChids()):
            #print(cd)
            chid = 'A'

            if cd == None:
                _ag = _ag_all
            else:
                _ag = _ag_all.select('chid ' + cd)
                
            for i in np.unique(_ag.getResindices()):
                c = _ag.select('resindex ' + str(i))
                coords.extend(c.getCoords())
                chids.extend([chid for x in range(len(c))])
                names.extend(c.getNames())
                resnames.extend(c.getResnames())
                resnums.extend([resnum for x in range(len(c))])
                resnum += 1

    ag.setCoords(np.array(coords))
    ag.setChids(chids)
    ag.setNames(names)
    ag.setResnames(resnames)
    ag.setResnums(resnums)
    return ag


def cal_dist_info(target, loop):
    '''
    Generate dist map infomation target and loop.
    The current dist is not exactly real radius.
    '''

    target_coords = target.select('bb and name N CA C').getCoords()
    target_coords_r = target_coords.reshape(int(target_coords.shape[0]/3), 9)
    loop_coords = loop.select('bb and name N CA C').getCoords()
    loop_coords_r = loop_coords.reshape(int(loop_coords.shape[0]/3), 9)

    neigh_y = NearestNeighbors(radius= 1.5) 
    neigh_y.fit(target_coords_r)

    x_in_y = neigh_y.radius_neighbors(loop_coords_r)

    chidnames = target.select('bb and name CA').getChids()
    resnums = target.select('bb and name CA').getResnums()
    resnames = target.select('bb and name CA').getResnames()
    target_ids = [(chidnames[id[0]], resnums[id[0]], resnames[id[0]]) if len(id) > 0 else ('', '', '') for id in x_in_y[1]]
    
    
    loopchidnames = loop.select('bb and name CA').getChids()
    loopresnums = loop.select('bb and name CA').getResnums()
    loopresnames = loop.select('bb and name CA').getResnames()
    loop_ids = [(loopchidnames[i], loopresnums[i], loopresnames[i]) for i in range(len(loopchidnames))]
    
    dists = [dist[0]  if len(dist) > 0 else -1 for dist in x_in_y[0]]
    return target_ids, loop_ids, dists


def plot_dist_info(target_ids, loop_ids, dists, filepath):
    '''

    '''
    #fig, (ax) =plt.subplots(1, 1, figsize=(6, 2))
    fig = plt.figure()
    ax = fig.add_subplot()
    data = []
    data.append(target_ids)
    data.append(loop_ids)
    data.append([round(f, 2) for f in dists])

    ax.set_axis_off()
    df=pd.DataFrame(data).T
    ax.axis('tight')
    ax.axis('off')
    cols = ['target', 'loop', 'dist']
    tab = ax.table(cellText=df.values, colLabels= cols, rowLabels=list(range(1, df.shape[0]+1)), cellLoc='center', loc='upper left')
    
    tab.set_fontsize(10)
    #tab.scale(1, 2)
    tab.auto_set_font_size(False)
    ax.set_ylabel("Legend", fontsize = 12)

    fig.tight_layout()
    plt.tight_layout()
    plt.savefig(filepath+'_info.png')
    plt.close()
    return


def auto_pick_pos(target_ids, loop_ids, dists, front =7, back = 7):
    '''
    
    '''
    front_min = back_min = 2.0
    front_min_id = back_min_id = 0

    for i in range(front):
        if dists[i] >= 0 and dists[i] < front_min:
            front_min = dists[i]
            front_min_id = i
    for i in range(len(dists)-back, len(dists)):
        if dists[i] >= 0 and dists[i] < back_min:
            back_min = dists[i]
            back_min_id = i

    target_pos = (target_ids[front_min_id], target_ids[back_min_id])
    loop_pos = (loop_ids[front_min_id+1], loop_ids[back_min_id-1])

    return target_pos, loop_pos


def auto_generate_sels(target, loops, outdir, target_start, target_end, front=7, back=7):
    '''
    Calculate the distance between N CA C, extract the position with min distance of N CA C. 
    '''
    structs = []
    sels = []
    if target_start == '':
        target_sel = target.select('bb and name CA')[0]
        sels.append((target_sel.getChid(), target_sel.getResnum(), target_sel.getResname()))
    else:
        target_starts = target_start.split(',')
        sels.append((target_starts[0], target_starts[1], target_starts[2]))
    for loop in loops:
        structs.append(target)
        structs.append(loop)
        target_ids, loop_ids, dists = cal_dist_info(target, loop)

        plot_dist_info(target_ids, loop_ids, dists, outdir + loop.getTitle())

        target_pos, loop_pos = auto_pick_pos(target_ids, loop_ids, dists, front, back)
        sels.append(target_pos[0])
        sels.append(loop_pos[0])
        sels.append(loop_pos[1])
        sels.append(target_pos[1])

    structs.append(target)
    if target_end == '':
        target_sel = target.select('bb and name CA')[-1]
        sels.append((target_sel.getChid(), target_sel.getResnum(), target_sel.getResname()))
    else:
        target_ends = target_end.split(',')
        sels.append((target_ends[0], target_ends[1], target_ends[2]))

    sels_order = [(sels[i*2], sels[i*2+1]) for i in range(len(loops)*2+1)]
    return structs, sels_order


def generate_ags(structs, sels_order):
    '''
    
    '''
    ags = []
    for i in range(0, len(structs)):
        #TO DO: generate selection.
        sel_str = 'protein and chid ' + sels_order[i][0][0] + ' and resnum ' + str(sels_order[i][0][1]) + 'to' + str(sels_order[i][1][1])
        ag = structs[i].select(sel_str)
        ags.append(ag)
    return ags


def connect_struct(outdir, title, targetpath, looppaths, target_start = '', target_end = ''):
    '''
    
    '''
    target = pr.parsePDB(targetpath)
    loops = [pr.parsePDB(looppath) for looppath in looppaths]

    structs, sels = auto_generate_sels(target, loops, outdir, target_start, target_end)
        
    ags = generate_ags(structs, sels)

    combined_ag = combine_ags_into_one_chain(ags, title)

    pr.writePDB(outdir + title, combined_ag)


def get_user_sels(user_sel, NumberOfloops):
    #all_sels = user_sel.strip()
    all_sels = user_sel.replace(' ','').strip().split('\n')
    #print(all_sels)
    sels = []
    for i in range(2*NumberOfloops + 1):
        sels.append(((all_sels[i].split(',')[1], all_sels[i].split(',')[2]), (all_sels[i].split(',')[3], all_sels[i].split(',')[4])))
    return sels


def print_auto_sels(structs, sels, NumberOfloops):
    auto_sels = ''
    for i in range(2*NumberOfloops + 1):
        #print(structs[i].getTitle())
        #print(sels[i])
        auto_sels += structs[i].getTitle() + ',' + sels[i][0][0] + ',' + str(sels[i][0][1]) + ',' + sels[i][1][0] + ',' + str(sels[i][1][1]) + '\n'
    print(auto_sels)
    return auto_sels

def test():
    test_dir = '/mnt/e/DesignData/Metalloprotein/SAHA_Vorinostat/run_design_cgs3/parametric_bundles/param_ala/loop_connect/'
    #target = pr.parsePDB(test_dir + '00009.f63440efff7e.allbb_ala_min_ala_0001.pdb') 
    #loop = pr.parsePDB(test_dir + 'loop_ss_20220721-224244/ABC_frt/_cent_/A-30-36-A-39-45_cent_rg_5_clu_121.pdb')
    targetpath = test_dir + '00009.f63440efff7e.allbb_ala_min_ala_0001.pdb'
    looppaths = [
        test_dir + 'A-104-110-A-113-119_cent_rg_5_clu_66.pdb',
        test_dir + 'A-138-144-A-5-11_cent_rg_4_clu_14.pdb',
        test_dir + 'A-30-36-A-39-45_cent_rg_5_clu_121.pdb'
    ]

    outdir = test_dir
    title = 'combined'
    target_start = 'A,75,ALA'
    target_end = 'A,74,ALA'
    connect_struct(outdir, title, targetpath, looppaths, target_start, target_end)


