'''
Loop secondary structure in a user defined way. 
'''

import os
import sys
import gzip
import shutil
import numpy as np
import prody as pr
import tarfile

from dataclasses import dataclass

from metalprot.basic import prody_ext

from ..basic import cluster_loops, plot, struct_analysis
from ..master import query, extract_master


@dataclass
class Struct_info:
    loop_info: str
    loop_len: int 
    cent_pdb:str
    clust_num: int   
    min_rmsd: float
    median_rmsd: float
    redundancy:str = "Unknown" 
    clust_num_2nd: int = 0


def loop_search_query_search(loop_workdir, loop_query, loop_range, loop_target_list, rmsdCut, master_query_loop_top, loop_query_len = 14, rm_duplicate = True):
    """find loops with MASTER"""
    gapLen = str(loop_range[0]) + '-' + str(loop_range[1])
    loop_outfile = loop_workdir + 'stdout'
    print('Querying MASTER for loops of length {} to {}.'.format(
            str(loop_range[0]), str(loop_range[1])))
    #print('Database: {}'.format(loop_target_list))
    query.master_query_loop(loop_query, loop_target_list, 
                            rmsdCut=rmsdCut, topN=master_query_loop_top,
                            gapLen=gapLen, outdir=loop_workdir, 
                            outfile=loop_outfile)

    print('Sorting loop PDBs by loop length.' + loop_workdir)
    # sort PDBs into directories by loop length.
    # rm duplicates.
    seq_set = dict()
    for path in os.listdir(loop_workdir):
        if '.pdb' in path and ('match' in path or 'wgap' in path):
            # with open(loop_workdir + path, 'r') as f:
            #     res_ids = set([int(line[23:26]) for line in f.read().split('\n') if line[:4] == 'ATOM'])
            #     # subtract query ends from loop length
            #     l = len(res_ids) - loop_query_len
            pdb = pr.parsePDB(loop_workdir + path)
            seq = ''.join(pdb.select('name CA').getSequence())

            if seq in seq_set.keys():
                seq_set[seq].append(path)
            else:
                seq_set[seq] = [path]

            if rm_duplicate and len(seq_set[seq]) > 1:
                os.remove(loop_workdir + path)
                continue

            l = len(pdb.select('name CA')) - loop_query_len          
            l_dir = loop_workdir + str(l)
            # create a directory for the loop length if necessary
            #print(l_dir)
            os.makedirs(l_dir, exist_ok= True)
            os.rename(loop_workdir + path, l_dir + '/' + os.path.basename(path))

    ### write duplicate summary.
    if rm_duplicate:
        with open(loop_workdir + '_dup_summary.txt', 'w') as f:
            for k in seq_set.keys():
                f.write(k + '\t' + '\t'.join(seq_set[k]) + '\n')
                
    return  


def ext_sel(dir, target, lo):
    '''
    lo = ('A', 18, 24, 'B', 1, 7)
    # Chain A loop B. [(A, pos1, pos2, B, pos1, pos2)]
    '''
    resnum_a = list(range(lo[1],  lo[2] + 1, 1))
    a = target.select('chid ' + lo[0] + ' and resnum ' + ' '.join([str(x) for x in resnum_a]))

    resnum_b = list(range(lo[4], lo[5] + 1, 1))
    b = target.select('chid ' + lo[3] + ' and resnum ' + ' '.join([str(x) for x in resnum_b]))

    name = '-'.join([str(x) for x in lo])
    ab = prody_ext.combine_ags_into_one_chain([a, b], title = name)

    pr.writePDB(dir + name + '.pdb', ab)

    lo_query_len = len(resnum_a) + len(resnum_b)
    return lo_query_len


def ext_lo_query_len(lo_name):
    '''
    lo = ('A', 18, 24, 'B', 1, 7)
    '''
    lo = lo_name.split('-')
    resnum_a = list(range(int(lo[1]), int(lo[2]) + 1, 1))
    resnum_b = list(range(int(lo[4]), int(lo[5]) + 1, 1))
    lo_query_len = len(resnum_a) + len(resnum_b)
    return lo_query_len


def _get_top_cluster_summary(topo_dir, cluster_count_cut, loop_range, select_min_rmsd_pdb):
    '''
    After search & cluster loops, calculate the clustered loop info. 
    Copy the centroid or min-rmsd loop into self.workdir/loops_X_Y folder for futher consideration.
    Plot the phi/psi, sequence logo, hydrophobicity.
    '''
    _infos = []
    for _lo_dir_name in os.listdir(topo_dir):
        lo_dir = topo_dir + _lo_dir_name + '/'
        if not os.path.isdir(lo_dir):
            continue

        for l in range(loop_range[0], loop_range[1] + 1):
            subdir = lo_dir + '{}/clusters/1'.format(str(l))    
                        
            try:
                loop_pdbs = os.listdir(subdir)
            except:
                loop_pdbs = []
            #print(subdir + '---' + str(len(loop_pdbs)))

            if len(loop_pdbs) > 1:       
                loop_rmsds, loop_seqs, loop_pdss = extract_master._get_pdbs_master_info(lo_dir + '/match.txt', lo_dir + '/seq.txt', loop_pdbs)
                _cent_pdb = ''  

                #Copy centroid pdb and plot   
                if len(loop_pdbs) >= cluster_count_cut:

                    _cent_dir = topo_dir[:-1] + '_cent_info/'
                    if not os.path.exists(_cent_dir):
                        os.mkdir(_cent_dir)
                    _cent_name = _lo_dir_name + '_cent_rg_' + str(l) + '_clu_' + str(len(loop_pdbs))
                    _cent_pdb = _cent_dir + _cent_name + '.pdb'

                    if select_min_rmsd_pdb:
                        in_pdb = subdir + '/' + loop_pdbs[np.argmin(loop_rmsds)]
                    else:
                        in_pdb = subdir + '/' + [lpdb for lpdb in loop_pdbs if 'centroid' in lpdb][0]

                    with gzip.open(in_pdb, 'rb') as f_in:
                        with open(_cent_pdb, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)                  

                    phi, psi, _sel_seq = struct_analysis.cal_phipsi(_cent_pdb)
                    lo_query_len = ext_lo_query_len(_lo_dir_name)
                    plot._plot_all(_cent_dir + _cent_name, loop_seqs, loop_rmsds, l, lo_query_len, phi, psi, _sel_seq)     

                #Add summary                
                info = Struct_info(loop_info = _lo_dir_name, loop_len = l, clust_num = len(loop_pdbs), min_rmsd = min(loop_rmsds), median_rmsd = np.median(loop_rmsds), cent_pdb = _cent_pdb)              
                
                _infos.append(info)      

    return _infos

def run_loop_ss(outdir, target_file, loop_topo_sels, para):
    os.makedirs(outdir, exist_ok=True)
    target = pr.parsePDB(target_file)

    for topo in loop_topo_sels.keys():
        loop_sels = loop_topo_sels[topo]

        #Create topo_dir
        topo_dir = outdir + topo + '/'
        os.makedirs(topo_dir, exist_ok=True)

        for chidchid in loop_sels.keys():
            loops = loop_sels[chidchid]
            for lo in loops:
 
                #Create lo_dir
                name = '-'.join([str(x) for x in lo])
                lo_dir = topo_dir + name + '/'
                os.makedirs(lo_dir, exist_ok=True)
                # Extract sel lo.
                lo_query_len = ext_sel(lo_dir, target, lo)
                lo_query = lo_dir + name + '.pdb'
                loop_search_query_search(lo_dir, lo_query, para.loop_range, para.loop_target_list, para.rmsdCut, para.master_query_loop_top, loop_query_len = lo_query_len, rm_duplicate = para.rm_duplicate)
                log = ''
                cluster_loops.run_cluster(lo_dir + '/', log, lo_query_len,  para.cluster_rmsd , outfile= lo_dir + 'stdout')


        _get_top_cluster_summary(topo_dir, para.cluster_count_cut, para.loop_range, para.select_min_rmsd_pdb)

        #zip the file and delete the unzipped pdbs.
        tar = tarfile.open(outdir + topo + ".tar.gz", "w:gz")
        tar.add(topo_dir, arcname=topo)
        tar.close()
        shutil.rmtree(topo_dir)

    return 



