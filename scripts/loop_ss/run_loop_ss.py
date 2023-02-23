import os
import sys
from rvrsprot.loop import loop_ss
import datetime


'''
python /mnt/e/GitHub_Design/Rvrsprot/scripts/loop_ss/run_loop_ss.py

'''

class Para():
    workdir = '/mnt/e/DesignData/bpp_fluo_comb/fluo/output1_09_f63440_nick_ala/sel/'
    #outdir = 'loop_ss_20220401-164040/'
    outdir = None
    target_file = 'clean.pdb'


    # loop_topo_sels = {
    #     # Toplogy
    #     'ABC_frt' : {
    #         'A-B': [('A', 21, 3, 'B', 4, 3), ('A', 18, 3, 'B', 7, 3)], # Chain A loop B. [(A, pos, +-len, B, pos, +-len)]
    #         'B-C': [('B', 21, 3, 'C', 6, 3), ('B', 18, 3, 'C', 9, 3)],
    #         'C-D': [('C', 21, 3, 'D', 5, 3), ('C', 18, 3, 'D', 8, 3)],
    #         'D-E': [('D', 22, 3, 'E', 10, 3)],
    #         'E-F': [('E', 21, 3, 'F', 4, 3), ('E', 18, 3, 'F', 7, 3)],
    #         'F-A': [('F', 21, 3, 'A', 4, 3), ('F', 18, 3, 'A', 7, 3)],
    #     },
    #     'ABC_rev' : {
    #         'A-F': [('A', 21, 3, 'F', 5, 3), ('A', 18, 3, 'F', 8, 3)],
    #         'F-E': [('F', 21, 3, 'E', 4, 3), ('F', 18, 3, 'E', 7, 3)],
    #         'E-D': [('E', 21, 3, 'D', 10, 3)],
    #         'D-C': [('D', 21, 3, 'C', 6, 3), ('D', 18, 3, 'C', 9, 3)],
    #         'C-B': [('C', 21, 3, 'B', 6, 3), ('C', 18, 3, 'B', 9, 3)],
    #         'B-A': [('B', 21, 3, 'A', 4, 3), ('B', 18, 3, 'A', 7, 3)],
    #     },
    # }

    loop_topo_sels = {
        # Toplogy
        'ABC_frt' : {
            'A-B': [('A', 15, 21, 'A', 130, 136), ('A', 18, 24, 'A', 127, 133)], # Chain A loop B. [(A, pos, +-len, B, pos, +-len)]
        },
    }

    ### Master search para
    loop_range = [3, 18]
    #loop_target_list = '/mnt/e/DesignData/Database/Qbits/pds_list_2p5.txt' # 
    loop_target_list = '/mnt/e/DesignData/Database/master_db/list'
    rmsdCut = 1.5
    master_query_loop_top = 500
    cluster_rmsd = 1.0
    rm_duplicate = True

    ### Output para
    cluster_count_cut = 10
    select_min_rmsd_pdb= False



def main():
    para = Para()
    time_tag = datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 
    outdir = para.workdir + 'loop_ss_' + time_tag + '/'
    if para.outdir:
        outdir = para.workdir + para.outdir
    os.makedirs(outdir, exist_ok=True)
    loop_ss.run_loop_ss(outdir, para.workdir + para.target_file, para.loop_topo_sels, para)

    return

if __name__=='__main__':
    main()