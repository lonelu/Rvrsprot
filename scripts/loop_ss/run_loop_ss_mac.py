import os
import sys
from rvrsprot.loop import loop_ss
import datetime


'''
python /Users/lonelu/GitHub_Design/Rvrsprot/scripts/loop_ss/run_loop_ss_mac.py 

'''

class Para():
    workdir = '/Users/lonelu/DesignData/LigandBB/BASF_input_backbones/bb_ALA/'
    #outdir = 'loop_ss_20220401-164040/'
    outdir = None
    target_file = '1f4m_CCCP_ALA.pdb'

    loop_topo_sels = {
        # Toplogy
        'ABC_frt' : {
            'A-B': [('A', 36, 42, 'B', 17, 23)], # Chain A loop B. [(A, pos, +-len, B, pos, +-len)]
        },
    }

    ### Master search para
    loop_range = [12, 13]
    loop_target_list = '/Users/lonelu/Database/master_db/list'
    rmsdCut = 0.6
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