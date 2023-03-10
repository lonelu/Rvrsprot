'''
Rvrsprot loops connect auto.
'''
import os
import sys
import prody as pr
from rvrsprot.loop import connect_loops
import datetime
import re


'''
conda activate ...
python /mnt/e/GitHub_Design/Rvrsprot/scripts/loop_ss/run_loop_insert.py
'''

class Para():
    workdir = '/mnt/e/DesignData/bpp_fluo_comb/fluo/output1_09_f63440_nick_ala/sel/loop/'
    outdir = None

    target_file = 'bb_prep_unloop.pdb'
    #The loop_files need to be ordered from N-terminal to C-terminal.
    loop_files = ['A-62-68-A-5-11_cent_rg_13_clu_90.pdb',
                    'A-18-24-A-127-133_cent_rg_4_clu_23.pdb', 
                    'A-138-142-A-78-84_cent_rg_12_clu_108.pdb']

    title ='combined'  #The output pdb name.
    NumberOfloops = 3
    target_start ='A,47,GLY'
    target_end='A,101,GLY'

    user_define_connection=True 
    #The user_sel must be split by line '\n' or ';' or '.'
    user_sel ='''
bb_prep_unloop,A,47,A,65
A-62-68-A-5-11_cent_rg_13_clu_90,A,65,A,86
bb_prep_unloop,A,11,A,23
A-18-24-A-127-133_cent_rg_4_clu_23,D,560,D,569
bb_prep_unloop,A,132,A,143
A-138-142-A-78-84_cent_rg_12_clu_108,A,472,A,483
bb_prep_unloop,A,79,A,101
    '''

def main():
    para = Para()
    time_tag = datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 
    outdir = para.workdir + 'loop_insert_' + time_tag + '/'
    if para.outdir:
        outdir = para.workdir + para.outdir
    os.makedirs(outdir, exist_ok=True)

    target = pr.parsePDB(para.workdir + para.target_file)
    loops = [pr.parsePDB(para.workdir + loop_file) for loop_file in para.loop_files]

    structs, sels = connect_loops.auto_generate_sels(target, loops, outdir, para.target_start, para.target_end)

    auto_sels_str = connect_loops.print_auto_sels(structs, sels, para.NumberOfloops)
    with open(outdir + 'auto_sel.txt', 'w') as f:
        f.write(auto_sels_str)

    if para.user_define_connection:
        sels = connect_loops.get_user_sels(para.user_sel, para.NumberOfloops)
        with open(outdir + 'user_sel.txt', 'w') as f:
            f.write(para.user_sel)

    ags = connect_loops.generate_ags(structs, sels)
    combined_ag = connect_loops.combine_ags_into_one_chain(ags, para.title)
    pr.writePDB(outdir + para.title, combined_ag)
    return

if __name__=='__main__':
    main()
