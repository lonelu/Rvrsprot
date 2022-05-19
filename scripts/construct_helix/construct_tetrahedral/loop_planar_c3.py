'''
Using the best vdMs, one could generate planar C3. 
Then Master Search the loops.
'''
import os
from rvrsprot.loop import loop_ss
import multiprocessing as mp

'''
python /mnt/e/GitHub_Design/Rvrsprot/scripts/construct_helix/construct_tetrahedral/loop_planar_c3.py
'''

class para:
    def  __init__(self) -> None:
        self.loop_topo_sels = {
            # Toplogy
            'ABC_frt' : {
                'A-B': [('A', 8, 2, 'B', 8, 2)], # Chain A loop B. [(A, pos, +-len, B, pos, +-len)]
                'A-B': [('A', 8, 2, 'C', 8, 2)],
            },
        }

        ### Master search para
        self.loop_range = [3, 15]
        self.loop_target_list = '/mnt/e/DesignData/Database/Qbits/pds_list_2p5.txt' # gpu:'/mnt/e/GitHub_Design/master_db/list'
        self.rmsdCut = 0.6
        self.master_query_loop_top = 500
        self.cluster_rmsd = 1.0
        self.rm_duplicate = True

        ### Output para
        self.cluster_count_cut = 10
        self.select_min_rmsd_pdb= False


def main():
    workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c3/c3_vc0/'

    files = os.listdir(workdir)
    _para = para() 

    pool = mp.Pool(6)
    [pool.apply_async(loop_ss.run_loop_ss, args=(workdir + 'loop_ss_' + file.split('.')[0] + '/', workdir + file, _para.loop_topo_sels, _para)) for file in files]
    pool.close()
    pool.join()

    return


if __name__=='__main__':
    main()
