## activate the conda environment. conda activate env_smallprot
## please check/change the parameters in 'parameter_loop_truc.ini' 
## you can use ipython or just 'python run.py'
import sys
sys.path.append(r'/mnt/e/GitHub_Design/Qbits')
sys.path.append(r'/mnt/e/GitHub_Design/smallprot')

print('Thanks for using smallprot!')

from smallprot import smallprot, smallprot_config

para = smallprot_config.Parameter(
        num_iter = 3, 
        top = 5,      
        master_query_top = 200,
        screen_compactness = False,
        rmsdCut = 1.0,
        qbits_rmsd = 1.5,
        qbits_window = 10,
        secstruct = None,
        min_nbrs = 1,     
        lowest_rmsd_loop = True,
        ###Database
        database='/mnt/e/GitHub_Design/Qbits/database',
        loop_target_list='/mnt/e/GitHub_Design/master_db/list',   
        ###For loop searching     
        master_query_loop_top = 200,
        max_nc_dist = 16.0,
        loop_query_win =7,   
        min_loop_length = 3,
        max_loop_length=20,
        cluster_count_cut=20,
        loop_distance_cut=15,
        construct_keep = 0
)


seed_pdb = '/mnt/e/DesignData/nina/seed.pdb'
query_pdb = None
exclusion_pdb = None

workdir = '/mnt/e/DesignData/nina/output_test_build1/'


hhh = smallprot.SmallProt(seed_pdb, query_pdb, exclusion_pdb,  workdir, para)
# n_truncations=list(range(20))
# c_truncations=list(range(10))
# n_truncations = [16, 17, 18, 19]
# c_truncations = [1, 2, 3, 4, 5, 6, 7, 8, 9]
n_truncations = [16]
c_truncations = [1]
direction=[2, 1, 0, 3]
# direction = []
hhh.loop_structure(direction = direction, n_truncations = n_truncations, c_truncations = c_truncations)