## activate the conda environment. conda activate env_smallprot 
## you can use ipython or just 'python run.py'
import sys
sys.path.append(r'/mnt/e/GitHub_Design/Qbits')
sys.path.append(r'/mnt/e/GitHub_Design/smallprot')

print('Thanks for using smallprot!')

from smallprot import extend_sse, smallprot_config

para = smallprot_config.Parameter(
        num_iter = 2, 
        top = 5,      
        master_query_top = 500,
        screen_compactness = False,
        rmsdCut = 1.0,
        qbits_rmsd = 1.5,
        qbits_window = 10,
        secstruct = 'H',
        min_nbrs = 1,     
        lowest_rmsd_loop = True,
        ###Database
        database='/mnt/e/GitHub_Design/Qbits/database',
        loop_target_list='/mnt/e/GitHub_Design/master_db/list',   
        ###For loop searching     
        master_query_loop_top = 200,
        max_nc_dist = 15.0,
        loop_query_win =7,   
        min_loop_length = 3,
        max_loop_length=20,
        cluster_count_cut=20,
        loop_distance_cut=15,
        construct_keep = 1
)

seed_pdb = '/mnt/e/DesignData/smallprot/lee/seed_ab.pdb'
query_pdb = None
exclusion_pdb = None

workdir = '/mnt/e/DesignData/smallprot/lee/output_extend'

hhh = extend_sse.ExtendSSE(seed_pdb, query_pdb, exclusion_pdb,  workdir, para)

extension_lengths = [(12, 0), (0, 12)]
hhh.get_extention(extension_lengths)
