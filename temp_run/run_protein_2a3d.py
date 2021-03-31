## activate the conda environment. conda activate env_smallprot
## please check/change the parameters in 'parameter_loop_truc.ini' 
## you can use ipython or just 'python run.py'
import sys
sys.path.append(r'/mnt/e/GitHub_Design/Qbits')
sys.path.append(r'/mnt/e/GitHub_Design/smallprot')

print('Thanks for using smallprot!')

from smallprot import smallprot, smallprot_config

para = smallprot_config.Parameter(
        num_iter = 2, 
        top = 20,      
        master_query_top = None,
        screen_compactness = False,
        rmsdCut = 1.0,
        qbits_rmsd = 1.5,
        qbits_window = 10,
        secstruct = 'H',
        min_nbrs = 1,     
        lowest_rmsd_loop = True,
        ###Database
        database='/home/gpu/Nick/Qbits/database',
        loop_target_list='/home/gpu/Nick/master_db/list', 
        ###For loop searching     
        master_query_loop_top = 200,
        max_nc_dist = 16.0,
        loop_query_win =7,   
        min_loop_length = 3,
        max_loop_length=15,
        cluster_count_cut=20,
        loop_distance_cut=15,
        construct_keep = 0
)

seed_pdb = '/mnt/e/DesignData/smallprot/2a3d/seed1.pdb'
query_pdb = None
exclusion_pdb = None

workdir = '/mnt/e/DesignData/smallprot/2a3d/output_build/'

hhh = smallprot.SmallProt(seed_pdb, query_pdb, exclusion_pdb,  workdir, para)
#n_truncations=list(range(5))
#c_truncations=list(range(5))
#hhh.build_protein(n_truncations = n_truncations, c_truncations = c_truncations)
hhh.build_protein()