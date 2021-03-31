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

seed_pdb = '/mnt/e/DesignData/smallprot/trihelix/seed2.pdb'
query_pdb = None
exclusion_pdb = None

workdir = '/mnt/e/DesignData/smallprot/trihelix/output_top20'

hhh = smallprot.SmallProt(seed_pdb, query_pdb, exclusion_pdb,  workdir, para)

hhh.build_protein()

'''
hhh.pre_build_pdbs_summary[0]

def _write_protein_summary(filename, pdb_paths, infos):
    with open(filename, 'w') as f:
        f.write('pdb_path\tca_in_ahull_ratio\tmedian_min_dist_diff\n')
        for i in range(len(pdb_paths)):
            median_min_dist = '\t'.join([str(d) for d in infos[i].median_min_dists])
            f.write(pdb_paths[i] + '\t' + str(infos[i].ahull_in_ratio) +  '\t' + median_min_dist + '\n')

filename = hhh.workdir + '/protein_build_summary.txt' 
hhh._write_protein_summary(filename, hhh.pre_build_pdbs, hhh.pre_build_pdbs_summary)


'''