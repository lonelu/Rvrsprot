## activate the conda environment. conda activate env_smallprot
## please check/change the parameters in 'parameter_loop_truc.ini' 
## you can use ipython or just 'python run_rocker.py'

print('Thanks for using smallprot!')

from smallprot import smallprot 

query_pdb = None
seed_pdb = '/mnt/e/GitHub_Design/smallprot/data/rocker/seed_correct.pdb'
exclusion_pdb = None
workdir = '/mnt/e/GitHub_Design/smallprot/data/rocker/output_build'
para = '/mnt/e/GitHub_Design/smallprot/parameter_loop_truc_rocker.ini'

hhh = smallprot.SmallProt(query_pdb, seed_pdb, exclusion_pdb, workdir, para)
# n_truncations=list(range(20))
# c_truncations=list(range(10))
n_truncations = [1, 2, 3]
c_truncations = [1, 2]
hhh.loop_seed_single_structure(n_truncations = n_truncations, c_truncations = c_truncations)