## activate the conda environment. conda activate env_smallprot
## please check/change the parameters in 'parameter_loop_truc.ini' 
## you can use ipython or just 'python run.py'
import sys
sys.path.append(r'/mnt/e/GitHub_Design/Qbits')

print('Thanks for using smallprot!')

from smallprot import smallprot 

seed_pdb = '/mnt/e/GitHub_Design/smallprot/data/nina/seed.pdb'
query_pdb = None
exclusion_pdb = None

workdir = '/mnt/e/GitHub_Design/smallprot/data/nina/output_test_build1/'
para = '/mnt/e/GitHub_Design/smallprot/parameter_loop_truc.ini'

hhh = smallprot.SmallProt(seed_pdb, query_pdb, exclusion_pdb,  workdir, para)
# n_truncations=list(range(20))
# c_truncations=list(range(10))
# n_truncations = [16, 17, 18, 19]
# c_truncations = [1, 2, 3, 4, 5, 6, 7, 8, 9]
n_truncations = [16]
c_truncations = [1]
direction=[2, 1, 0, 3]
# direction = None
hhh.loop_structure(direction = direction, n_truncations = n_truncations, c_truncations = c_truncations)