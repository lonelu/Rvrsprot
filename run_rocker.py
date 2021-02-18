## activate the conda environment. conda activate env_smallprot
## please check/change the parameters in 'parameter_loop_truc.ini' 
## you can use ipython or just 'python run_rocker.py'

print('Thanks for using smallprot!')

from smallprot import smallprot 
hhh = smallprot.SmallProt('/mnt/e/GitHub_Design/smallprot/parameter_loop_truc_rocker.ini')
# n_truncations=list(range(20))
# c_truncations=list(range(10))
n_truncations = [1]
c_truncations = [1]
hhh.loop_seed_single_structure(n_truncations = n_truncations, c_truncations = c_truncations)