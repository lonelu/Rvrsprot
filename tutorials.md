## please check/change the parameters in 'parameter_loop_truc.ini' 

## change directory to the smallprot folder
## cd ~/Lei/smallprot

from smallprot import smallprot 
hhh = smallprot.SmallProt('parameter_loop_truc.ini')
# n_truncations=list(range(20))
# c_truncations=list(range(10))
# n_truncations = [16, 17, 18, 19]
# c_truncations = [1, 2, 3, 4, 5, 6, 7, 8, 9]
n_truncations = [16]
c_truncations = [1]
hhh.loop_seed_single_structure(direction=[2, 1, 0, 3], n_truncations = n_truncations, c_truncations = c_truncations)
# hhh.loop_seed_single_structure( n_truncations = n_truncations, c_truncations = c_truncations)