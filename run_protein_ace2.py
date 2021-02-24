## activate the conda environment. conda activate env_smallprot
## please check/change the parameters in 'parameter_loop_truc.ini' 
## you can use ipython or just 'python run.py'

print('Thanks for using smallprot!')

from smallprot import smallprot 

seed_pdb = '/mnt/e/GitHub_Design/smallprot/data/ace2_input/query.pdb'
query_pdb = None
exclusion_pdb = '/mnt/e/GitHub_Design/smallprot/data/ace2_input/exclusion.pdb'

workdir = '/mnt/e/GitHub_Design/smallprot/data/ace2_input/output_build_3helix_test/'
para = '/mnt/e/GitHub_Design/smallprot/parameter_loop_truc.ini'

hhh = smallprot.SmallProt(seed_pdb, query_pdb, exclusion_pdb,  workdir, para)
#n_truncations=list(range(5))
#c_truncations=list(range(5))
#hhh.build_protein(n_truncations = n_truncations, c_truncations = c_truncations)
hhh.build_protein()