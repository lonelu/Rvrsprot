# create the environment
> cd /mnt/e/GitHub_Design/smallprot
> conda env create -f environment_linux.yml
> conda activate env_smallprot
> ipython
>> import smallprot

# However, the Qbits is not able to be imported. Add the PYTHONPATH
> cd ~
> nano .bashrc
# Add to the bottom of the .bashrc file 
 export PYTHONPATH="/mnt/e/GitHub_Design/Qbits/"
# Add path of master and qbits. (If not working, put the 'export...' lines in .bashrc)
> echo $PATH
> export PATH=$PATH:/mnt/e/GitHub_Design/master-v1.6/bin
> export PATH=$PATH:/mnt/e/GitHub_Design/Qbits


# restart the terminal.
> echo $PYTHONPATH


# Test the ace2 input
> cd smallprot
> ipython
>> from smallprot import smallprot 

## build protein. 
>> sm = smallprot.SmallProt(query_pdb="ace2_input/query.pdb", exclusion_pdb="ace2_input/exclusion.pdb", workdir="/mnt/e/GitHub_Design/zips/smallprot/ace2_input/output_test/")

## build loop
>> hh = smallprot.SmallProt(seed_pdb="loop/seed.pdb", workdir="/mnt/e/GitHub_Design/zips/smallprot/loop/output/")
>> hh.loop_seed_structure()


## Test Para
from smallprot import smallprot_config
smallprot_config.writeConfig() 


## build protein2.
from smallprot import smallprot 
sm = smallprot.SmallProt('parameter.ini')
sm.para.workdir
# sm.build_protein_deprecate()
sm.build_protein()
# sm.build_protein_parallel()

## build loop2
from smallprot import smallprot 
hh = smallprot.SmallProt('parameter_loop.ini')
# hh.loop_seed_structure()
hh.loop_seed_single_structure(direction=[2, 1, 0, 3])

## test build loop trunc
from smallprot import smallprot 
hhh = smallprot.SmallProt('parameter_loop_truc.ini')
n_truncations=[13, 14, 15, 16, 17, 18, 19]
c_truncations=[0, 1, 2, 3, 4, 5]
hhh.loop_seed_single_structure(direction=[2, 1, 0, 3], n_truncations = n_truncations, c_truncations = c_truncations)