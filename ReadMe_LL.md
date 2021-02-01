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
>> sm = smallprot.SmallProt(query_pdb="ace2_input/query.pdb", exclusion_pdb="ace2_input/exclusion.pdb", workdir="/mnt/e/GitHub_Design/smallprot/ace2_input/output/")

## build loop
>> hh = smallprot.SmallProt(seed_pdb="huong/seed.pdb", workdir="/mnt/e/GitHub_Design/smallprot/huong/output/")
>> hh.loop_seed_structure()


## Test Para
from smallprot import smallprot_config
smallprot_config.writeConfig() 

## build protein2.
from smallprot import smallprot 
sm = smallprot.SmallProt('parameter.ini')
sm.para.workdir
# sm.build_protein_deprecate()
# sm.build_protein()
sm.build_protein_parallel()