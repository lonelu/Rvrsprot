# Smallprot for loop building
## The program has been tested on Helix bundles.

### Example of usage 
Please check files in './smallprot/gpu_run'
You will find example run files such as 'run_loop.py'

Copy the 'run_loop.py' into your data folder or where you want.
Change parameters and filepaths of the 'run_loop.py' file.
You can put the program in your local computer, but I found the program runs much faster at gpu.

### About the parameters. 
Here is all the parameters.
para = smallprot_config.Parameter(
    ###Database
        database='/home/gpu/Nick/Qbits/database',
        loop_target_list='/home/gpu/Nick/master_db/list',   
    ###Loop searching     
        master_query_loop_top = 200,
        max_nc_dist = 16.0,
        loop_query_win =7,   
        min_loop_length = 3,
        max_loop_length=20,
        select_min_rmsd_pdb = True,
    ###Loop searching -> cutoff 
        cluster_count_cut=20,
        loop_distance_cut=16,
        construct_keep = 0
)

database: The database used for master and qbits search. You don't want to change this in general. 
loop_target_list: The loop database used for master and qbits search. You don't want to change this in general. 

master_query_loop_top: The top queries used from master search. You want to use a larger number if the hits are low.
max_nc_dist: If two helix has a distance larger than 16, they are not supposed to be loop-able. 
loop_query_win: the size of the query to be selected from the helix.
min_loop_length - max_loop_length: length of loops for the master query.
select_min_rmsd_pdb: if True, the program will build complete structure with the min rmsd loop candidates. Otherwise, with centroid loop candidates.

cluster_count_cut: if a loop candidate with clustered number smaller than 20, they won't be selected for complete structure building.
loop_distance_cut: distance limitation for the loops at the same side. For example, for a 4 helix bundle looped in order A-B-C-D, loop_AB and loop_CD is at the same side and their distance should be less than 16. 
construct_keep: keep == 0, keep min; keep == -1, keep loop sequence as much as possible; keep == 1, keep seed sequence as much as possible. 
### Start to run the script.
Activite the python conda environment.
> conda activate env_smallprot
Run the script. It generally take ~10-30 mins depend on your setting.
> python run_loop.py


