import os

workdir = '/mnt/e/DesignData/smallprot_loops/huong/'
#outdir = 'loop_ss_20220401-100308/'
outdir = None
target_file = 'seed.pdb'


loop_topo_sels = {
    # Toplogy
    'ABC_frt' : {
        'A-B': [('A', 21, 3, 'B', 4, 3), ('A', 18, 3, 'B', 7, 3)], # Chain A loop B. [(A, pos, +-len, B, pos, +-len)]
        'B-C': [('B', 21, 3, 'C', 6, 3), ('B', 18, 3, 'C', 9, 3)],
        'C-D': [('C', 21, 3, 'D', 5, 3), ('C', 18, 3, 'D', 8, 3)],
        'D-E': [('D', 22, 3, 'E', 10, 3)],
        'E-F': [('E', 21, 3, 'F', 4, 3), ('E', 18, 3, 'F', 7, 3)],
        'F-A': [('F', 21, 3, 'A', 4, 3), ('F', 18, 3, 'A', 7, 3)],
    },
    'ABC_rev' : {
        'A-F': [('A', 21, 3, 'F', 5, 3), ('A', 18, 3, 'F', 8, 3)],
        'F-E': [('F', 21, 3, 'E', 4, 3), ('F', 18, 3, 'E', 7, 3)],
        'E-D': [('E', 21, 3, 'D', 10, 3)],
        'D-C': [('D', 21, 3, 'C', 6, 3), ('D', 18, 3, 'C', 9, 3)],
        'C-B': [('C', 21, 3, 'B', 6, 3), ('C', 18, 3, 'B', 9, 3)],
        'B-A': [('B', 21, 3, 'A', 4, 3), ('B', 18, 3, 'A', 7, 3)],
    },
}

# loop_topo_sels = {
#     # Toplogy
#     'ABC_frt' : {
#         'A-B': [('A', 21, 3, 'B', 4, 3), ('A', 18, 3, 'B', 7, 3)], # Chain A loop B. [(A, pos, +-len, B, pos, +-len)]
#     },
# }

### Master search para
loop_range = [10, 20]
loop_target_list = '/mnt/e/GitHub_Design/master_db/list'
rmsdCut = 0.6
master_query_loop_top = 500
cluster_rmsd = 1.0

### Output para
cluster_count_cut = 10
select_min_rmsd_pdb= False
