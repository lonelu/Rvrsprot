import os
import configparser

""" Parameters for Smallprot
    query_pdb : str
        Path to PDB file of query structure from which to generate the 
        initial Qbit rep in the first design iteration.
    seed_pdb : str, optional
        Path to PDB file of seed structure to include both in the query 
        and in the final output protein structures.
    exclusion_pdb : str, optional
        Path to PDB file of structure to be used to define a region of 
        volume exclusion in all Qbits searches.
    workdir : str, optional
        Path at which to create a working directory where the MASTER 
        query PDB file and output file structure will be created. (If 
        None, they will be created in the directory of seed_pdb.)
    num_iter : int, optional
        Number of iterations of MASTER and Qbits to be used in finding 
        secondary structure to add to the designed protein.
    max_nc_dist : float, optional
        Maximum distance between N and C termini of adjacent secondary 
        structural elements of the designed protein.
    screen_compactness : bool, optional
        If True, ensure that all output structures satisfy the alpha hull-
        based compactness criterion in the pdbutils submodule.
    rmsdCut : float, optional
        rmsdCut for MASTER queries.
    qbits_rmsd : float, optional
        RMSD threshold for qbits searches.
    qbits_window : int, optional
        Window size for qbits searches.
    secstruct : str, optional
        DSSP code for allowable secondary structure in the designed 
        protein.  If None, all secondary structure is allowable.
    min_nbrs : int, optional
        Minimum number of neighbors for a residue in a qbit rep.
    min_loop_length : int, optional
        Minimum length of loops in the final structure. Default: 3.
    max_loop_length : int, optional
        Maximum length of loops in the final structure. Default: 20.
    lowest_rmsd_loop : bool, optional
        If True, extract the lowest-RMSD loop to the query from each 
        cluster instead of the cluster centroid.
    database : str, optional
        Path to database folder of pickled Prody objects for use in 
        MASTER queries and Qbits searches.
    target_list : str, optional
        Filename of target list within the database of objects to 
        be used in MASTER queries for SSEs.
    loop_target_list : str, optional
        Filename of target list within the database of objects to 
        be used in MASTER queries for loops.
"""
class Parameter:
    def __init__(self, num_iter = 3, top = 5, master_query_top = 200, screen_compactness = False, rmsdCut = 1.0, 
    qbits_rmsd = 1.5, qbits_window = 10, secstruct = None, min_nbrs = 1, lowest_rmsd_loop = True, 
    database='/mnt/e/GitHub_Design/Qbits/database', loop_target_list='/mnt/e/GitHub_Design/master_db/list', 
    master_query_loop_top = 200, max_nc_dist = 15.0, loop_query_win =7, min_loop_length = 3, max_loop_length=20, 
    cluster_count_cut=20, loop_distance_cut=15, construct_keep = 0):
        self.num_iter = num_iter  
        self.top = top      
        self.master_query_top = master_query_top
        self.screen_compactness = screen_compactness
        self.rmsdCut = rmsdCut
        self.qbits_rmsd = qbits_rmsd
        self.qbits_window = qbits_window
        self.secstruct = secstruct
        self.min_nbrs = min_nbrs       
        self.lowest_rmsd_loop = lowest_rmsd_loop
        ###Database
        self.database=database
        self.loop_target_list=loop_target_list   
        ###For loop searching     
        self.master_query_loop_top = master_query_loop_top
        self.max_nc_dist = max_nc_dist
        self.loop_query_win =loop_query_win  
        self.min_loop_length = min_loop_length
        self.max_loop_length=max_loop_length
        self.cluster_count_cut=cluster_count_cut
        self.loop_distance_cut=loop_distance_cut
        self.construct_keep = construct_keep
#Write and Read Parameters used for Smallprot

def writeConfig(folderPath = ''):
    config = configparser.ConfigParser()
    config['Smallprot'] = {
                        'num_iter': '3',   
                        'top': '5', 
                        'master_query_top': 'master_query_top',                  
                        'screen_compactness': 'false',
                        'rmsdCut': '1.0',
                        'qbits_rmsd': '1.5',
                        'qbits_window': '10',
                        'secstruct': 'None',
                        'min_nbrs': '1',
                        'lowest_rmsd_loop': 'true',
                        'database': '/mnt/e/GitHub_Design/Qbits/database',
                        'loop_target_list': '/mnt/e/GitHub_Design/master_db/list',  
                        'master_query_loop_top': 'master_query_loop_top',
                        'max_nc_dist': '15',
                        'loop_query_win': '7',                      
                        'min_loop_length': '3',
                        'max_loop_length': '20',
                        'cluster_count_cut': '20',
                        'loop_distance_cut': '15',
                        'construct_keep': '0'
                        }

    with open(folderPath + 'parameter.ini', 'w') as configfile:
        config.write(configfile)


def writeConfig(filePath, para):
    config = configparser.ConfigParser()
    config['Smallprot'] = {
                        'num_iter': str(para.num_iter),   
                        'top':str(para.top),
                        'master_query_top':str(para.master_query_top),              
                        'screen_compactness': str(para.screen_compactness),
                        'rmsdCut': str(para.rmsdCut),
                        'qbits_rmsd': str(para.qbits_rmsd),
                        'qbits_window': str(para.qbits_window),
                        'secstruct': str(para.secstruct),
                        'min_nbrs': str(para.min_nbrs),
                        'lowest_rmsd_loop': str(para.lowest_rmsd_loop),
                        'database': str(para.database),
                        'loop_target_list': str(para.loop_target_list), 
                        'master_query_top':str(para.master_query_loop_top), 
                        'max_nc_dist': str(para.max_nc_dist),
                        'loop_query_win': str(para.loop_query_win),                      
                        'min_loop_length': str(para.min_loop_length),
                        'max_loop_length': str(para.max_loop_length),
                        'cluster_count_cut': str(para.cluster_count_cut),
                        'loop_distance_cut': str(para.loop_distance_cut),
                        'construct_keep': str(para.construct_keep)
                        }

    with open(filePath, 'w') as configfile:
        config.write(configfile)


def readConfig(filePath = 'parameter.ini'):
    cfg = configparser.ConfigParser()
    if not os.path.exists(filePath):
        print('Config file is not existed in the current folder.')
    cfg.read(filePath)
    para = Parameter()
    para.num_iter = cfg.getint('Smallprot','num_iter')    
    para.num_iter = cfg.getint('Smallprot','top') 
    para.master_query_top = cfg.getint('Smallprot', 'master_query_top')
    para.screen_compactness = cfg.getboolean('Smallprot','screen_compactness')
    para.rmsdCut = cfg.getfloat('Smallprot','rmsdCut')
    para.qbits_rmsd = cfg.getfloat('Smallprot','qbits_rmsd')
    para.qbits_window = cfg.getint('Smallprot','qbits_window')
    para.min_nbrs = cfg.getint('Smallprot','min_nbrs')
    para.lowest_rmsd_loop = cfg.getboolean('Smallprot','lowest_rmsd_loop')
    para.database = cfg['Smallprot']['database'] if cfg['Smallprot']['database']!='None' else None
    para.loop_target_list = cfg['Smallprot']['loop_target_list'] if cfg['Smallprot']['loop_target_list']!='None' else None 
    para.master_query_loop_top = cfg.getint('Smallprot', 'master_query_loop_top')
    para.max_nc_dist = cfg.getfloat('Smallprot','max_nc_dist')
    para.loop_query_win = cfg.getint('Smallprot', 'loop_query_win')  
    para.min_loop_length = cfg.getint('Smallprot','min_loop_length')
    para.max_loop_length = cfg.getint('Smallprot','max_loop_length')
    para.cluster_count_cut = cfg.getint('Smallprot','cluster_count_cut')
    para.loop_distance_cut = cfg.getint('Smallprot','loop_distance_cut')
    para.construct_keep = cfg.getint('Smallprot', 'construct_keep')
    return para










