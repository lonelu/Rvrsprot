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
    def __init__(self):
        self.query_pdb = None
        self.seed_pdb = None
        self.exclusion_pdb = None
        self.workdir = None
        self.num_iter = 3
        self.max_nc_dist = 15.0
        self.screen_compactness = False
        self.rmsdCut = 1.0
        self.qbits_rmsd = 1.5
        self.qbits_window = 10
        self.secstruct = None
        self.min_nbrs = 1
        self.min_loop_length = 3
        self.max_loop_length=20
        self.lowest_rmsd_loop = True
        self.database='/mnt/e/GitHub_Design/Qbits/database'
        self.loop_target_list='/mnt/e/GitHub_Design/master_db/list'

#Write and Read Parameters used for Smallprot

def writeConfig(folderPath = ''):
    config = configparser.ConfigParser()
    config['Smallprot'] = {'query_pdb': 'None',
                        'seed_pdb': 'None',
                        'exclusion_pdb': 'None',
                        'workdir': 'None',
                        'num_iter': '3',
                        'max_nc_dist': '15',
                        'screen_compactness': 'false',
                        'rmsdCut': '1.0',
                        'qbits_rmsd': '1.5',
                        'qbits_window': '10',
                        'secstruct': 'None',
                        'min_nbrs': '1',
                        'min_loop_length': '3',
                        'max_loop_length': '20',
                        'lowest_rmsd_loop': 'true',
                        'database': '/mnt/e/GitHub_Design/Qbits/database',
                        'loop_target_list': '/mnt/e/GitHub_Design/master_db/list'}

    with open(folderPath + 'parameter.ini', 'w') as configfile:
        config.write(configfile)


def readConfig(filePath = 'parameter.ini'):
    cfg = configparser.ConfigParser()
    if not os.path.exists(filePath):
        print('Config file is not existed in the current folder.')
    cfg.read(filePath)
    para = Parameter()
    para.query_pdb = cfg['Smallprot']['query_pdb'] if cfg['Smallprot']['query_pdb']!='None' else None
    para.seed_pdb = cfg['Smallprot']['seed_pdb'] if cfg['Smallprot']['seed_pdb']!='None' else None
    para.exclusion_pdb = cfg['Smallprot']['exclusion_pdb'] if cfg['Smallprot']['exclusion_pdb']!='None' else None
    para.workdir = cfg['Smallprot']['workdir'] if cfg['Smallprot']['workdir']!='None' else None
    para.num_iter = cfg.getint('Smallprot','num_iter')
    para.max_nc_dist = cfg.getfloat('Smallprot','max_nc_dist')
    para.screen_compactness = cfg.getboolean('Smallprot','screen_compactness')
    para.rmsdCut = cfg.getfloat('Smallprot','rmsdCut')
    para.qbits_rmsd = cfg.getfloat('Smallprot','qbits_rmsd')
    para.qbits_window = cfg.getint('Smallprot','qbits_window')
    para.workdir = cfg['Smallprot']['workdir'] if cfg['Smallprot']['workdir']!='None' else None
    para.min_nbrs = cfg.getint('Smallprot','min_nbrs')
    para.min_loop_length = cfg.getint('Smallprot','min_loop_length')
    para.max_loop_length = cfg.getint('Smallprot','max_loop_length')
    para.lowest_rmsd_loop = cfg.getboolean('Smallprot','lowest_rmsd_loop')
    para.database = cfg['Smallprot']['database'] if cfg['Smallprot']['database']!='None' else None
    para.loop_target_list = cfg['Smallprot']['loop_target_list'] if cfg['Smallprot']['loop_target_list']!='None' else None
    return para










