import os
import sys
from rvrsprot.vdm_rev import reverse_vdm


def vdm_db_key(path_to_database):
    cg_aa_lists = []
    for folder in os.listdir(path_to_database):
        if not os.path.isdir(path_to_database + folder):
            continue
        for file in os.listdir(path_to_database + folder):
            if not os.path.isfile(path_to_database + folder + '/' + file):
                continue
            cg_aa_lists.append((folder, file))
    return cg_aa_lists


def generate_rvdm():
    ind = sys.argv[1]-1
    #path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'
    path_to_database='/wynton/home/degradolab/lonelu/GitHub_Design/Combs2_library/vdMs/'
    cg_aa_lists = vdm_db_key(path_to_database)

    cg = cg_aa_lists[ind][0] 
    aa = cg_aa_lists[ind][1] 

    #outdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/vdMs/' + cg + '/'
    outdir = '/wynton/home/degradolab/lonelu/DesignData/Database/rvdMs/' + cg + '/'
    os.makedirs(outdir, exist_ok=True)
    reverse_vdm.generate_rev_vdm_lib(path_to_database, outdir, cg, aa)

    return 

if __name__=='__main__':
    generate_rvdm()

# Generate coo HIS.
'''
path_to_database='/mnt/e/DesignData/Combs/Combs2_database/'

cg = 'coo'
aa = 'HIS'

outdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/vdMs/coo/'

generate_rev_vdm_lib(path_to_database, outdir, cg, aa)

outdir_pdb = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/coo_his_pdb_100/'
write_vdm('/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/vdMs/', cg, aa, outdir_pdb, total = 100)

outdir_pdb_ori = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/coo_his_pdb_100_origin/'
write_vdm(path_to_database, cg, aa, outdir_pdb_ori, total = 100)
'''

# Generate hid ASP.
'''
path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'

cg = 'hid'
aa = 'ASP'

import os
outdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/vdMs/hid/'
os.makedirs(outdir, exist_ok = True)

generate_rev_vdm_lib(path_to_database, outdir, cg, aa)

outdir_pdb = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/hid_ASP_pdb_100/'
os.makedirs(outdir_pdb, exist_ok = True)
write_vdm('/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/vdMs/', cg, aa, outdir_pdb, total = 100)


outdir_pdb_ori = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/hid_ASP_pdb_100_origin/'
os.makedirs(outdir_pdb_ori, exist_ok = True)
write_vdm(path_to_database, cg, aa, outdir_pdb_ori, total = 100)

'''

# Generate hid GLU.
'''
path_to_database='/mnt/e/DesignData/Combs/Combs2_database/vdMs/'

cg = 'hid'
aa = 'GLU'

import os
outdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/vdMs/hid/'
os.makedirs(outdir, exist_ok = True)

generate_rev_vdm_lib(path_to_database, outdir, cg, aa)

outdir_pdb = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/hid_GLU_pdb_100/'
os.makedirs(outdir_pdb, exist_ok = True)
write_vdm('/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/vdMs/', cg, aa, outdir_pdb, total = 100)


outdir_pdb_ori = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/hid_GLU_pdb_100_origin/'
os.makedirs(outdir_pdb_ori, exist_ok = True)
write_vdm(path_to_database, cg, aa, outdir_pdb_ori, total = 100)

'''

