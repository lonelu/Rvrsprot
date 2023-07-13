# Generate hid ASP helix.

import os
from rvrsprot.combs import reverse_vdm

path_to_database='/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/vdMs/'
helix_path =  '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/15mer_ALA.pdb'

helix_cg_path =  '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/AAA_cluster_0_sc_0.64.pdb'
helix_cg_sel = 'name CG ND1 CD2 CE1 NE2'
vdm_cg_sel = 'chid Y and name CG ND1 CD2 CE1 NE2'

cg = 'hid'
aa = 'GLU'
helix_resind = 7
outdir_vh = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/hid_GLU_H/'
os.makedirs(outdir_vh, exist_ok = True)
reverse_vdm.generate_vdm2H(outdir_vh, path_to_database, cg, aa, helix_path, helix_resind, helix_cg_path, helix_cg_sel, vdm_cg_sel)

cg = 'hid'
aa = 'ASP'
helix_resind = 7
outdir_vh = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/_rvrs_vdm/hid_ASP_H/'
os.makedirs(outdir_vh, exist_ok = True)
reverse_vdm.generate_vdm2H(outdir_vh, path_to_database, cg, aa, helix_path, helix_resind, helix_cg_path, helix_cg_sel, vdm_cg_sel)

