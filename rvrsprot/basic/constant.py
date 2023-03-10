import numpy as np
from pathlib import Path

# hydrophobicity from https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html

hydro_scale = ['Kyte-Doolittle scale', 'ww', 'hh', 'mf', 'tt']

hydro_dict = {'LYS': [-3.9, -0.99, 2.71, 5.39, -3.46],
          'ASP': [-3.5, -1.23, 3.49, 2.95, -3.27],
          'PHE': [2.8, 1.13, -0.32, -2.20, 1.98],
          'ASN': [-3.5, -0.42, 2.05, 3.47, -1.62],
          'GLN': [-3.5, -0.58, 2.36, 3.01, -1.84],
          'ALA': [1.8, -0.17, 0.11, 0.0, 0.38],
          'ARG': [-4.5, -0.81, 2.58, 3.71, -2.57],
          'THR': [-0.7, 0.14, 0.52, 1.78, -0.32],
          'GLY': [-0.4, -0.01, 0.74, 1.72, -0.19],
          'TYR': [-1.3, 0.94, 0.68, -1.09, 0.49],
          'LEU': [3.8, 0.56, -0.55, -1.81, 1.82],
          'VAL': [4.2, -0.07, -0.31, -0.78, 1.46],
          'GLU': [-3.5, -2.02, 2.68, 1.64, -2.90],
          'PRO': [-1.6, -0.45, 2.23, -1.52, -1.44],
          'SER': [-0.8, -0.13, 0.84, 1.83, -0.53],
          'CYS': [2.5, 0.24, -0.13, 0.49, -0.30],
          'MET': [1.9, 0.23, -0.10, -0.76, 1.40],
          'TRP': [-0.9, 1.85, 0.30, -0.38, 1.53],
          'ILE': [4.5, 0.31, -0.60, -1.56, 1.97],
          'HIS': [-3.2, -0.96, 2.06, 4.76, -1.44],
          'MSE': [1.9, 0.23, -0.10, -0.76, 1.40],
          'HSE': [-3.2, -0.96, 2.06, 4.76, -1.44],
          }


# Based on paper "Intrinsic Secondary Structure Propensities of the Amino Acids, Using Statistical Matrices: Comparison With Experimental Scales"
propensity_scale = ['alpha_helix', 'ww', 'hh', 'mf', 'tt']

propensity_dict = {'LYS': [0.615, 0, 0, 0, 0],
          'ASP': [0.870, 0, 0, 0, 0],
          'PHE': [0.706, 0, 0, 0, 0],
          'ASN': [0.906, 0, 0, 0, 0],
          'GLN': [0.594, 0, 0, 0, 0],
          'ALA': [0.432, 0, 0, 0, 0],
          'ARG': [0.503, 0, 0, 0, 0],
          'THR': [0.884, 0, 0, 0, 0],
          'GLY': [1.162, 0, 0, 0, 0],
          'TYR': [0.778, 0, 0, 0, 0],
          'LEU': [0.494, 0, 0, 0, 0],
          'VAL': [0.706, 0, 0, 0, 0],
          'GLU': [0.167, 0, 0, 0, 0],
          'PRO': [1.945, 0, 0, 0, 0],
          'SER': [0.928, 0, 0, 0, 0],
          'CYS': [0.877, 0, 0, 0, 0],
          'MET': [0.444, 0, 0, 0, 0],
          'TRP': [0.690, 0, 0, 0, 0],
          'ILE': [0.566, 0, 0, 0, 0],
          'HIS': [0.802, 0, 0, 0, 0],
          'MSE': [0.444, 0, 0, 0, 0],
          'HSE': [0.802, 0, 0, 0, 0],
          }


one_letter_code = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                   'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                   'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                   'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
                   'MSE': 'm', 'ANY': '.', 'FE': 'fe', 'ZN': 'zn', 'HEM': 'h', 'SEP':'s', 'TPO':'t'}

inv_one_letter_code = {
    'C': 'CYS',
    'D': 'ASP',
    'S': 'SER',
    'Q': 'GLN',
    'K': 'LYS',
    'I': 'ILE',
    'P': 'PRO',
    'T': 'THR',
    'F': 'PHE',
    'N': 'ASN',
    'G': 'GLY',
    'H': 'HIS',
    'L': 'LEU',
    'R': 'ARG',
    'W': 'TRP',
    'A': 'ALA',
    'V': 'VAL',
    'E': 'GLU',
    'Y': 'TYR',
    'M': 'MET',
    'm': 'MSE',
    '.': 'ANY',
    'fe': 'FE',
    'zn': 'ZN',
    'h': 'HEM',
    's': 'SEP',
    't': 'TPO'}

resnames_aa_20 = ['CYS', 'ASP', 'SER', 'GLN', 'LYS',
                   'ILE', 'PRO', 'THR', 'PHE', 'ASN',
                   'GLY', 'HIS', 'LEU', 'ARG', 'TRP',
                   'ALA', 'VAL', 'GLU', 'TYR', 'MET',
                   'MSE', 'CSO', 'TPO', 'SEP', 'TYS', 'HIP', 'NEP', 'PTR', 'SEC']

def read_apble(filepath):
    apble_dict ={}
    with open(filepath) as file_in:
        lines = file_in.readlines()
        count = -1
        table = []
        key = ''
        for line in lines:
            if line.split('\t')[0] in resnames_aa_20:
                key = line.split('\t')[0]
                count = 0
                continue
            if count >= 0 and count <= 35:
                count += 1
                table.append([s for s in line.split('\t')[0:36]])
                if count == 35:
                    apble_dict[key] = np.array(table, dtype = object) 
            else:
                table = []
                key = ''
                count = -1
    return apble_dict



APBEL_DICT = read_apble(Path(__file__).parent.parent / 'constants/APBLE.txt')
            
                