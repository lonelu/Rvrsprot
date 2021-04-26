import logomaker
import matplotlib.pyplot as plt
from fpdf import FPDF
import pandas as pd
import numpy as np
import qbits
from smallprot import constant
from smallprot import extract_master
import os


def _plot_log(fig, ax, seqs):
    seqs_one_letter = extract_master._get_loop_candidate_seqs_one_letter(seqs)
    df = logomaker.alignment_to_matrix(sequences=seqs_one_letter, to_type='counts',
                                            characters_to_ignore='-', pseudocount=0.01)
    logo = logomaker.Logo(df,
                        ax = ax,
                        #font_name='DejaVu Sans',
                        color_scheme='NajafabadiEtAl2017',
                        vpad=.1,
                        width=.8)
  
    ax.set_ylabel('Count')
    xs = list(range(1, len(df)+1))
    ax.set_xticks(range(len(xs)))
    ax.set_xticklabels(xs)


def _plot_hydro(fig, ax, seqs, seqlen, loop_query_win):
    all_hydro=[]
    for seq in seqs:
        hydro = [constant.hydro_dict[res] for res in seq]
        all_hydro.append(hydro)

    means = np.mean(all_hydro, 0)
    sds = np.std(all_hydro, 0)

    x = list(range(1, seqlen + 2*loop_query_win+1))
    #You can plot different hydro_scales or all of them.
    #for i in range(len(constant.hydro_scale)):
    i = 0  
    ax.set_xlabel('AA', fontsize = 12)
    ax.set_ylabel(constant.hydro_scale[i], fontsize = 12)
    ax.errorbar(x, means[:, i], yerr = sds[:, i])
    ax.set_xticks(x)


def _plot_propensity(fig, ax, seqs, seqlen, loop_query_win):

    all_propen = []
    for seq in seqs:
        propen = [constant.propensity_dict[res] for res in seq]
        all_propen.append(propen)

    means = np.mean(all_propen, 0)
    sds = np.std(all_propen, 0)

    x = list(range(1, seqlen + 2*loop_query_win+1))
    #for i in range(len(constant.propensity_scale)):
    #fig = plt.figure(figsize=(10, 8))
    #ax = fig.add_subplot(111)            
    ax.set_xlabel('AA', fontsize = 18)
    ax.set_ylabel(constant.propensity_scale[0], fontsize = 18)

    ax.errorbar(x, means[:, 0], yerr = sds[:, 0])
    plt.savefig(filepath+'_helix_propen.png')
    plt.clf()

def _plot_phipsi(fig, ax, phi, psi, seqlen):
    #plt.figure(figsize=(10, 3.5))
    x = list(range(1, len(phi) + 1))
    ax.plot(x, phi, label = 'phi')
    ax.plot(x, psi, label = 'psi')
    ax.legend()
    ax.set_xticks(x)
    ax.set_ylabel("Angle", fontsize = 12)

def _plot_table(fig, ax, seq, phi, psi):
    data = []
    data.append([qbits.constants.one_letter_code[s] for s in seq])
    data.append([round(f, 1) if f else None for f in phi])
    data.append([round(f, 1) if f else None for f in psi])
    character = []
    for i in range(len(seq)):
        s = seq[i]
        if phi[i] and psi[i]:  
            phi_ind = int((phi[i]+ 180)/10 -1)    
            psi_ind = int(35 - (psi[i]+ 180)/10)  
            if s not in constant.APBEL_DICT.keys():
                s = 'ALA'     
            character.append(constant.APBEL_DICT[s][psi_ind, phi_ind])
        else:
            character.append(None)
    data.append(character)

    ax.set_axis_off()
    df=pd.DataFrame(data)
    ax.axis('tight')
    ax.axis('off')
    rows = ['AA', 'phi', 'psi', 'CC']
    tab = ax.table(cellText=df.values, colLabels= list(range(1, len(seq)+1)), rowLabels=rows, cellLoc='center', loc='upper left')
    fig.tight_layout()
    tab.set_fontsize(20)
    tab.scale(1, 2)
    tab.auto_set_font_size
    ax.set_ylabel("Legend", fontsize = 12)

def _plot_rmsds(fig, ax, loop_rmsds):
    ax.hist(loop_rmsds, bins = 10, color = 'blue', edgecolor = 'black')  
    # Title and labels
    ax.set_xlabel('RMSD', size = 12)
    ax.set_ylabel('Count', size= 12)


def _plot_all(filepath, loop_seqs, loop_rmsds, seqlen, loop_query_win, phi, psi, seq):
    fig, (ax1, ax2, ax3, ax4, ax5) =plt.subplots(5, 1, figsize=(15, 17.5))
    _plot_rmsds(fig, ax1, loop_rmsds)
    _plot_phipsi(fig, ax2, phi, psi, seqlen)
    _plot_table(fig, ax3, seq, phi, psi)
    _plot_log(fig, ax4, loop_seqs)
    _plot_hydro(fig, ax5, loop_seqs, seqlen, loop_query_win)
    plt.tight_layout()
    plt.savefig(filepath+'_info.png')
    plt.close()

    ### How to generate pdf file.
    # pdf = FPDF('P', 'mm', (210, 320))
    # pdf.set_margins(left= -0.5, top = 0)
    # pdf.add_page()
    # pdf.set_xy(0, 0)
    # pdf.set_font('Arial', 'B', 12)
    # pdf.cell(30)
    # pdf.cell(80, 10, filepath.split('/')[-1], 0, 2, 'C')
    # pdf.cell(-30)
    # pdf.cell(-80, 0, " ", 0, 2, 'C')
    # #pdf.image(filepath + '_ramachandran.png', x = None, y = None, w = 0, h = 0, type = '', link = '')
    # pdf.image(filepath + '_phipsi.png', x = None, y = None, w = 200, h = 70, type = '', link = '')
    # pdf.image(filepath + '_seqTable.png', x = None, y = None, w = 200, h = 70, type = '', link = '')
    # pdf.image(filepath + '_logo.png', x = None, y = None, w = 200, h = 70, type = '', link = '')
    # pdf.image(filepath + '_hydro.png', x = None, y = None, w = 200, h = 70, type = '', link = '')  
    # pdf.output(filepath + '_report.pdf', 'F')
    # #os.remove(filepath + '_.png')

def _plot_contact_map(filepath):
    pdb_prody = pr.parsePDB(file_path)
    ca = pdb_prody.select('name CA')
    pr.showContactMap(ca, cmap = 'Reds')
    plt.savefig(filepath[:-4] + '.png')
    plt.close()

