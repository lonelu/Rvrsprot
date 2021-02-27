import logomaker
import matplotlib.pyplot as plt
from fpdf import FPDF
import pandas as pd
import numpy as np
import qbits
from smallprot import constant
import os


def _plot_log( seqfile, seqlen, filepath):
    with open(seqfile, 'r') as f:
        lines = f.read().split('\n')
    all_seqs = []
    for line in lines:
        if len(line) > 0:
            seq = ''
            for res in line.split(' '):
                if len(res) == 3 and res[0].isalpha():
                    seq += qbits.constants.one_letter_code[res]
                elif len(res) == 4 and res[0] == '[':
                    seq += qbits.constants.one_letter_code[res[1:]]
                elif len(res) == 4 and res[-1] == ']':
                    seq += qbits.constants.one_letter_code[res[:-1]]
            all_seqs.append(seq)
    seqs = []
    for s in all_seqs:
        if len(s) == 14 + seqlen:            
            seqs.append(s)

    plt.figure(figsize=(10, 3.5))
    df = logomaker.alignment_to_matrix(sequences=seqs, to_type='counts',
                                            characters_to_ignore='-', pseudocount=0.01)
    logo = logomaker.Logo(df,
                        font_name='Arial',
                        color_scheme='NajafabadiEtAl2017',
                        vpad=.1,
                        width=.8)
    #logo.style_xticks(anchor=1, spacing=1)      
    logo.ax.set_ylabel('Count')
    #logo.ax.set_xlim([0, len(df)+1])
    xs = list(range(1, len(df)+1))
    logo.ax.set_xticks(range(len(xs)))
    logo.ax.set_xticklabels(xs)
    #logo.fig.savefig(filepath) 
    plt.tight_layout() 
    plt.savefig(filepath + '_logo.png')
    plt.close()

def _plot_hydro(seqfile, seqlen, filepath, loop_query_win):
    with open(seqfile, 'r') as f:
        lines = f.read().split('\n')
    all_hydro = []
    for line in lines:
        if len(line) > 0:
            hydro = []
            for res in line.split(' '):
                if len(res) == 3 and res[0].isalpha():
                    hydro.append(constant.hydro_dict[res]) 
                elif len(res) == 4 and res[0] == '[':
                    hydro.append(constant.hydro_dict[res[1:]])
                elif len(res) == 4 and res[-1] == ']':
                    hydro.append(constant.hydro_dict[res[:-1]])
            if len(hydro) == seqlen + 2*loop_query_win:
                all_hydro.append(hydro)
    all_hydro_arr = np.array(all_hydro)

    means = np.mean(all_hydro_arr, 0)
    sds = np.std(all_hydro_arr, 0)

    x = list(range(1, seqlen + 2*loop_query_win+1))
    #You can plot different hydro_scales or all of them.
    #for i in range(len(constant.hydro_scale)):
    i = 0
    fig = plt.figure(figsize=(10, 3.5))
    ax = fig.add_subplot(111)            
    ax.set_xlabel('AA', fontsize = 12)
    ax.set_ylabel(constant.hydro_scale[i], fontsize = 12)
    ax.errorbar(x, means[:, i], yerr = sds[:, i])
    #plt.savefig(filepath+'_hydro_'+str(i)+'.png')
    plt.xticks(x)
    plt.tight_layout()
    plt.savefig(filepath+'_hydro.png')
    plt.close()

def _plot_propensity(seqfile, seqlen, filepath, loop_query_win):
    with open(seqfile, 'r') as f:
        lines = f.read().split('\n')
    all_propen = []
    for line in lines:
        if len(line) > 0:
            propen = []
            for res in line.split(' '):
                if len(res) == 3 and res[0].isalpha():
                    propen.append(constant.propensity_dict[res]) 
                elif len(res) == 4 and res[0] == '[':
                    propen.append(constant.propensity_dict[res[1:]])
                elif len(res) == 4 and res[-1] == ']':
                    propen.append(constant.propensity_dict[res[:-1]])
            if len(propen) == seqlen + 2*loop_query_win:
                all_propen.append(propen)
    all_propen_arr = np.array(all_propen)
    means = np.mean(all_propen_arr, 0)
    sds = np.std(all_propen_arr, 0)

    x = list(range(1, seqlen + 2*loop_query_win+1))
    #for i in range(len(constant.propensity_scale)):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)            
    ax.set_xlabel('AA', fontsize = 18)
    ax.set_ylabel(constant.propensity_scale[0], fontsize = 18)

    ax.errorbar(x, means[:, 0], yerr = sds[:, 0])
    plt.savefig(filepath+'_helix_propen.png')
    plt.close()

def _plot_phipsi(phi, psi, seqlen, filepath):
    plt.figure(figsize=(10, 3.5))
    x = list(range(1, len(phi) + 1))
    plt.plot(x, phi, label = 'phi')
    plt.plot(x, psi, label = 'psi')
    plt.legend()
    plt.xticks(x)
    plt.tight_layout()
    plt.ylabel("Angle", fontsize = 12)
    plt.savefig(filepath + '_phipsi.png')
    plt.close()

    # plt.plot(phipsi[2:-2:2], phipsi[3:-1:2], 's', color='red', markersize=5, markerfacecolor='white')
    # plt.xlim(-180, 180)
    # plt.ylim(-180, 180)
    # xticks = [-180, -135, -90, -45, 0, 45, 90, 135, 180]
    # yticks = [-180, -135, -90, -45, 0, 45, 90, 135, 180]
    # plt.xticks(xticks)
    # plt.yticks(yticks)
    # plt.savefig(filepath + '_ramachandran.png')
    # plt.close()

def _plot_table(seq, phi, psi, filepath):
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

    fig, ax =plt.subplots(figsize=(10, 3.5))
    ax.set_axis_off()
    df=pd.DataFrame(data)
    ax.axis('tight')
    ax.axis('off')
    rows = ['AA', 'phi', 'psi', 'CC']
    tab = ax.table(cellText=df.values, colLabels= list(range(1, len(seq)+1)), rowLabels=rows, cellLoc='center', loc='upper left')
    fig.tight_layout()
    tab.set_fontsize(18)
    tab.scale(1, 2)
    tab.auto_set_font_size
    ax.set_ylabel("Legend", fontsize = 12)
    plt.savefig(filepath+ '_seqTable.png')
    plt.close()

def _plot_all(filepath, seqfile, seqlen, loop_query_win, phi, psi, seq):
    _plot_log(seqfile, seqlen, filepath)
    _plot_hydro(seqfile, seqlen, filepath, loop_query_win)
    _plot_phipsi(phi, psi, seqlen, filepath)
    _plot_table(seq, phi, psi, filepath)

    pdf = FPDF('P', 'mm', (210, 320))
    pdf.set_margins(left= -0.5, top = 0)
    pdf.add_page()
    pdf.set_xy(0, 0)
    pdf.set_font('arial', 'B', 12)
    pdf.cell(30)
    pdf.cell(80, 10, filepath.split('/')[-1], 0, 2, 'C')
    pdf.cell(-30)
    pdf.cell(-80, 0, " ", 0, 2, 'C')
    #pdf.image(filepath + '_ramachandran.png', x = None, y = None, w = 0, h = 0, type = '', link = '')
    pdf.image(filepath + '_phipsi.png', x = None, y = None, w = 200, h = 70, type = '', link = '')
    pdf.image(filepath + '_seqTable.png', x = None, y = None, w = 200, h = 70, type = '', link = '')
    pdf.image(filepath + '_logo.png', x = None, y = None, w = 200, h = 70, type = '', link = '')
    pdf.image(filepath + '_hydro.png', x = None, y = None, w = 200, h = 70, type = '', link = '')  
    

    pdf.output(filepath + '_report.pdf', 'F')
    #os.remove(filepath + '_ramachandran.png')
    os.remove(filepath + '_phipsi.png')
    os.remove(filepath + '_logo.png')
    os.remove(filepath + '_hydro.png')
    os.remove(filepath + '_seqTable.png')

