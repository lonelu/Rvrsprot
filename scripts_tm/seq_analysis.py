import sys
import numpy as np
from metalprot.basic import constant
import logomaker
import matplotlib.pyplot as plt

'''
python /mnt/e/GitHub_Design/Metalprot/scripts_tm/seq_analysis.py /mnt/e/DesignData/tm/ test_0.pdb.pds.opm.seq
'''

def _get_seq_rmsd(seqfile):
    with open(seqfile, 'r') as f:
        lines = f.read().split('\n')
    all_seqs = []
    all_rmsds = []
    for line in lines:
        if len(line) > 0:
            seq =[]
            for res in line.split(' '):
                if res.replace('.','',1).isdigit():
                    all_rmsds.append(float(res))
                if len(res) == 3 and res[0].isalpha():
                    seq.append(res)
                elif len(res) == 4 and res[0] == '[':
                    seq.append(res[1:])
                elif len(res) == 4 and res[-1] == ']':
                    seq.append(res[:-1])
            all_seqs.append(seq)
    return all_seqs, all_rmsds

def _get_loop_candidate_seqs_one_letter(seqs):
    seqs_one_letter = []
    for s in seqs:      
        seqs_one_letter.append(''.join([constant.one_letter_code[res] for res in s]))
    return seqs_one_letter

def _get_pdbs_master_info(matchfile, seqfile, loop_pdbs):
    with open(matchfile, 'r') as f:
        ## A possible bug.
        # rmsds = [float([val for val in line.split(' ') if val != ''][0]) 
        #             for line in f.read().split('\n') if len(line) > 0]
        all_pdss = [[val for val in line.split('/') if '.pds' in val][0] 
                    for line in f.read().split('\n') if len(line) > 0]
    all_seqs, all_rmsds = _get_seq_rmsd(seqfile)

    loop_rmsds = []
    loop_seqs = []
    loop_pdss = []
    for loop_pdb in loop_pdbs:
        if 'wgap' in loop_pdb:
            idx = int(loop_pdb.split('_')[-1][4:-7]) - 1
        else:
            idx = int(loop_pdb.split('_')[-1][5:-7]) - 1
        loop_rmsds.append(all_rmsds[idx])
        loop_seqs.append(all_seqs[idx])
        loop_pdss.append(all_pdss[idx])
    return loop_rmsds, loop_seqs, loop_pdss
    #return loop_pdbs[np.argmin(loop_rmsds)]

def _plot_log(fig, ax, seqs):
    seqs_one_letter = _get_loop_candidate_seqs_one_letter(seqs)
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
    return df


def write_top_aa(workdir, df, top = 5, filename = '_top_seq.txt'):
    '''
    For the logo frequency, output the top aa.
    '''
    with open(workdir + filename, 'w') as f:
        for i in range(df.shape[0]):
            idx = np.argsort(df.iloc[i])
            aas = np.flip(df.columns[idx][-top:])
            f.write(str(i) + '\t' + ''.join([a for a in aas]) + '\n')
    return


def main(workdir, seqfile, out_name = '_info.png'):
    all_seqs, all_rmsds = _get_seq_rmsd(workdir + seqfile)
    fig, (ax1) =plt.subplots(1, 1, figsize=(15, 6))
    df = _plot_log(fig, ax1, all_seqs)
    plt.tight_layout()
    plt.savefig(workdir + out_name )
    plt.close()
    write_top_aa(workdir, df, top = 5, filename = '_top_seq.txt')

if __name__=='__main__':
    main(sys.argv[1], sys.argv[2])
