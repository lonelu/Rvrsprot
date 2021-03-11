import numpy as np
import qbits

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

def _get_loop_candidate_seqs(all_seqs, loop_len, loop_query_wins):
    seqs = []
    for s in all_seqs:
        if len(s) == sum(loop_query_wins) + loop_len:            
            seqs.append(s)
    return seqs

def _get_loop_candidate_seqs_one_letter(seqs):
    seqs_one_letter = []
    for s in seqs:      
        seqs_one_letter.append(''.join([qbits.constants.one_letter_code[res] for res in s]))
    return seqs_one_letter

def _cal_loop_candidate_rmsd(seqfile, loop_len, loop_query_wins):
    all_seqs, all_rmsds = _get_seq_rmsd(seqfile)
    rmsds = []
    for i in range(len(all_seqs)):
        if len(all_seqs[i]) == sum(loop_query_wins) + loop_len:            
            rmsds.append(all_rmsds[i])  
    return min(rmsds), np.median(rmsds)