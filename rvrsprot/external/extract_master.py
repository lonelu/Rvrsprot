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

def _get_loop_candidate_seqs_one_letter(seqs):
    seqs_one_letter = []
    for s in seqs:      
        seqs_one_letter.append(''.join([qbits.constants.one_letter_code[res] for res in s]))
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