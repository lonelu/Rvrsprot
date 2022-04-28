from pymol import cmd
import glob
import os
import shutil
import shlex
import subprocess
import prody as pr

def prepare_fasta_by_prody():
    '''
    prepare fasta file of each chain each tm protein.
    It turned out the prody have troubles to extract fasta. Using pymol function is better. 
    '''
    workdir = '/mnt/e/DesignData/tm/database/'
    pdbs = []
    failed = []
    for file in os.listdir(workdir + 'tms_filter3.5/'):
        if '.pdb' not in file:
            continue
        try:
            pdb = pr.parsePDB(workdir + 'tms_filter3.5/' + file)
            pdbs.append(pdb)

        except:
            failed.append(file)

    with open(workdir + 'allseq.fasta', 'w') as f:
        for pdb in pdbs:
            for chid in pdb.select('protein').getChids():
                c = pdb.select('protein and chid ' + chid)
                f.write('>' + pdb.getTitle() + '_' + chid + '\n')
                f.write(c.select('name CA').getSequence() + '\n')


def prepare_fasta_by_pymol(workdir, outdir):
    '''
    prepare fasta file of each chain each tm protein by pymol.
    '''
    file_dict = {}
    all_files = []

    for _file in glob.glob(os.path.join(workdir, '*.pdb')):
        all_files.append(_file)

    for _file in all_files:
        file_dict[_file] = []
        cmd.load(_file)
        objs = cmd.get_names()
        
        for obj in objs:
            obj_prot = obj + '_prot'
            count = cmd.count_atoms('model %s and polymer.protein' % (obj))
            if count <= 10:
                continue

            cmd.create(obj_prot ,'model %s and polymer.protein' % (obj))
            
            for ch in cmd.get_chains(obj_prot):
                if len(ch) <= 0:
                    print('Chain has no alphabeta: ' + obj)

                name = obj_prot + '_' + ch
                _count = cmd.count_atoms('model %s and chain %s and polymer.protein' % (obj_prot, ch))
                file_dict[_file].append(name + '#' + str(count))
                if _count < 10:
                    continue
                cmd.create(name, 'model %s and chain %s and polymer.protein' % (obj_prot, ch))
                cmd.save(_file.split('.')[0] + '_' + ch + '.fasta', name)
        cmd.delete("all")
    
    with open(outdir + '_summary.tsv', 'w') as f:
        for k in file_dict.keys():
            f.write(k + '\t' + '\t'.join(file_dict[k]) + '\n')
    return


def comb_fasta(workdir, outdir):
    '''
    Combine all the fasta files extracted by pymol.
    '''
    with open(outdir + '_all_seq.fasta', 'ab') as wfd:
        for file in os.listdir(workdir):
            if not '.fasta' in file:
                continue
            with open(workdir + file, 'rb') as fd:
                shutil.copyfileobj(fd, wfd) 
    return


def comb_sep_fasta(workdir, outdir):
    '''
    Combine the fasta files with the same id.
    '''
    os.makedirs(outdir, exist_ok=True)
    db_dict = {}
    for file in os.listdir(workdir):
        if not '.fasta' in file:
            continue
        key = file.split('_')[0]
        if key in db_dict.keys():
            db_dict[key].append(workdir + file)
        else:
            db_dict[key] = [workdir + file]

    for key in db_dict.keys():
        count = len(db_dict[key])
        if count > 1:
            count = 2
        outdir_count = outdir + str(count) + '/'
        os.makedirs(outdir_count, exist_ok=True)
        with open(outdir_count + key + '.fasta', 'ab') as wfd:
            for file in db_dict[key]:
                with open(file, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd) 

    # with open(outdir + 'summary.tsv', 'w') as f:
    #     for key in db_dict.keys():
    #         f.write(key + '\t' + str(len(db_dict[key])) + '\n')
    return


def run_mmseqs(workdir, outdir):
    os.makedirs(outdir, exist_ok=True)
    outfile = outdir + 'mmseq_out.txt'
    for file in os.listdir(workdir):
        if not '.fasta' in file:
            continue
        cmd = 'mmseqs easy-cluster -c 0.9 ' + workdir + file + ' ' + outdir + file.split('.')[0] + ' ' + outdir + 'temp'
        with open(outfile, 'a') as f:
            subprocess.run(shlex.split(cmd), stdout=f)
    return


def extract_mmseq(workdir, outdir):
    '''
    summarize mmseq out information.
    '''
    info = {}
    for file in os.listdir(workdir):
        if not '.tsv' in file:
            continue
        chids = {}
        with open(workdir + file, 'r') as f:
            for line in f.readlines():
                k = line.split('\t')[0]
                if k in chids.keys():
                    chids[k] += 1
                else:
                    chids[k] = 1
        info[file.split('_')[0]] = chids

    with open(outdir + '_2_mmseq_summary.tsv', 'w') as f:
        for key in info.keys():
            chids = '\t'.join([k + '\t' + str(info[key][k]) for k in info[key].keys()])
            count = len(info[key].keys())
            f.write(key + '\t' + str(count) + '\t' + chids + '\n')
    return


def unique_chid_redundant_chid(workdir, workdir_single, outdir):
    '''
    After mmseq search, we want to have a table to keep the unique chids and redundant chid for sequence analysis.
    '''
    unique = set()
    redundance = set()
    for file in os.listdir(workdir):
        if not '.tsv' in file:
            continue
        
        with open(workdir + file, 'r') as f:
            for line in f.readlines():
                k1 = line.split('\t')[0].split('_')[0] + '_' + line.split('\t')[0].split('_')[-1]
                k2 = line.split('\t')[1].split('_')[0] + '_' + line.split('\t')[1].split('_')[-1][:-1]
                unique.add(k1)
                if not k2 in unique:
                    redundance.add(k2)

    #Add single chid into unique.
    db_dict = {}
    for file in os.listdir(workdir_single):
        if not '.fasta' in file:
            continue
        key = file.split('_')[0]
        if key in db_dict.keys():
            db_dict[key].append(file)
        else:
            db_dict[key] = [file]

    for key in db_dict.keys():
        count = len(db_dict[key])
        if count == 1:
            unique.add(db_dict[key][0].split('.')[0])


    with open(outdir + '_unique.tsv', 'w') as f:
        for k in unique:
            f.write(k + '\n')
    with open(outdir + '_redundant.tsv', 'w') as f:
        for k in redundance:
            f.write(k + '\n')    
    return


workdir = '/mnt/e/DesignData/tm/database/tms_filter35/'  

outdir = workdir + 'seq/'
cmd.cd(outdir)
os.makedirs(outdir, exist_ok=True)
print(outdir)
prepare_fasta_by_pymol(workdir, outdir)

comb_fasta(outdir, workdir)

comb_sep_fasta(outdir, workdir + 'seq_combinations/')

run_mmseqs(workdir + 'seq_combinations/2/', workdir + 'seq_combinations/2_mmseq/')
extract_mmseq(workdir + 'seq_combinations/2_mmseq/', workdir + 'seq_combinations/')


unique_chid_redundant_chid(workdir + 'seq_combinations/2_mmseq/', workdir + 'seq_combinations/1/', workdir + 'seq_combinations/')