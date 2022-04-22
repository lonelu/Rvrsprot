import os
import shlex
import shutil
import subprocess
import qbits

createPDS = os.path.dirname(__file__) + '/createPDS'
master = os.path.dirname(__file__) + '/master'

def master_query_depre(pdb_path, targetList, rmsdCut=1., topN=None, 
                 outfile=None, clobber=False, outdir = None):
    """Execute a MASTER query for a given PDB file.

    Parameters
    ----------
    pdb_path : str
        Path to PDB file for which to execute the MASTER query.
    targetList : str
        Path to file containing paths to PDS files for target structures.
    rmsdCut : float, optional
        RMSD threshold at which to permit MASTER matches.
    topN : int, optional
        Number of top database hits to include in the output.  If None, 
        an unlimited number of hits are permitted.
    outfile : str, optional
        Path to a file to which MASTER output should be redirected.
    clobber : bool, optional
        If True, clobber MASTER output if it already exists.
    """
    cwd = os.getcwd()
    pdb_path = os.path.realpath(pdb_path)
    pdb_dir = os.path.dirname(pdb_path)
    if not clobber:
        if os.path.exists(pdb_dir + '/match.txt'):
            return
    os.chdir(pdb_dir)
    cmd = createPDS + ' --type query --pdb {}'.format(os.path.basename(pdb_path))
    if outfile:
        with open(outfile, 'a') as f:
            subprocess.run(shlex.split(cmd), stdout=f)
    else:
        subprocess.run(shlex.split(cmd))
    if not outdir:
        outdir = 'match_'
    pre, ext = os.path.splitext(pdb_path)
    pds_path = pre + '.pds'
    if topN is not None:
        cmd = (master  +' --query {} --targetList {} --topN {} --rmsdCut {} '
               '--seqOut seq.txt --matchOut match.txt --structOut {}').format(pds_path, 
                targetList, str(topN),str(rmsdCut), outdir)
    else:
        cmd = (master + ' --query {} --targetList {} --rmsdCut {} '
               '--seqOut seq.txt --matchOut match.txt --structOut {}').format(pds_path, 
                                                               targetList, str(rmsdCut), outdir)
    # make sure the same query has not been made
    parent_seed = None
    parent_dir = os.path.dirname(pdb_dir)
    while 'seed.pdb' in os.listdir(parent_dir):
        with open(parent_dir + '/seed.pdb', 'r') as f0:
            with open(pdb_path, 'r') as f1:
                if f0.read() == f1.read():
                    parent_seed = parent_dir + '/seed.pdb'
                    break
        parent_dir = os.path.dirname(parent_dir)
    if parent_seed is not None:
        for filename in ['/seed.pds', '/match.txt', '/seq.txt']:
            shutil.copyfile(parent_dir + filename, pdb_dir + filename)
    else:
        if outfile:
            with open(outfile, 'a') as f:
                subprocess.run(shlex.split(cmd), stdout=f)
        else:
            subprocess.run(shlex.split(cmd))
    os.chdir(cwd)


def master_query(outdir, pdb_path, targetList, rmsdCut=1., topN=None, outfile=None):
    """Execute a MASTER query for a given PDB file.

    Parameters
    ----------
    pdb_path : str
        Path to PDB file for which to execute the MASTER query.
    targetList : str
        Path to file containing paths to PDS files for target structures.
    rmsdCut : float, optional
        RMSD threshold at which to permit MASTER matches.
    topN : int, optional
        Number of top database hits to include in the output.  If None, 
        an unlimited number of hits are permitted.
    outfile : str, optional
        Path to a file to which MASTER output should be redirected.
    clobber : bool, optional
        If True, clobber MASTER output if it already exists.
    """
    pds_path = outdir + os.path.basename(pdb_path) + '.pds'

    cmd = createPDS + ' --type query --pdb {} --pds {}'.format(pdb_path, pds_path)
    if outfile:
        with open(outfile, 'a') as f:
            subprocess.run(shlex.split(cmd), stdout=f)
    else:
        subprocess.run(shlex.split(cmd))
    
    if topN is not None:
        cmd = (master + ' --query {} --targetList {} --topN {} --rmsdCut {} '
               '--seqOut {} --matchOut {} --structOut {}').format(pds_path, 
                targetList, str(topN),str(rmsdCut), outdir + 'seq.txt', outdir + 'match.txt', outdir)
    else:
        cmd = (master + ' --query {} --targetList {} --rmsdCut {} '
               '--seqOut {} --matchOut {} --structOut {}').format(pds_path, 
                targetList, str(rmsdCut), outdir + 'seq.txt', outdir + 'match.txt', outdir)
    if outfile:
        with open(outfile, 'a') as f:
            subprocess.run(shlex.split(cmd), stdout=f)
    else:
        subprocess.run(shlex.split(cmd))
    return 


def master_query_loop(pdb_path, targetList, rmsdCut=1., gapLen=10,  
                      topN=None, outdir=None, outfile=None):
    """Execute a MASTER query for a given PDB file.

    Parameters
    ----------
    pdb_path : str
        Path to PDB file for which to execute the MASTER query.
    targetList : str
        Path to file containing paths to PDS files for target structures.
    rmsdCut : float, optional
        RMSD threshold at which to permit MASTER matches.
    gapLen : int or str, optional
        Maximum length of gap between discontiguous segments in MASTER targets.
    topN : int, optional
        Number of top database hits to include in the output.  If None, 
        an unlimited number of hits are permitted.
    outdir : str, optional
        Directory in which to output PDB files for loop structures.
    outfile : str, optional
        Path to a file to which MASTER output should be redirected.
    """
    cwd = os.getcwd()
    pdb_path = os.path.realpath(pdb_path)
    pdb_dir = os.path.dirname(pdb_path)
    os.chdir(pdb_dir)
    cmd = createPDS + ' --type query --pdb {}'.format(os.path.basename(pdb_path))
    if outfile:
        with open(outfile, 'a') as f:
            subprocess.run(shlex.split(cmd), stdout=f)
    else:
        subprocess.run(shlex.split(cmd))
    if not outdir:
        outdir = 'loops_' + str(gapLen)
    pre, ext = os.path.splitext(pdb_path)
    pds_path = pre + '.pds'
    if topN is not None:
        cmd = (master + ' --query {} --targetList {} --topN {} '
               '--rmsdCut {} --gapLen {} --outType wgap '
               '--matchOut match.txt --seqOut seq.txt '
               '--structOut {}').format(pds_path, targetList, 
                                        str(topN), str(rmsdCut), 
                                        str(gapLen), outdir)
    else:
        cmd = (master + ' --query {} --targetList {} '
               '--rmsdCut {} --gapLen {} --outType wgap '
               '--matchOut match.txt --seqOut seq.txt '
               '--structOut {}').format(pds_path, targetList, 
                                        str(rmsdCut), str(gapLen), 
                                        outdir)
    if outfile:
        with open(outfile, 'a') as f:
            subprocess.run(shlex.split(cmd), stdout=f)
    else:
        subprocess.run(shlex.split(cmd))
    os.chdir(cwd)


def qbits_search(query_pdb_path, query_full_pdb_path, query_dir, chains_dict_path, 
                 outdir, window_size=10, rmsd=1.5, top=5, sec_struct=None, 
                 antiparallel=False, min_nbrs=1, contiguous=False):
    """Parse MASTER matches and search for neighbors to generate qbit reps.

    Parameters
    ----------
    query_pdb_path : str
        Path to PDB file of query structure.
    query_full_pdb_path : str
        Path to PDB file of full structure including query structure.
    chains_dict_path : str
        Path to pickled dict containing chain information for qbits search.
    outdir : str
        Directory in which to write output files from the qbits search.
    window_size : int, optional
        Window size (in amino acid residues) for qbits search.
    rmsd : float, optional
        RMSD cutoff at which to include neighbors in the qbits search.
    top : int, optional
        Number of top qbit reps to output.
    sec_struct : str, optional
        String index for secondary structural elements to which to 
        restrict the qbits search (e.g. "H" for alpha-helices).
    antiparallel : bool, optional
        Keep only antiparallel window fragments in the Qbits search.
    min_nbrs : int, optional
        Minimum number of neighbors for a residue to include it in the 
        Qbit reps.
    contiguous : bool, optional
        If True, limit to contiguous qbit reps.
    """
    query_pdb_path = os.path.realpath(query_pdb_path)
    query_full_pdb_path = os.path.realpath(query_full_pdb_path)
    #query_dir = os.path.dirname(query_pdb_path)
    match_path = query_dir + 'match.txt'
    seq_path = query_dir + 'seq.txt'
    # parse MASTER matches
    p = qbits.parse.Parse(query_full_pdb_path, query_pdb_path, match_path, 
                          seq_path)
    p.parse(outdir=outdir + '/', show_parsing_progress=False)
    # search through matches for neighbors
    qs = qbits.search.Qsearch(window_size=window_size, rmsd=rmsd, 
                              min_nbrs=min_nbrs)
    qs.sift(p, chains_dict_path, include_contiguous=True, 
            contiguous_only=contiguous, filter_ss=bool(sec_struct), 
            sec_struct=sec_struct, antiparallel=antiparallel, 
            min_nbrs=min_nbrs)
    if contiguous:
        qs.get_contiguous_qbit_reps(p, top=top, outdir=outdir + '/qbit_reps/')
    else:
        qs.get_qbit_reps(p, top=top, outdir=outdir + '/qbit_reps/')

