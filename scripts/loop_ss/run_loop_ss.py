import os
import sys
import importlib
from rvrsprot.loop import loop_ss
import datetime 


def main():
    path = sys.argv[1]
    para = importlib.machinery.SourceFileLoader('para', path).load_module()
    print('Task: ' + path)

    time_tag = datetime.datetime.now().strftime('%Y%m%d-%H%M%S') 
    outdir = para.workdir + 'loop_ss_' + time_tag + '/'
    if para.outdir:
        outdir = para.workdir + para.outdir
    os.makedirs(outdir, exist_ok=True)
    loop_ss.run_loop_ss(outdir, para.workdir + para.target_file, para.loop_topo_sels, para)

    return

if __name__=='__main__':
    main()