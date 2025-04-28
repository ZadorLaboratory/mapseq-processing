#!/usr/bin/env python
#
# Top level processing script. 
# Enter into each subdirectory and run all MAPseq processing steps. 
# Assume all defaults/configs are correct. 
#
import argparse
import logging
import os
import sys

from configparser import ConfigParser

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import *
from mapseq.utils import *
from mapseq.stats import *  


STEPLIST=[ 'process_fastq_pairs',
           'aggregate_reads',
           'filter_split' ,
           'align_collapse',
           'make_readtable',
           'make_vbctable',
           'filter_vbctable',
           #'make_matrices'
           ]


DIRMAP = { 'process_fastq_pairs': 'reads',
           'aggregate_reads'    : 'aggregated',
           'filter_split'       : 'filtered',
           'align_collapse'     : 'collapsed',
           'make_readtable'     : 'readtable',
           'make_vbctable'      : 'vbctable',
           'filter_vbctable'    : 'vbcfiltered',
           #'make_matrices'      : 'matrices'
    }


def process_mapseq_all(config_file, sampleinfo_file, infiles , outdir=None, force=False ):    
    '''    
    performs end-to-end default processing. 
    executes each pipeline script in an external process. 
    
    '''
    logging.info(f'{len(infiles)} input files. config={config_file} sampleinfo={sampleinfo_file}, outdir={outdir}, force={force}')
    cp = ConfigParser()
    cp.read(config_file)
    expid = cp.get('project','project_id')
  
    logging.debug(f'exe={sys.executable} sys.argv={sys.argv}')
    (dirpath, base, ext) = split_path(sys.argv[0])
    logging.debug(f'script_dir = {dirpath}')

    for step in STEPLIST:
        sname = DIRMAP[step]
        logging.debug(f'handling step={step} sname={sname}')
        outfile = os.path.join(outdir, f'{sname}.out/{expid}.{sname}.tsv')
    
        # define infile
        if sname != 'reads':
            instep = STEPLIST [ STEPLIST.index(step) - 1 ]
            insname = DIRMAP[instep] 
            infile = os.path.join( outdir, f'{insname}.out/{expid}.{insname}.parquet')
        
        log_file = os.path.join(outdir, f'{step}.log')
        cmd = [ os.path.join(dirpath, f'{step}.py'),
           '-d',
           '-c', config_file, 
           '-L', log_file,
           '-o', outfile
           ]
        if sname == 'readtable':
            cmd.append('-s')
            cmd.append(sampleinfo_file)
                    
        if sname == 'reads':
            for fn in infiles:
                cmd.append(fn)
        else:
            cmd.append(infile)
            if not os.path.exists(infile):
                logging.error(f'required infile {infile} does not exist. Exitting.')
                sys.exit(1)
        logging.debug(f"made command={' '.join(cmd)}")

        if os.path.exists(outfile):
            logging.info(f'outfile={outfile} exists. Continue...')
        else:
            logging.debug(f'will run {cmd} Logging to ')
            start = dt.datetime.now()
            cp = run_command_shell(cmd)
            end = dt.datetime.now()
            elapsed =  end - start
            logging.debug(f"ran cmd='{cmd}' return={cp.returncode} {elapsed.seconds} seconds.")        
        

            
             
def old_code(config_file, sampleinfo_file, infiles , outdir=None, force=False ):
    
    # process_fasq_pairs.py

    log_file = os.path.join(outdir, f'process_fastq_pairs.log')
    cmd = [ os.path.join(dirpath, 'process_fastq_pairs.py'),
           '-d',
           '-c', config_file, 
           '-L', log_file ,
           '-o', outfile
           ]

    for fn in infiles:
        cmd.append( fn)
    logging.debug(f'will run {cmd} Logging to ')
    start = dt.datetime.now()
    cp = run_command_shell(cmd)
    end = dt.datetime.now()
    elapsed =  end - start
    logging.debug(f"ran cmd='{cmd}' return={cp.returncode} {elapsed.seconds} seconds.")
    
    # aggregate_reads.py
    
    # filter_split.py
    
    # align_collapse.py
    
    # make_readtable.py
    
    # make vbctable.py
    
    # filter_vbctable.py
    
    # make_matrices.py








if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARN)
    
    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')

    parser.add_argument('-c','--config', 
                        metavar='config',
                        required=True,
                        default=os.path.expanduser('~/git/mapseq-processing/etc/mapseq.conf'),
                        type=str, 
                        help='config file.')
    
    parser.add_argument('-s','--sampleinfo', 
                        metavar='sampleinfo',
                        required=True,
                        default=None,
                        type=str, 
                        help='XLS sampleinfo file. ')  

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. output file base dir if not given.')
    
    parser.add_argument('-f','--force', 
                    action="store_true", 
                    default=False, 
                    help='Recalculate even if output exists.') 

    parser.add_argument('infiles' ,
                        metavar='infiles', 
                        type=str,
                        nargs='*',
                        default=None, 
                        help='Initial FASTQ input to process. ')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        loglevel = 'debug'
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   
        loglevel = 'info'

    logging.debug(f'indirs={args.infiles}')
    
    cp = ConfigParser()
    cp.read(args.config)
    cdict = format_config(cp)
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infiles={args.infiles}')

    # set outdir / outfile
    outdir = os.path.abspath('./')
    if args.outdir is not None:
        outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    
    datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    sh = StatsHandler(outdir=outdir, datestr=datestr)

    process_mapseq_all( config_file=args.config, 
                        sampleinfo_file=args.sampleinfo, 
                        infiles=args.infiles , 
                        outdir=outdir, 
                        force=args.force
                        )

 
 
   