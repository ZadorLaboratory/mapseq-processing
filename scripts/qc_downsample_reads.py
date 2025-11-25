#!/usr/bin/env python
#
# Top level processing script. 
# Enter into each subdirectory and run all MAPseq processing steps. 
# Assume all defaults/configs are correct. 
#
# After creating reads, perform downsampling on reads into sub-directories. 
# Then run the rest of the pipeline for subdir. 
#
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


STEPLIST=[ 'reads'      ,
           'aggregated',
           'filtered',
           'readtable' ,
           'collapsed' ,
           'vbctable'  ,
           'vbcfiltered',
           'matrices' 
         ]

STEPMAP={ 'reads'       : 'process_fastq_pairs',
          'aggregated'  : 'aggregate_reads',
          'filtered'    : 'filter_split' ,
          'readtable'   : 'make_readtable',
          'collapsed'   : 'align_collapse',
          'vbctable'    : 'make_vbctable',
          'vbcfiltered' : 'filter_vbctable',
          'matrices'    : 'make_matrices'
        }


DIRMAP = { 'process_fastq_pairs': 'reads',
           'aggregate_reads'    : 'aggregated',
           'filter_split'       : 'filtered',
           'make_readtable'     : 'readtable',
           'align_collapse'     : 'collapsed',
           'make_vbctable'      : 'vbctable',
           'filter_vbctable'    : 'vbcfiltered',
           'make_matrices'      : 'matrices'
    }


def qc_process_downsampled_all( config_file, 
                                sampleinfo_file, 
                                infiles , 
                                outdir=None, 
                                force=False,
                                halt=None ):
    '''
    
    
    '''






def process_mapseq_all(config_file, 
                       sampleinfo_file, 
                       infiles , 
                       outdir=None, 
                       force=False,
                       halt=None ):    
    '''    
    performs end-to-end default processing. 
    executes each pipeline script in an external process. 
    
    '''
    global STEPLIST
    
    logging.info(f'{len(infiles)} input files. config={config_file} sampleinfo={sampleinfo_file}, outdir={outdir}, force={force}')
    cp = ConfigParser()
    cp.read(config_file)
    project_id = cp.get('project','project_id')

    if outdir is None:
        outdir = os.path.abspath('./')
  
    logging.debug(f'exe={sys.executable} sys.argv={sys.argv}')
    (dirpath, base, ext) = split_path(sys.argv[0])
    logging.debug(f'script_dir = {dirpath}')

    if halt is not None:
        newstep = []
        for step in STEPLIST:
            if step != halt:
                newstep.append(step)
            elif step == halt:
                newstep.append(step)
                logging.debug(f'found halting step {halt}. breaking. ')
                break
        logging.debug(f'new STEPLIST={STEPLIST}')
        STEPLIST = newstep

    for step in STEPLIST:
        runstep = True
        soutdir = None
        soutfile = None
        sprog = STEPMAP[step]
        sname = DIRMAP[sprog]
        logging.debug(f'handling step={step} sprog={sprog} sname={sname}')
        
        if sname == 'matrices':
            soutdir = os.path.join(outdir, f'{sname}.out/')
        else:
            soutfile = os.path.join(outdir, f'{sname}.out/{project_id}.{sname}.tsv')
    
        # define infile
        if sname != 'reads':
            instep = STEPLIST [ STEPLIST.index(step) - 1 ]
            inprog = STEPMAP[instep]
            insname = DIRMAP[inprog]
            infile = os.path.join( outdir, f'{insname}.out/{project_id}.{insname}.tsv')
        
        log_file = os.path.join(outdir, f'{step}.log')
        cmd = [ os.path.join(dirpath, f'{sprog}.py'),
           '-d',
           '-c', config_file, 
           '-L', log_file,
           ]

        if soutfile is not None:
            cmd.append('-o')
            cmd.append(soutfile)
        
        if soutdir is not None:
            cmd.append('-O')
            cmd.append(soutdir)            
        
        if sname == 'readtable':
            cmd.append('-s')
            cmd.append(sampleinfo_file)
                    
        if sname == 'reads':
            cmd.append('-s')
            cmd.append(sampleinfo_file)
            for fn in infiles:
                cmd.append(fn)
        else:
            cmd.append(infile)
            if not os.path.exists(infile):
                logging.error(f'required infile {infile} does not exist. Exitting.')
                sys.exit(1)
        logging.debug(f"made command={' '.join(cmd)}")

        
        if sname == 'matrices':
            logging.debug(f'make_matrices. outfile not known, proceed...')
            runstep = True
        elif os.path.exists(soutfile):
            logging.info(f'soutfile={soutfile} exists. runstep=False')
            runstep = False
            
        if runstep:
            logging.debug(f'will run {cmd} Logging to ')
            start = dt.datetime.now()
            cp = run_command_shell(cmd)
            end = dt.datetime.now()
            elapsed =  end - start
            logging.debug(f"ran cmd='{cmd}' return={cp.returncode} {elapsed.seconds} seconds.")        
        else:
            logging.debug(f'Output exists, skipping.')
        

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

    parser.add_argument('-H','--halt', 
                        metavar='halt',
                        required=False,
                        default=None,
                        type=str, 
                        help=f'Stage name to stop after: {STEPLIST}')

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
    
    qc_process_downsampled_all( config_file=args.config, 
                                sampleinfo_file=args.sampleinfo, 
                                infiles=args.infiles , 
                                outdir=outdir, 
                                force=args.force, 
                                halt = args.halt
                                )

    logging.info('Done process_all.') 
 
   