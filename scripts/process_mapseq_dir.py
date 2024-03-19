#!/usr/bin/env python
#
#
#

import argparse
import logging
import os
import sys
import traceback

from configparser import ConfigParser

import pandas as pd

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import *
from mapseq.utils import *
from mapseq.stats import *


def get_ordered_fastq_filestring(fastqdir):
    '''
    XXXXX needs to be fixed to handle multi-pair readfiles. 
    
    '''
    logging.debug(f'figuring fastq order in {fastqdir} ...')
    fqfiles = glob.glob(f'{fastqdir}/*fastq*')
    fqfiles.sort()
    s = ''
    for f in fqfiles:
        s += f' {f} '
    return s

def get_bcfasta_filestring(fastadir):    
    logging.debug(f'getting BC*.fasta in {fastadir} ...')
    fafiles = glob.glob(f'{fastadir}/BC*.fasta')
    s = ''
    for f in fafiles:
        s += f' {f} '
    return s


def process_mapseq_dir(dirpath, expid, loglevel, conffile=None):
    '''
    E.g. 
   
    process_fastq_pairs.py -v -o fastq.2.out/M253.all.fasta fastq/M253_CR_S1_R* 
    align_collapse.py -v -b 30 -m 3 -O collapse.2.out fastq.2.out/M253.all.fasta   
    process_fasta.py -v -s M253_sampleinfo.xlsx -O fasta.out collapse.out/M253.all.collapsed.fasta 
    process_ssifasta.py -v -s M253_sampleinfo.xlsx -o ssifasta.out/M253.merged.all.tsv  -O ssifasta.out
    process_merged.py -v -s M253_sampleinfo.xlsx -e M253.default -O merged.out ssifasta.out/M253.all.tsv
    
    
    '''
    logging.info(f'dirpath={dirpath} expid={expid}, loglevel={loglevel}, conffile={conffile}')
    
    dirpath = os.path.abspath(dirpath)
    if not os.path.exists(dirpath):
        sys.exit(f'Experiment directory {dirpath} does not exist.')
    logging.info(f'processing experiment dir: {dirpath}')

    datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    logging.info(f'datestring for this run is {datestr}')
   
    cp = ConfigParser()
    if conffile is None:
        cp = get_default_config()
    else:
        cp.read(conffile) 
    
    cdict = format_config(cp)    
    logging.debug(f'Running with config. {args.config}: {cdict}')
    
    rundir = f'{dirpath}/run.{datestr}.out'
    os.makedirs(rundir, exist_ok=True)
    logging.debug(f'making rundir: {rundir} ')
    configfile = f'{rundir}/mapseq.conf'
    configfile = write_config( cp, configfile, timestamp=True, datestring=datestr)
    logging.info(f'wrote {configfile} for usage by all.')      
    
    samplefile = f'{dirpath}/{expid}_sampleinfo.xlsx'
    if not os.path.exists(samplefile):
        sys.exit(f'sampleinfo file {samplefile} does not exist.') 

    # all scripts...
    scriptroot = f'{gitpath}/scripts'
    
    
    #
    # process_fastq_pairs
    #
    prog = f'{scriptroot}/process_fastq_pairs.py'
    # fastq_pairs
    logging.info(f'running process_fastq_pairs.py ... ')
    fqfiles = get_ordered_fastq_filestring(f'{dirpath}/fastq')
    outdir =  f'{rundir}/fastq.out'
    outfile = f'{outdir}/{expid}.all.fasta'
    logfile = f'{outdir}/fastq.log'
    logging.debug(f'prog={prog} logfile={logfile} configfile={configfile} outdir={outdir} fqfiles={fqfiles}')
    cmd = [ prog ,
            '-d' ,
            '-L', logfile , 
            '-c', configfile , 
            '-o' , outfile,
            fqfiles
           ] 
    try:
        result = run_command_shell(cmd)           
    except NonZeroReturnException:
        logging.warning(f'NonZeroReturn Exception: {cmd}')
        sys.exit(result.returncode)
    if result.stderr is not None:
        logging.warn(f"got stderr: {result.stderr}")    
    logging.info(f'done process_fastq_pairs.py')
    
    
    #
    # align_collapse
    #
    prog = f'{scriptroot}/align_collapse.py'
    # fastq_pairs
    logging.info(f'running align_collapse.py ... ')
    infile = f'{rundir}/fastq.out/{expid}.all.fasta' 

    outdir =  f'{rundir}/collapse.out'
    logfile = f'{outdir}/collapse.log'
    logging.debug(f'prog={prog} logfile={logfile} configfile={configfile} outdir={outdir} ')
    cmd = [ prog ,
            '-d' ,
            '-L', logfile , 
            '-c', configfile , 
            '-O' , outdir,
            infile
           ] 
    try:
        result = run_command_shell(cmd)           
    except NonZeroReturnException:
        logging.warning(f'NonZeroReturn Exception: {cmd}')
        sys.exit(result.returncode)
    if result.stderr is not None:
        logging.warn(f"got stderr: {result.stderr}")    
    logging.info(f'done align_collapse.py')    


    #
    # process_fasta.py
    #
    prog = f'{scriptroot}/process_fasta.py'
    # fastq_pairs
    logging.info(f'running process_fasta.py ... ')
    
    infile = f'{rundir}/collapse.out/{expid}.all.fasta' 
    outdir =  f'{rundir}/fasta.out'
    logfile = f'{outdir}/fasta.log'
    logging.debug(f'prog={prog} logfile={logfile} configfile={configfile} infile={infile} outdir={outfile} ')
    cmd = [ prog ,
            '-d' ,
            '-L', logfile , 
            '-c', configfile , 
            '-O' , outdir,
            infile
           ] 
    try:
        result = run_command_shell(cmd)           
    except NonZeroReturnException:
        logging.warning(f'NonZeroReturn Exception: {cmd}')
        sys.exit(result.returncode)
    if result.stderr is not None:
        logging.warn(f"got stderr: {result.stderr}")    
    logging.info(f'done process_fasta.py')       


    #
    # process_ssifasta.py
    #
    prog = f'{scriptroot}/process_ssifasta.py'
    logging.info(f'running process_ssifasta.py ... ')
    
    indir = f'{rundir}/fasta.out'
    bcfasta_str = get_bcfasta_filestring(indir) 
    outdir =  f'{rundir}/ssifasta.out'
    outfile = f'{rundir}/ssifasta.out/{expid}.all.tsv'
    logfile = f'{outdir}/ssifasta.log'
    logging.debug(f'prog={prog} logfile={logfile} configfile={configfile} indir={indir} outfile={outfile} outdir={outfile} ')
    cmd = [ prog ,
            '-d' ,
            '-L', logfile , 
            '-c', configfile , 
            '-O' , outdir,
            '-o', outfile , 
            bcfasta_str
           ] 
    try:
        result = run_command_shell(cmd)           
    except NonZeroReturnException:
        logging.warning(f'NonZeroReturn Exception: {cmd}')
        sys.exit(result.returncode)
    if result.stderr is not None:
        logging.warn(f"got stderr: {result.stderr}")    
    
    logging.info(f'done process_ssifasta.py')     
    

    #
    # process_merged.py
    #
    prog = f'{scriptroot}/process_merged.py'
    logging.info(f'running process_merged.py ... ')   
    infile = f'{rundir}/ssifasta.out/{expid}.all.tsv'
    outdir =  f'{rundir}/merged.out'
    logfile = f'{outdir}/merged.log'
    logging.debug(f'prog={prog} logfile={logfile} configfile={configfile} infile={infile} outdir={outfile} ')
    cmd = [ prog ,
            '-d' ,
            '-L', logfile ,
            '-e', expid ,  
            '-c', configfile , 
            '-O' , outdir, 
            infile
           ] 
    try:
        result = run_command_shell(cmd)           
    except NonZeroReturnException:
        logging.warning(f'NonZeroReturn Exception: {cmd}')
        sys.exit(result.returncode)
    if result.stderr is not None:
        logging.warn(f"got stderr: {result.stderr}")    
    
    logging.info(f'done process_merged.py')   
    
    
    
    
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
                        required=False,
                        default=os.path.expanduser('~/git/mapseq-processing/etc/mapseq.conf'),
                        type=str, 
                        help='config file.')    

    parser.add_argument('-e','--expid', 
                    metavar='expid',
                    required=True,
                    type=str, 
                    help='Explicitly provided experiment id, e.g. M205')

#    parser.add_argument('-O','--outdir', 
#                    metavar='outdir',
#                    required=False,
#                    default=None, 
#                    type=str, 
#                    help='outdir. input file base dir if not given.')     

#    parser.add_argument('-s','--sampleinfo', 
#                        metavar='sampleinfo',
#                        required=True,
#                        default=None,
#                        type=str, 
#                        help='XLS sampleinfo file. ')

    parser.add_argument('indir',
                        metavar='indir',
                        nargs ="?",
                        type=str,
                        help='Standardized MAPseq data dir, with fastq/ subdir and <EXPID>_sampleinfo.xlsx')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    cp.read(args.config)
    cdict = format_config(cp)    
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infiles={args.indir}')
    
    loglevel = logging.getLogger().getEffectiveLevel()
    
    process_mapseq_dir( args.indir, args.expid, loglevel=loglevel, conffile=args.config)
    
    
    







        
