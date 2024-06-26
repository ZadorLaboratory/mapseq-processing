#!/usr/bin/env python
#
#   processes per-barcode fasta files. 
#
#    awk "NR%2==0" BC${BCidx[$i]} | grep -v N | cut -b 1-44 | sort | uniq -c | sort -nr > Mseq204_YC_inj2processedBC${BCidx[$i]}.txt
#    #split output files into two files per index, one that is containing the read counts of each unique sequnce, the other the unique sequences themselves.
#    awk '{print $1}' Mseq204_YC_inj2processedBC${BCidx[$i]}.txt > Mseq204_YC_inj2BC${BCidx[$i]}_counts.txt
#    awk '{print $2}' Mseq204_YC_inj2processedBC${BCidx[$i]}.txt > Mseq204_YC_inj2_BC${BCidx[$i]}seq.txt
#
#   BC1.fasta  ->   BC1_processed.tsv
#
#    import resource, sys
#    resource.setrlimit(resource.RLIMIT_STACK, (2**29,-1))
#    sys.setrecursionlimit(10**6)
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

    parser.add_argument('-L','--logfile', 
                    metavar='logfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Logfile for subprocess.')


    parser.add_argument('-t','--threads', 
                        metavar='threads',
                        required=False,
                        default=None,
                        type=int, 
                        help='Handle each input file concurrently in a separate process.')

    parser.add_argument('-s','--sampleinfo', 
                        metavar='sampleinfo',
                        required=True,
                        default=None,
                        type=str, 
                        help='XLS sampleinfo file. ')

    parser.add_argument('-e','--expid', 
                    metavar='expid',
                    required=False,
                    default='M001',
                    type=str, 
                    help='explicitly provided experiment id')

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='out file for all merged info')  

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. input file base dir if not given.')     

    parser.add_argument('infiles',
                        metavar='infiles',
                        nargs ="+",
                        type=str,
                        help='SSI fasta file(s)')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    cp.read(args.config)   
    cdict = format_config(cp)
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infiles={args.infiles}')
    
    outdir = None
    outfile = None
    if args.outdir is not None:
        outdir = os.path.abspath(args.outdir)
        logging.debug(f'making missing outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)
    else:
        afile = args.infiles[0]
        filepath = os.path.abspath(afile)    
        dirname = os.path.dirname(filepath)
        outdir = dirname
    
    if args.outfile is not None:
        outfile = args.outfile
    else:
        outfile = f'{outdir}/{args.expid}.all.tsv'
    
    #sampdf = load_sample_info(config, args.sampleinfo)
    #logging.debug(f'\n{sampdf}')    
    # sample info passed as filename rather than dataframe to function
    #    because it needs to be a command line argument to sub-processes. 
    
    if args.logfile is not None:
        log = logging.getLogger()
        FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
        formatter = logging.Formatter(FORMAT)
        logStream = logging.FileHandler(filename=args.logfile)
        logStream.setFormatter(formatter)
        log.addHandler(logStream)
    
    
    outdf = process_ssifasta_files(cp, args.sampleinfo, args.infiles, numthreads=args.threads, outdir=outdir)
    logging.info(f'saving output to {outfile}')
    outdf.to_csv(outfile, sep='\t')
    
    # create filtered subset based on reads/UMI ratio. 
    dir, base, ext = split_path(outfile)
    tdfoutfile = f'{dir}/{base}.filtered{ext}'
    logging.info(f'saving read thresholded output to {tdfoutfile}')
    tdf = read_threshold_all(cp, outdf )
    tdf.to_csv(tdfoutfile, sep='\t')
    
    
    
    
    
    
