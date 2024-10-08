#!/usr/bin/env python
import argparse
import logging
import os
import sys

from configparser import ConfigParser

import pandas as pd

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import *
from mapseq.barcode import *
from mapseq.utils import *
from mapseq.stats import *


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
                        help='out file.')    

    parser.add_argument('-L','--logfile', 
                    metavar='logfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Logfile for subprocess.')

    parser.add_argument('-a','--aligner', 
                    metavar='aligner',
                    required=False,
                    default=None, 
                    type=str, 
                    help='aligner tool  [bowtie | bowtie2]')

    parser.add_argument('-m','--max_mismatch', 
                        metavar='max_mismatch',
                        required=False,
                        default=None,
                        type=int, 
                        help='Max mismatch for barcode/SSI matching.')

    parser.add_argument('-b','--seq_length', 
                        metavar='seq_length',
                        required=False,
                        default=30,
                        type=int, 
                        help='Length of leading barcode to collapse [30]')

    parser.add_argument('-r','--recursion', 
                        metavar='recursion',
                        required=False,
                        default=20000,
                        type=int, 
                        help='Max recursion. Handle larger input to collapse() System default ~3000.')

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. input file base dir if not given.')   

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Combined out TSV for info. E.g. EXP.collapsed.tsv')

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')
   
    parser.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='TSV with sequence and read_count columns')
        
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    cp.read(args.config)
    
    if args.aligner is not None:
        cp.set('collapse','tool', args.aligner)

    if args.max_mismatch is not None:
        cp.set('collapse','max_mismatch', str(args.max_mismatch))    

    if args.seq_length is not None:
        cp.set('collapse','seq_length', str(args.seq_length) )
    
    cdict = format_config(cp)
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infiles={args.infile}')
    
    # set recursion
    logging.debug(f'recursionlimit = {sys.getrecursionlimit()}')
    if args.recursion is not None:
        rlimit = int(args.recursion)
        logging.info(f'set new recursionlimit={rlimit}')
        sys.setrecursionlimit(rlimit)
       
    # set outdir / outfile
    outdir = os.path.abspath('./')
    outfile = f'{outdir}/collapsed.tsv'
    if args.outdir is None:
        if args.outfile is not None:
            logging.debug(f'outdir not specified. outfile specified.')
            outfile = os.path.abspath(args.outfile)
            filepath = os.path.abspath(outfile)    
            dirname = os.path.dirname(filepath)
            filename = os.path.basename(filepath)
            (base, ext) = os.path.splitext(filename)   
            head = base.split('.')[0]
            outdir = dirname
        else:
            logging.debug(f'outdir/file not specified.')        
    else:
        # outdir specified
        outdir = os.path.abspath(args.outdir)
        if args.outfile is not None:
            logging.debug(f'outdir specified. outfile specified.')
            outfile = os.path.abspath(args.outfile)
        else:
            logging.debug(f'outdir specified. outfile not specified.')
            outfile = f'{outdir}/collapsed.tsv'

    logging.debug(f'making missing outdir: {outdir} ')
    os.makedirs(outdir, exist_ok=True)
    logging.info(f'outdir={outdir} outfile={outfile}')
    logging.debug(f'infile = {args.infile}')
    
    if args.logfile is not None:
        log = logging.getLogger()
        FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
        formatter = logging.Formatter(FORMAT)
        logStream = logging.FileHandler(filename=args.logfile)
        logStream.setFormatter(formatter)
        log.addHandler(logStream)
    
    logging.info(f'loading infile {args.infile}')
    df = load_df(args.infile)
    logging.debug(f'calling colllapse seq_length={args.seq_length} max_mismatch={args.max_mismatch}') 
    collapse_df = align_collapse_df(df, 
                      seq_length=args.seq_length, 
                      max_mismatch=args.max_mismatch, 
                      outdir=outdir, 
                      datestr=args.datestr,
                      cp=cp )
    logging.info(f'Writing collapse DF to {outfile}')
    collapse_df.to_csv(outfile, sep='\t')
    logging.info(f'Wrote DF to {outfile}')
    
         