#!/usr/bin/env python
#
#   merges and analyzes per-barcode dataframe files. 
#   outputs normalized barcode matrix. 
#
#

import argparse
import logging
import os
import sys
import traceback

from configparser import ConfigParser

import pandas as pd

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)
gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)


from cshlwork.utils import write_config, merge_tsvs
from mapseq.core import *

    
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

    parser.add_argument('-C', '--combined', 
                        action="store_true",
                        default=True, 
                        dest='combined', 
                        help='make combined multi-plot PDF file. ')
    
    parser.add_argument('-c','--config', 
                        metavar='config',
                        required=False,
                        default=os.path.expanduser('~/git/mapseq-processing/etc/mapseq.conf'),
                        type=str, 
                        help='config file.')    

    parser.add_argument('-r','--recursion', 
                        metavar='recursion',
                        required=False,
                        default=100000,
                        type=int, 
                        help='Max recursion. Handle larger input to collapse() Default is ~3000.')

    parser.add_argument('-e','--expid', 
                    metavar='expid',
                    required=False,
                    default='MAPSeq Experiment',
                    type=str, 
                    help='Explicitly provided experiment id, e.g. M205')

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. input file base dir if not given.')     

    parser.add_argument('-s','--sampleinfo', 
                        metavar='sampleinfo',
                        required=True,
                        default=None,
                        type=str, 
                        help='XLS sampleinfo file. ')

    parser.add_argument('infiles',
                        metavar='infiles',
                        nargs ="+",
                        type=str,
                        help='"all" TSV from process_ssifasta. columns=(sequence, counts, type, label, brain)')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    cp.read(args.config)
    cdict = {section: dict(cp[section]) for section in cp.sections()}
    
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infiles={args.infiles}')
    
    if args.recursion is not None:
        sys.setrecursionlimit(int(args.recursion))
    
        
    outdir = None
    if args.outdir is not None:
        outdir = os.path.abspath(args.outdir)
        logging.debug(f'making missing outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)
    else:
        afile = args.infiles[0]
        filepath = os.path.abspath(afile)    
        dirname = os.path.dirname(filepath)
        outdir = dirname

    if args.outdir is not None:
        cfilename = f'{args.outdir}/process_merged.config.txt'
    else:
        afile = args.infiles[0]
        filepath = os.path.abspath(afile)    
        dirname = os.path.dirname(filepath)
        cfilename = f'{dirname}/process_merged.config.txt'
    
    write_config(cp, cfilename, timestamp=True)        
    
    sampdf = load_sample_info(cp, args.sampleinfo)
    logging.debug(f'\n{sampdf}')
        
    # create and handle 'real' 'spikein' and 'normalized' barcode matrices...
    process_merged(cp, args.infiles, outdir, expid=args.expid, recursion = args.recursion, combined_pdf=args.combined )
    
    
    
              
    