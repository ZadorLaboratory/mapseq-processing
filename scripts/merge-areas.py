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

from cshlwork.utils import write_config, merge_tsvs

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

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
    
    parser.add_argument('-c','--config', 
                        metavar='config',
                        required=False,
                        default=os.path.expanduser('~/git/mapseq-processing/etc/mapseq.conf'),
                        type=str, 
                        help='config file.')    
    
    parser.add_argument('-o','--outprefix', 
                    metavar='outprefix',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outfile prefix. stdout if not given.')  

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
    cdict = {section: dict(cp[section]) for section in cp.sections()}
    
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infiles={args.infiles}')
        
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
        cfilename = f'{args.outdir}/merge_areas.config.txt'
    else:
        afile = args.infiles[0]
        filepath = os.path.abspath(afile)    
        dirname = os.path.dirname(filepath)
        cfilename = f'{dirname}/merge_areas.config.txt'
    
    write_config(cp, cfilename, timestamp=True)        
    (bcmdf, sbcmdf) = process_merge_areas(cp, args.infiles, outdir)
    if args.outprefix is None:
        print(bcmdf)
        print(sbcmdf)
    else:
        bcmdf.to_csv(f'{args.outprefix}.bcm.tsv', sep='\t')
        sbcmdf.to_csv(f'{args.outprefix}.sbcm.tsv', sep='\t')    
    
          
    