#!/usr/bin/env python
import argparse
import logging
import os
import sys

from configparser import ConfigParser

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

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

    parser.add_argument('-i','--inj_min_umi', 
                        required=False,
                        default=None,
                        type=int, 
                        help='Minimum injection UMIs for inclusion.')
    
    parser.add_argument('-t','--target_min_umi', 
                        required=False,
                        default=None,
                        type=int, 
                        help='Minimum target UMIs for inclusion.')

    parser.add_argument('-e','--expid', 
                    metavar='expid',
                    required=False,
                    default='M001',
                    type=str, 
                    help='Explicitly provided experiment id, e.g. M205')

    parser.add_argument('-L','--logfile', 
                    metavar='logfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Logfile for subprocess.')


    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. input file base dir if not given.')   

#    parser.add_argument('-o','--outfile', 
#                    metavar='outfile',
#                    required=False,
#                    default=None, 
#                    type=str, 
#                    help='Full dataset table TSV') 

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')
   
    parser.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='Single TSV or parquet file')
        
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    cp.read(args.config)
       
    cdict = format_config(cp)
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infiles={args.infile}')
          
    # set outdir / outfile
    if args.outdir is None:
        outdir = os.path.abspath('./')
    else:
        outdir = os.path.abspath(args.outdir)
 
    os.makedirs(outdir, exist_ok=True)
    logging.info(f'handling {args.infile} to outdir {outdir}')    
    logging.debug(f'infile = {args.infile}')
    
    if args.logfile is not None:
        log = logging.getLogger()
        FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
        formatter = logging.Formatter(FORMAT)
        logStream = logging.FileHandler(filename=args.logfile)
        logStream.setFormatter(formatter)
        log.addHandler(logStream)    
    
    if args.infile.endswith('.tsv'):
        logging.info(f'loading {args.infile} as tsv')
        df = load_readtable(args.infile) 
    elif args.infile.endswith('.parquet'):
        logging.info(f'loading {args.infile} as parquet')
        df = pd.read_parquet(args.infile)
    else:
        logging.error('input file must have relevant extension .tsv or .parquet')
        sys.exit(1)
    
    logging.debug(f'loaded. len={len(df)} dtypes = {df.dtypes}') 
    process_make_matrices_pd(df,
                                   exp_id = args.expid,  
                                   inj_min_umi = args.inj_min_umi,
                                   target_min_umi = args.target_min_umi,
                                   outdir=outdir, 
                                   datestr=args.datestr,
                                   label_column='label',
                                   cp=cp)
    
    logging.info(f'Made matrices in {outdir}...')

    