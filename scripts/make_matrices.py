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
                        help='Minimum single target UMIs for including all.')

    parser.add_argument('-T','--target_min_umi_absolute', 
                        required=False,
                        default=None,
                        type=int, 
                        help='Minimum target UMIs absolute.')    
    
    parser.add_argument('-e','--expid', 
                    metavar='expid',
                    required=False,
                    default=None,
                    type=str, 
                    help='Explicitly provided tag.')

    parser.add_argument('-l','--label', 
                    metavar='label',
                    required=False,
                    default='label',
                    type=str, 
                    help='Field to use as column label, e.g. [ label | ourtube | region ]')


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
    o = cp.read(args.config)
    if len(o) < 1:
        logging.error(f'No valid configuration. {args.config}')
        sys.exit(1)
       
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
    
    logging.info(f'loading {args.infile}') 
    df = load_mapseq_df( args.infile, fformat='vbctable', use_dask=False)
    logging.debug(f'loaded. len={len(df)} dtypes =\n{df.dtypes}') 

    if args.datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    else:
        datestr = args.datestr
    sh = StatsHandler(outdir=outdir, datestr=datestr)
    
    logging.debug(f'loaded. len={len(df)} dtypes = {df.dtypes}') 

    process_make_matrices(df,
                               exp_id = args.expid,  
                               outdir=outdir, 
                               label_column=args.label,
                               cp=cp)

    logging.info(f'Made matrices in {outdir}...')

    