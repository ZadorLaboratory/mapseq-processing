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

    parser.add_argument('-s','--sampleinfo', 
                        metavar='sampleinfo',
                        required=False,
                        default=None,
                        type=str, 
                        help='Optional XLS sampleinfo file for per-sample changes.')

    parser.add_argument('-S','--samplesheet', 
                        metavar='samplesheet',
                        required=False,
                        default='Sample information',
                        type=str, 
                        help='XLS sheet tab name.')

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. input file base dir if not given.')   
   
    parser.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='Single TSV or parquet read-oriented file')
        
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

    if args.sampleinfo is not None:
        logging.debug(f'loading sample DF...')
        sampdf = load_sample_info(args.sampleinfo, args.samplesheet, cp)
        logging.debug(f'\n{sampdf}')
        sampdf.to_csv(f'{outdir}/sampleinfo.tsv', sep='\t')
    else:
        logging.info(f'No sampleinfo given. BCXXX column labels. No per-sample spikeins.')
        sampdf = None    
       
    logging.info(f'loading {args.infile}') 
    df = load_mapseq_df( args.infile, fformat='reads', use_dask=False)
    logging.debug(f'loaded. len={len(df)} dtypes =\n{df.dtypes}') 
    
    logging.info(f'dropping sequence column...')
    df = df.drop('sequence', axis=1)

    sldf = qc_make_readmatrix( df,
                        sampdf=sampdf,
                        outdir=outdir, 
                        cp=cp)
    outfile = os.path.join(outdir, 'EXP.source_reads.tsv')
    sldf.to_csv(outfile, sep='\t')
        
    logging.info(f'Made QC read matrix in {outfile}...')

    