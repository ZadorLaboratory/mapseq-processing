#!/usr/bin/env python
#
# Consume vbctable, sampleinfo, and control.real TSVs 
# Check for control VBC presence in injection sites in vbctable.  
#
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

    parser.add_argument('-C','--controls', 
                        metavar='controls',
                        required=True,
                        default=None,
                        type=str, 
                        help='Mapseq controls vbctable.')

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=True,
                    default=None, 
                    type=str, 
                    help='Input matrix modified with added columns for control sites. ')   
  
    parser.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='Real barcode matrix (rbcm) DF, with VBC sequence index.')
        
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
    outdir = os.path.abspath('./')
    outfile = f'{outdir}/readtable.tsv'          
    if args.outfile is not None:
        logging.debug(f'outfile specified.')
        outfile = os.path.abspath(args.outfile)
        filepath = os.path.abspath(outfile)    
        dirname = os.path.dirname(filepath)
        filename = os.path.basename(filepath)
        (base, ext) = os.path.splitext(filename)   
        head = base.split('.')[0]
        outdir = dirname
        logging.debug(f'outdir set to {outdir}')
    else:
        logging.debug(f'outdir/file not specified.')      
    os.makedirs(outdir, exist_ok=True)
    logging.info(f'handling {args.infile} to outdir {outfile}')    
       
    logging.info(f'loading controls DF {args.controls} ')
    cdf = load_mapseq_df(args.controls, fformat='vbctable')
    logging.debug(f'loaded. len={len(cdf)} dtypes =\n{cdf.dtypes}')

    logging.info('loading matrix file...')
    mdf = load_mapseq_matrix_df( args.infile)
    logging.debug(f'loaded. len={len(mdf)} dtypes =\n{mdf.dtypes}') 
    project_id = cp.get('project','project_id') 

    outdf = qc_matrix_vbc_overlap(mdf, cdf, cp=cp )
    logging.debug(f'got outdf: \n{outdf}')      
    outdf.to_csv(outfile, sep='\t')    
    logging.info(f'Wrote new matrix file: {outfile} ')




   