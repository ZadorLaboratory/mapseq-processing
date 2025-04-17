#!/usr/bin/env python
#
# Make shoulder plots from READTABLE OR AGGREGATED
#
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
from mapseq.plotting import *


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

    parser.add_argument('-b','--barcodes', 
                        metavar='barcodes',
                        required=False,
                        default=os.path.expanduser('~/git/mapseq-processing/etc/barcode_v2.txt'),
                        type=str, 
                        help='barcode file space separated: label sequence')

    parser.add_argument('-f','--format', 
                        metavar='format',
                        required=False,
                        default='readtable',
                        type=str, 
                        help='input dataframe type [ aggregated | readtable ]')

    #parser.add_argument('-s','--sampleinfo', 
    #                    metavar='sampleinfo',
    #                    required=True,
    #                    default=None,
    #                    type=str, 
    #                    help='XLS sampleinfo file. ')

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
                        help='Single TSV or Parquet MAPseq readtable (with read_count and site columns')
        
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
    if args.outdir is not None:
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

    if args.datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    else:
        datestr = args.datestr
    sh = StatsHandler(outdir=outdir, datestr=datestr)

    #logging.debug(f'loading sample DF...')
    #sampdf = load_sample_info(cp, args.sampleinfo)
    #logging.debug(f'\n{sampdf}')
    #sampdf.to_csv(f'{outdir}/sampleinfo.tsv', sep='\t')

    logging.info(f'loading {args.infile}') 
    #if args.format == 'aggregated':
    #    df = load_mapseq_df( args.infile, fformat='aggregated', use_dask=False)
    #    logging.debug(f'loaded. len={len(df)} dtypes =\n{df.dtypes}') 
    #    df = set_siteinfo(df, sampdf, cp=cp)
    #    logging.debug(f'set siteinfo on dataframe: {df}')
    
    df = load_mapseq_df( args.infile, fformat=args.format, use_dask=False)   
    logging.debug(f'loaded. len={len(df)} dtypes = {df.dtypes}') 
    make_counts_plots(df, outdir=outdir, groupby='site', column='read_count', cp=cp )
    
    logging.info(f'Plots written to {outdir}')
   