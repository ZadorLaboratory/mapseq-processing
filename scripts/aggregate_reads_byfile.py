#!/usr/bin/env python
#
# Aggregate XYZ.reads.tsv by file, 
# w/ parallel output to XYZ.aggregated.tsv 
#
# Needed to overcome Novaseq experiment size which is too
# large for global aggregation. 
#
import argparse
import logging
import os
import random
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

    parser.add_argument('-L','--logfile', 
                    metavar='logfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Logfile for subprocess.')
  

    parser.add_argument('-k', '--use_dask', 
                        action="store_true", 
                        dest='use_dask',
                        help='Use Dask subsystem.')

    parser.add_argument('-t','--dask_temp', 
                    metavar='dask_temp',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Fast, roomy storage for DASK temp files.')

    parser.add_argument('-m','--min_reads', 
                        metavar='min_reads',
                        required=False,
                        default=None,
                        type=int, 
                        help='Min reads to retain initial full read.')

    parser.add_argument('-C','--column', 
                    metavar='column',
                    required=False,
                    default=['sequence','source'], 
                    type=str, 
                    help='Read column to aggregate.')

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=True,
                    default=None, 
                    type=str, 
                    help='outdir. required for parallel operation.')

    parser.add_argument('infiles' ,
                        metavar='infiles', 
                        type=str,
                        nargs='+',
                        default=None, 
                        help='Read format TSV or Parquet files. "read" required in name. ')
        
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
    
    outdir = os.path.abspath(args.outdir)
    logging.debug(f'making missing outdir: {outdir} ')
    os.makedirs(outdir, exist_ok=True)
    logging.info(f'outdir={outdir} ')
      
    logging.debug(f'infiles = {args.infiles}')
    
    if args.logfile is not None:
        log = logging.getLogger()
        FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
        formatter = logging.Formatter(FORMAT)
        logStream = logging.FileHandler(filename=args.logfile)
        logStream.setFormatter(formatter)
        log.addHandler(logStream)

    if args.min_reads is None:
        min_reads= int( cp.get('fastq','min_reads') )
    else:
        min_reads = int(args.min_reads)

    if args.datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    else:
        datestr = args.datestr
    sh = StatsHandler(outdir=outdir, datestr=datestr) 

    for infile in args.infiles:
        logging.info(f'handling infile={infile}')
        infile = os.path.abspath(infile)
        dirpath, filename = os.path.split(infile)
        base, ext = os.path.splitext(filename)
        outbase = base.replace('reads','aggregated')
        outfile = os.path.join(outdir, f'{outbase}.tsv')
        logging.info(f'calculated outfile = {outfile}')

        df = load_mapseq_df( infile, fformat='reads', use_dask=args.use_dask, use_arrow=True)
        logging.debug(f'loaded. len={len(df)} dtypes =\n{df.dtypes}') 
        df = aggregate_reads( df, 
                            column=args.column,
                            outdir=outdir,
                            min_reads = min_reads,
                            use_dask = args.use_dask, 
                            dask_temp = args.dask_temp,
                            cp=cp 
                            )
        logging.debug(f'got aggregated df len={len(df)}:\n{df} ')

        write_mapseq_df(df, outfile)
    logging.info('Done aggregate_reads')
    