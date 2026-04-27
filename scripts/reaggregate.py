#!/usr/bin/env python
#
# Take in one or more aggregated dataframe files
# Re-aggregate by sequence
# Output SINGLE fully aggregated dataframe. 
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

    parser.add_argument('-C','--column', 
                    metavar='column',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Column(s) to aggregate on.')

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=True, 
                    type=str, 
                    help='Full dataset table TSV') 

    parser.add_argument('infiles' ,
                        metavar='infiles', 
                        type=str,
                        nargs='+',
                        default=None, 
                        help='Read format TSV or Parquet files.')
        
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
    
    outpath, of = os.path.split(args.outfile)
    outdir = os.path.abspath(outpath)
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

    if args.datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    else:
        datestr = args.datestr
    sh = StatsHandler(outdir=outdir, datestr=datestr) 

    columns = None
    if args.column is not None:
        columns = [ x.strip() for x in args.columns.split(',') ] 
    else:
        columns = ['sequence']

    df = reaggregate_reads( args.infiles, 
                            column=args.column,
                            use_dask = args.use_dask, 
                            dask_temp = args.dask_temp,
                            cp=cp )
    logging.info(f'Got reaggregated DF len={len(df)}')

    write_mapseq_df(df, args.outfile)
    logging.info('Done aggregate_reads')
    