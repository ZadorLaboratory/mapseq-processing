#!/usr/bin/env python
import argparse
import logging
import os
import sys

from configparser import ConfigParser

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import json


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

    parser.add_argument('-C','--column', 
                    metavar='column',
                    required=False,
                    default='vbc_read', 
                    type=str, 
                    help='column to collapse [vbc_read]')
    
    parser.add_argument('-P','--parent_column', 
                    metavar='parent_column',
                    required=False,
                    default='read_count', 
                    type=str, 
                    help='how to choose sequence [read_count]')    

    parser.add_argument('-r','--recursion', 
                        metavar='recursion',
                        required=False,
                        default=20000,
                        type=int, 
                        help='Max recursion. Handle larger input to collapse() System default ~3000.')

    parser.add_argument('-M','--mcomponents', 
                        metavar='mcomponents',
                        required=True,
                        type=str, 
                        help='multi-components.json file')

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
                    help='Combined read, read_count TSV') 

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')
   
    parser.add_argument('fulldf',
                        metavar='fulldf',
                        type=str,
                        help='Single standard reads TSV/Parquet from process_fastq_pairs ')

    parser.add_argument('uniquedf',
                        metavar='uniquedf',
                        type=str,
                        help='Single unique TSV ')
        
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    cp.read(args.config)
       
    cdict = format_config(cp)
    logging.debug(f'Running with config. {args.config}: {cdict}')
    
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
            logging.debug(f'outdir set to {outdir}')
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
      
    logging.debug(f'fulldf = {args.fulldf}')
    logging.debug(f'uniquedf = {args.uniquedf}')
    logging.debug(f'mcomponents = {args.mcomponents}')

    
    if args.logfile is not None:
        log = logging.getLogger()
        FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
        formatter = logging.Formatter(FORMAT)
        logStream = logging.FileHandler(filename=args.logfile)
        logStream.setFormatter(formatter)
        log.addHandler(logStream)
    
    logging.info(f'loading {args.fulldf}') 
    if args.fulldf.endswith('.tsv'):
        logging.info(f'loading {args.fulldf} as tsv')
        fulldf = load_readstsv(args.fulldf) 
    elif args.fulldf.endswith('.parquet'):
        logging.info(f'loading {args.fulldf} as parquet')
        fulldf = pd.read_parquet(args.fulldf)
    else:
        logging.error('input file must have relevant extension .tsv or .parquet')
        sys.exit(1)
    logging.debug(f'loaded fulldf. len={len(fulldf)} dtypes = {fulldf.dtypes}')

    logging.info(f'loading {args.uniquedf}') 
    if args.uniquedf.endswith('.tsv'):
        logging.info(f'loading {args.uniquedf} as tsv')
        uniquedf = load_df(args.uniquedf) 
    elif args.uniquedf.endswith('.parquet'):
        logging.info(f'loading {args.uniquedf} as parquet')
        uniquedf = pd.read_parquet(args.uniquedf)
    else:
        logging.error('input file must have relevant extension .tsv or .parquet')
        sys.exit(1)
    logging.debug(f'loaded uniquedf. len={len(uniquedf)} dtypes = {uniquedf.dtypes}')


    logging.info(f'loading {args.mcomponents}')
    mcomponents = None
    with open(args.mcomponents, 'r') as fh:
        mcomponents = json.load(fh)
    logging.debug(f'mcomponents type = {type(mcomponents)}')

    if args.datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    else:
        datestr = args.datestr
    sh = StatsHandler(outdir=outdir, datestr=datestr) 
 
    
    
    df = collapse_only_pd(fulldf,
                     uniquedf,
                     mcomponents,   
                     column=args.column,
                     pcolumn=args.parent_column,
                     outdir=outdir, 
                     cp=cp)
    
    logging.info(f'Saving len={len(df)} as TSV to {outfile}...')
    df.to_csv(outfile, sep='\t')
    
    dir, base, ext = split_path(outfile)
    outfile = os.path.join(dir, f'{base}.parquet')
    logging.info(f'df len={len(df)} as parquet to {outfile}...')
    df.to_parquet(outfile)
    