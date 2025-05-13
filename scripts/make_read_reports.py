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

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Read table TSV aggregated to viral barcode TSV.') 

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')
   
    parser.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='Single readtable TSV or Parquet.')
        
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
    outdir = os.path.abspath('./')
    outfile = f'{outdir}/read.table.tsv'
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
            outfile = f'{outdir}/read.table.tsv'

    outdir = os.path.abspath(outdir)    
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
    df = load_mapseq_df( args.infile, fformat='readtable', use_dask=False)
    logging.debug(f'loaded. len={len(df)} dtypes =\n{df.dtypes}') 
       
    if args.datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    else:
        datestr = args.datestr


    logging.info(f'making read_count frequency plots...')
    make_counts_plots(df, 
                      outdir=outdir, 
                      groupby='label', column='read_count', cp=cp )
    
    make_read_report_xlsx(df,
                      outdir=outdir,
                      cp=cp)
    
    logging.info('Done make_samplereports.')