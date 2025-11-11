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

    parser.add_argument('-b','--barcodes', 
                        metavar='barcodes',
                        required=False,
                        default=None,
                        type=str, 
                        help='barcode file space separated: label sequence')

    parser.add_argument('-s','--sampleinfo', 
                        metavar='sampleinfo',
                        required=True,
                        default=None,
                        type=str, 
                        help='XLS or TSV sampleinfo file. ')

    parser.add_argument('-S','--samplesheet', 
                        metavar='samplesheet',
                        required=False,
                        default='Sample information',
                        type=str, 
                        help='XLS sheet tab name.')

    parser.add_argument('-t','--truncate', 
                        metavar='truncate',
                        required=False,
                        default=None,
                        type=int, 
                        help='Truncate vbc_read field. ')

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
                    help='Full dataset table TSV') 

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')
   
    parser.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='Single split/filtered TSV or parquet file.')
        
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    if args.logfile is not None:
        log = logging.getLogger()
        FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
        formatter = logging.Formatter(FORMAT)
        logStream = logging.FileHandler(filename=args.logfile)
        logStream.setFormatter(formatter)
        log.addHandler(logStream)

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
    outfile = f'{outdir}/readtable.tsv'
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
            outfile = f'{outdir}/readtable.tsv'

    outdir = os.path.abspath(outdir)    
    os.makedirs(outdir, exist_ok=True)

    logging.info(f'handling {args.infile} to outdir {outdir}')    
    logging.debug(f'infile = {args.infile}')
        
    logging.debug(f'loading sample DF...')
    sampdf = load_sample_info(args.sampleinfo, args.samplesheet, cp)
    logging.debug(f'\n{sampdf}')
    sampdf.to_csv(f'{outdir}/sampleinfo.tsv', sep='\t')
    
    logging.info(f'loading {args.infile}') 
    df = load_mapseq_df( args.infile, fformat='filtered', use_dask=False)
    logging.debug(f'loaded. len={len(df)} dtypes =\n{df.dtypes}') 
 
    if args.datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    else:
        datestr = args.datestr
    sh = StatsHandler(outdir=outdir, datestr=datestr)
        
    logging.debug(f'loaded. len={len(df)} dtypes = {df.dtypes}') 
    df = process_make_readtable_pd(df,
                                   sampdf, 
                                   args.barcodes, 
                                   outdir=outdir, 
                                   cp=cp)

    write_mapseq_df(df, outfile)
      
    logging.info(f'making read_count frequency plots...')
    make_counts_plots(df, 
                      outdir=outdir, 
                      groupby='label', 
                      column='read_count',
                      min_count = 1,
                      nranks=1000000, 
                      cp=cp )

    make_counts_plots(df, 
                      outdir=outdir, 
                      groupby='label', 
                      column='read_count',
                      min_count = 2, 
                      nranks=1000000,
                      cp=cp )    

    logging.info(f'making read_count frequency plots...')
    make_counts_plots(df, 
                      outdir=outdir, 
                      groupby='label', 
                      column='read_count',
                      min_count = 1,
                      nranks=None, 
                      cp=cp )

    make_counts_plots(df, 
                      outdir=outdir, 
                      groupby='label', 
                      column='read_count',
                      min_count = 2,
                      nranks=None, 
                      cp=cp )    



    
    logging.info(f'making overall read_count frequency plot...')        
    make_freqplot_single_sns(df, 
                           title='Overall read count frequency',  
                           outfile=os.path.join(outdir, 'readtable-frequency-plot.pdf'),
                           column='read_count',
                           scale='log10' )
    logging.info('Done.')

    
    make_read_report_xlsx(df,
                          outdir=outdir,
                          step='readtable',
                          cp=cp)

    if args.truncate is not None:
        N_TRUNC = int(args.truncate)
        logging.info(f'truncating vbc_read to {N_TRUNC}')
        df['vbc_read_full'] = df['vbc_read']
        df['vbc_read'] = df['vbc_read'].str.slice(0, N_TRUNC)
        dirname, base, ext = split_path(outfile)
        outfile = os.path.join(dirname, f'{base}.{N_TRUNC}.tsv')
        write_mapseq_df(df, outfile)

    
    logging.info('Done make_readtable.')
        