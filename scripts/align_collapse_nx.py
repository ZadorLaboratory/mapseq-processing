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
from mapseq.collapse import *


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

    parser.add_argument('-R','--min_reads', 
                        required=False,
                        default=None,
                        type=int, 
                        help='Minimum reads for inclusion.')
    

    parser.add_argument('-a','--aligner', 
                    metavar='aligner',
                    required=False,
                    default=None, 
                    type=str, 
                    help='aligner tool  [bowtie | bowtie2]')

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

    parser.add_argument('-G','--group_column', 
                    metavar='group_column',
                    required=False,
                    default='brain', 
                    type=str, 
                    help='separate groups to collapse in [brain]')

    parser.add_argument('-m','--max_distance', 
                        metavar='max_distance',
                        required=False,
                        default=None,
                        type=int, 
                        help='Max mismatch/edit distance for collapse.')


    parser.add_argument('-r','--max_recursion', 
                        metavar='max_recursion',
                        required=False,
                        default=None,
                        type=int, 
                        help='Max recursion. Handle larger input to collapse() System default ~3000.')

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

    parser.add_argument('-f','--force', 
                    action="store_true", 
                    default=False, 
                    help='Recalculate even if output exists.')  

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')
   
    parser.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='Single TSV with column to be collapsed.')
        
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
    
    # set recursion
    logging.debug(f'recursionlimit = {sys.getrecursionlimit()}')
    if args.max_recursion is not None:
        rlimit = int(args.max_recursion)
        logging.info(f'(from CLI) set new recursionlimit={rlimit}')
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

    if args.group_column == 'None':
        group_column = None
    else:
        group_column = args.group_column

    
    sh = StatsHandler(outdir=outdir, datestr=datestr) 
    df = align_collapse_nx_grouped(df, 
                        column=args.column,
                        pcolumn=args.parent_column,
                        gcolumn = group_column,
                        max_distance=args.max_distance,
                        max_recursion=args.max_recursion,
                        outdir=outdir,
                        force = args.force, 
                        min_reads = args.min_reads,                       
                        cp=cp
                        )

    write_mapseq_df(df, outfile)
    logging.info('Done align_collapse.')

    logging.info(f'Making read_count frequency plots...')
    make_counts_plots(df, 
                      outdir=outdir, 
                      groupby='label', column='read_count', cp=cp )
    logging.info(f'Making read report...')    
    make_read_report_xlsx(df,
                      outdir=outdir,
                      step='collapsed',
                      cp=cp)    
    logging.info('Done.')
    