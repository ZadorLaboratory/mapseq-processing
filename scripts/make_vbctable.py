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

    parser.add_argument('-i','--inj_min_reads', 
                        required=False,
                        default=None,
                        type=int, 
                        help='Minimum injection reads for inclusion.')
    
    parser.add_argument('-t','--target_min_reads', 
                        required=False,
                        default=None,
                        type=int, 
                        help='Minimum target reads for inclusion.')

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
                        help='Single table TSV or Parquet file of all read information.')
        
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
    
    if args.inj_min_reads is None:
        inj_min_reads = int(cp.get('vbctable','inj_min_reads'))
    else:
        inj_min_reads = args.inj_min_reads

    if args.target_min_reads is None:
        target_min_reads = int(cp.get('vbctable','target_min_reads'))
    else:
        target_min_reads = args.target_min_reads

     
    if args.infile.endswith('.tsv'):
        logging.info(f'loading {args.infile} as tsv')
        df = load_readtable(args.infile) 
    elif args.infile.endswith('.parquet'):
        logging.info(f'loading {args.infile} as parquet')
        df = pd.read_parquet(args.infile)
    else:
        logging.error('input file must have relevant extension .tsv or .parquet')
        sys.exit(1)
    
    logging.debug(f'loaded. len={len(df)} dtypes = {df.dtypes}') 
       
    # Make final VBC/UMI based table (each row is a neuron)
    logging.debug(f'args={args}')
    df = process_make_vbctable_pd(df,
                               outdir=outdir,
                               inj_min_reads = inj_min_reads,
                               target_min_reads = target_min_reads, 
                               datestr=args.datestr,
                               cp=cp)

    logging.debug(f'inbound df len={len(df)} columns={df.columns}')
    logging.info(f'Got dataframe len={len(df)} Writing to {outfile}')
    df.to_csv(outfile, sep='\t')
    
    dir, base, ext = split_path(outfile)
    outfile = os.path.join( dir , f'{base}.parquet')
    logging.info(f'df len={len(df)} as parquet to {outfile}...')
    df.to_parquet(outfile)