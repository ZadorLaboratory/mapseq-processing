#!/usr/bin/env python
#
# Make shoulder plot from arbitrary table. 
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
from mapseq.calibration import *


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

    parser.add_argument('-C','--column', 
                        metavar='column',
                        required=False,
                        default='umi_count',
                        type=str, 
                        help='Column to plot [ read_count | umi_count ]')


    parser.add_argument('-p','--proportion', 
                        metavar='proportion',
                        required=False,
                        default=0.99,
                        type=float, 
                        help='Threshold for proportion of data.')

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=True,
                    default=None, 
                    type=str, 
                    help='Single XLSX report filename.')

    parser.add_argument('-L','--logfile', 
                    metavar='logfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Logfile for subprocess.')

   
    parser.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='Single TSV or Parquet MAPseq readtable with relevant column(s) ')
        
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
    outfile = os.path.abspath(args.outfile)
    filepath = os.path.abspath(outfile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    head = base.split('.')[0]
    outdir = dirname 
    os.makedirs(outdir, exist_ok=True)
    logging.info(f'handling {args.infile} to {outfile}')    
    logging.debug(f'infile = {args.infile}')
    logging.info(f'loading {args.infile}') 

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
    df = load_mapseq_df( args.infile)   
    logging.debug(f'loaded. len={len(df)} dtypes = {df.dtypes}') 
    
    rdf = calib_viruslib_blacklist(df)
    
    from matplotlib.backends.backend_pdf import PdfPages as pdfpages

    title = 'viruslib dupes rates.'
    
    page_dims = (8, 6)
    
    pdffile = os.path.join(outdir, f'{base}.pdf')
    with pdfpages(outfile) as pdfpages:
        fig, axes = plt.subplots(figsize=page_dims)
        fig.suptitle(title)
        sns.scatterplot(ax=axes, x = rdf['max_umi'],y = rdf['p_dupes'])       
        pdfpages.savefig(fig)

    logging.info(f'Plot written to {args.outfile}')