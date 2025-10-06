#!/usr/bin/env python
import argparse
import itertools
import logging
import os
import sys
from configparser import ConfigParser
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.collapse import *
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

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Dataframe of component properties.') 

    parser.add_argument('-u', '--uniques',
                        metavar='uniques',
                        required=True,
                        type=str,
                        help='Dataframe TSV with vbc_read column and same index as components.txt')
    
    parser.add_argument('-C','--components',
                        metavar='components',
                        required=True,
                        type=str,
                        help='Text file with component sets. ')

    parser.add_argument('-e','--edges',
                        metavar='edges',
                        required=True,
                        type=str,
                        help='Text file with edge pairs.')

    parser.add_argument('-s','--column', 
                    metavar='column',
                    required=False,
                    default='vbc_read', 
                    type=str, 
                    help='Column that was collapsed.') 
        
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    
    if args.config is None:
        cp = get_default_config()
    else:
        cp = ConfigParser()
        o = cp.read(args.config)
        if len(o) < 1:
            logging.error(f'No valid configuration. {args.config}')
            sys.exit(1)
       
    cdict = format_config(cp)
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'uniques={args.uniques} components={args.components}')
          
    # set outdir / outfile
    outdir = os.path.abspath('./')
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

    outdir = os.path.abspath(outdir)    
    os.makedirs(outdir, exist_ok=True)
    logging.info(f'handling {args.uniques} {args.components}  to outdir {args.outfile}')    

    if args.datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    else:
        datestr = args.datestr
         
    sh = StatsHandler(outdir=outdir, datestr=datestr)
 
    comp_df = check_components(uniques_file=args.uniques, 
                               components_file=args.components,
                               edges_file = args.edges,
                               column=args.column, 
                               cp=cp )
   
    outfile = args.outfile
    comp_df.to_csv(outfile, sep='\t')
    logging.info('Done assessing components. ')
    
    