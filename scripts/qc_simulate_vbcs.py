#!/usr/bin/env python
import argparse
import itertools
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

from mapseq.collapse import *
from mapseq.core import *
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


    parser.add_argument('-n','--n_sequences', 
                        metavar='n_sequences',
                        required=False,
                        default=10,
                        type=int, 
                        help='Number of sequences to create. ')

    parser.add_argument('-b','--n_bases', 
                        metavar='n_bases',
                        required=False,
                        default=30,
                        type=int, 
                        help='Sequence length.')

    parser.add_argument('-C','--n_cycles', 
                        metavar='n_cycles',
                        required=False,
                        default=5,
                        type=int, 
                        help='Rounds of mutation. ')

    parser.add_argument('-M','--n_copies', 
                        metavar='n_copies',
                        required=False,
                        default=5,
                        type=int, 
                        help='Number of mutated copies to make per sequence/cycle ')

    parser.add_argument('-m','--max_mismatch', 
                        metavar='max_mismatch',
                        required=False,
                        default=3,
                        type=int, 
                        help='Max positions mutated per round.')

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=os.path.abspath('./component.info.tsv'), 
                    type=str, 
                    help='E.g. component assessment tsv.') 

        
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
    logging.info(f'handling {args.n_sequences} to outdir {outdir} outfile={outfile}')    

    if args.datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    else:
        datestr = args.datestr
         
    sh = StatsHandler(outdir=outdir, datestr=datestr)

    # generate/ write parents
    (parent_list, parent_df) = generate_random( n_bases=args.n_bases, 
                                                n_sequences = args.n_sequences)
    parent_df['parent_idx'] = parent_df.index
    outf = os.path.join(outdir, 'parent_sequences.tsv')
    logging.info(f'got {len(parent_df)} sequences. Writing to {outf} ')
    parent_df.to_csv(outf, sep='\t')
    
    # generate/write mutated    
    mutated_df = generate_mutated_df(parent_df, 
                                     n_copies=args.n_copies, 
                                     n_cycles=args.n_cycles, 
                                     max_mismatch=args.max_mismatch, 
                                     alphabet = 'AGCT')
    
    outf = os.path.join(outdir, 'mutated_sequences.tsv')
    logging.info(f'got {len(mutated_df)} sequences. Writing to {outf} ')
    mutated_df.to_csv(outf, sep='\t')
    
    cdf = align_collapse_pd(mutated_df, 
                         max_mismatch=args.max_mismatch, 
                         drop=False, 
                         force=True, 
                         column='sequence', 
                         pcolumn='count')
    
    
    outf = os.path.join(outdir, 'collapsed_sequences.tsv')
    cdf.to_csv(outf, sep='\t')
    logging.info(f'Done wrote {outf} ')
    
    