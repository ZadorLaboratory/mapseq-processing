#!/usr/bin/env python
#
#  Make Pandas Dataframe from paired-end sequence FASTQ files. 
#
#
import argparse
import logging
import os
import sys

from configparser import ConfigParser

import pandas as pd

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.utils import *


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
                        default=os.path.expanduser('~/git/mapseq-processing/etc/utils.conf'),
                        type=str, 
                        help='configuration file.')

    parser.add_argument('-C', '--noconcat', 
                        action="store_true", 
                        dest='noconcat',
                        default=False, 
                        help='Do not concat Read1 + Read2')
    
    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=True,
                    type=str, 
                    help='Combined read, read_count TSV')   
   
    parser.add_argument('infiles' ,
                        metavar='infiles', 
                        type=str,
                        nargs='+',
                        default=None, 
                        help='Read1 and Read2 [Read1B  Read2B ... ] fastq files')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    # check nargs. 
    if (len(args.infiles) < 2)  or (len(args.infiles) % 2 != 0 ):
        parser.print_help()
        print('error: the following arguments are required: 2 or multiple of 2 infiles')
        sys.exit(1)

    cp = ConfigParser()
    o = cp.read(args.config)
    if len(o) < 1:
        logging.error(f'No valid configuration. {args.config}')
        sys.exit(1)
    cdict = format_config(cp)
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infiles={args.infiles}')
    
    outfile = os.path.abspath(args.outfile)
    filepath = os.path.abspath(outfile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    head = base.split('.')[0]
    outdir = dirname

    logging.debug(f'making missing outdir: {outdir} ')
    os.makedirs(outdir, exist_ok=True)
    logging.info(f'outdir={outdir} outfile={outfile}')

    logging.info(f'handling {args.infiles[0]} and {args.infiles[1]} to outdir {args.outfile}')
    infilelist = package_pairs(args.infiles)   
    
    logging.debug(f'infilelist = {infilelist}')
           
    df = read_fastq_pairs(infilelist, noconcat=args.noconcat, cp=cp)
    logging.info(f'writing TSV to {outfile}')
    df.to_csv(outfile, sep='\t')
    
    of = os.path.join(dirname, f'{base}.fasta')
    logging.info(f'writing FASTA to {of}')
    write_fasta_from_df(df, of)
    
    logging.info('Done extract_fastq_pairs')  
    