#!/usr/bin/env python
import argparse
import logging
import os
import sys

from configparser import ConfigParser

import pandas as pd

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
                        default=os.path.expanduser('~/git/mapseq-processing/etc/barcode_v2.txt'),
                        type=str, 
                        help='barcode file space separated: label   sequence')

    parser.add_argument('-s','--sampleinfo', 
                        metavar='sampleinfo',
                        required=True,
                        default=None,
                        type=str, 
                        help='XLS sampleinfo file. ')

    parser.add_argument('-m','--max_mismatch', 
                        metavar='max_mismatch',
                        required=False,
                        default=0,
                        type=int, 
                        help='Max mismatch for barcode/SSI matching.')

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. input file base dir if not given.')   

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')

    parser.add_argument('-f','--force', 
                    action="store_true", 
                    default=False, 
                    help='Recalculate even if output exists.') 

    parser.add_argument('-p', '--countsplots', 
                        action="store_true", 
                        dest='countsplots',
                        default=False, 
                        help='create countsplots for all barcodes.' )

    parser.add_argument('-r', '--readtsvs', 
                        action="store_true", 
                        dest='readtsvs',
                        default=False, 
                        help='create read tsvs.' )   
   
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

    cp = ConfigParser()
    cp.read(args.config)
    cdict = format_config(cp)
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infiles={args.infiles}')


    # check nargs. 
    if len(args.infiles) < 2:
        parser.print_help()
        print('error: the following arguments are required: 2 infiles')
        sys.exit(1)

    if len(args.infiles) % 2 != 0:
        parser.print_help()
        print('error: the following arguments are required: 2 or multiple of 2 infiles')
        sys.exit(1)    
       
    # set outdir
    outdir = None
    if args.outdir is not None:
        outdir = os.path.abspath(args.outdir)
        logging.debug(f'making missing outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)
    else:
        afile = args.infiles[0]
        filepath = os.path.abspath(afile)    
        dirname = os.path.dirname(filepath)
        outdir = dirname

    sampdf = load_sample_info(cp, args.sampleinfo)
    logging.debug(f'\n{sampdf}')
    sampdf.to_csv(f'{outdir}/sampleinfo.tsv', sep='\t')
    rtlist = get_rtlist(sampdf)

    logging.debug(f'making barcodes with label list={rtlist}')
    bcolist = load_barcodes(cp, 
                            args.barcodes, 
                            labels=rtlist, 
                            outdir=outdir, 
                            eol=True, 
                            max_mismatch=args.max_mismatch)
    logging.info(f'made list of barcode handlers, length={len(bcolist)}')
    logging.debug(bcolist)
    logging.info(f'handling {args.infiles[0]} and {args.infiles[1]} to outdir {args.outdir} with countsplots={args.countsplots}')
           
    infilelist = package_pairfiles(args.infiles)
    
    logging.debug(f'infilelist = {infilelist}')
    process_fastq_pairs(cp, sampdf, infilelist, bcolist, outdir=args.outdir, force=args.force, countsplots=args.countsplots, datestr=args.datestr)
    
    