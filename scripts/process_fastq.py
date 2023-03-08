#!/usr/bin/env python
import argparse
import logging
import os
import sys

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)
gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from configparser import ConfigParser

import pandas as pd

from cshlwork.utils import JobRunner, JobStack, JobSet
from mapseq.core import load_sample_info, load_barcodes, process_fastq_pair, make_summaries  


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
                        help='barcode file space separated.')

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

    parser.add_argument('-o','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. cwd if not given.') 

    parser.add_argument('-f','--force', 
                    action="store_true", 
                    default=False, 
                    help='Recalculate even if output exists.') 

    parser.add_argument('infiles' ,
                        metavar='infiles', 
                        type=str,
                        nargs=2,
                        default=None, 
                        help='Read1 and Read2 fastq files')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    cp.read(args.config)
    cdict = {section: dict(cp[section]) for section in cp.sections()}
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infiles={args.infiles}')
    
    sampdf = load_sample_info(cp,args.sampleinfo)
    logging.debug(sampdf)
    rtlist = list(sampdf.rtprimer.dropna())
        
    bcolist = load_barcodes(cp, 
                            args.barcodes, 
                            labels=rtlist, 
                            outdir=args.outdir, 
                            eol=True, 
                            max_mismatch=args.max_mismatch)
    logging.debug(bcolist)
    process_fastq_pair(cp, args.infiles[0], args.infiles[1], bcolist, outdir=args.outdir, force=args.force)
    make_summaries(cp, bcolist)
    
    
    
    
    