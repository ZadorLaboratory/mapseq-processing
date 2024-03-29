#!/usr/bin/env python
#
#  merge paired-end fastq files. 
#  optionally move 3' barcode to beginning for processing by ultraplex. 
# 
import argparse
import logging
import os
import sys

from configparser import ConfigParser

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import *
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
                        default=os.path.expanduser('~/git/mapseq-processing/etc/mapseq.conf'),
                        type=str, 
                        help='out file.')    


    parser.add_argument('-n','--barcodelength', 
                        metavar='barcodelength',
                        required=False,
                        default=0,
                        type=int, 
                        help='Length of 3-prime barcode to pre-pend')

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. input file base dir if not given.')   

    parser.add_argument('-f','--force', 
                    action="store_true", 
                    default=False, 
                    help='Recalculate even if output exists.') 
   
    parser.add_argument('-t','--threads', 
                        metavar='threads',
                        required=False,
                        default=1,
                        type=int, 
                        help='Perform SSI matching in separate processes, one per SSI. Default 1. Negative numbers mean use all but that number.')    

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
    cdict = {section: dict(cp[section]) for section in cp.sections()}
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
    afile = args.infiles[0]
    filepath = os.path.abspath(afile)    
    outdir = os.path.dirname(filepath)
    if args.outdir is not None:
        outdir = args.outdir    

    cfilename = f'{outdir}/merge_ends.config.txt'
    write_config(cp, cfilename, timestamp=True)
    
    infilelist = package_pairfiles(args.infiles)
    
    logging.debug(f'infilelist = {infilelist}')
    
    merge_fastq_pairs(cp, infilelist, args.outdir)
