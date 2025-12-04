#!/usr/bin/env python
import argparse
import logging
import os
import sys

from configparser import ConfigParser

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.utils import *
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

    parser.add_argument('-S','--samplesheet', 
                        metavar='samplesheet',
                        required=False,
                        default='Sample information',
                        type=str, 
                        help='XLS sheet tab name.')

    parser.add_argument('-k','--kneed_sensitivity', 
                        metavar='kneed_sensitivity',
                        required=False,
                        default=None,
                        type=float, 
                        help='Kneed S value.')

    parser.add_argument('-p','--kneed_polynomial', 
                        metavar='kneed_polynomial',
                        required=False,
                        default=None,
                        type=int, 
                        help='Kneed polynomival value.')

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=True,
                    default=None, 
                    type=str, 
                    help='Single XLSX report filename.') 
   
    parser.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='Single QPCR .XLS file.')
        
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

    qpcr_check_wells(args.infile, 
                     outfile, 
                     sensitivity=args.kneed_sensitivity, 
                     polynomial=args.kneed_polynomial,
                     cp = cp)
    
    logging.info(f'Made QPCR report in {outfile}...')