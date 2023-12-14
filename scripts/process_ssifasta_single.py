#!/usr/bin/env python
#
# Single-CPU single-threaded barcode splitter
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

#from mapseq.utils import *
from mapseq.core import *
#from mapseq.barcode import *  


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
                        help='Config file to use. ')    

    parser.add_argument('-s','--sampleinfo', 
                        metavar='sampleinfo',
                        required=True,
                        default=None,
                        type=str, 
                        help='XLS sampleinfo file. ')

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=True,
                    default=None, 
                    type=str, 
                    help='Combined out TSV for all info for this FASTA, e.g. BC1.all.tsv') 

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
   
    parser.add_argument('infile' ,
                        metavar='infile', 
                        type=str,
                        nargs='?',
                        help='Barcode-specific FASTA file(s). Should be one.')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    if args.logfile is not None:
        log = logging.getLogger()
        FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
        formatter = logging.Formatter(FORMAT)
        logStream = logging.FileHandler(filename=args.logfile)
        logStream.setFormatter(formatter)
        log.addHandler(logStream)


    cp = ConfigParser()
    cp.read(args.config)
    cdict = format_config(cp)
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infile={args.infile}')
       
    # set outdir
    outdir = None
    if args.outdir is None: 
        outdir = os.path.dirname(os.path.abspath(args.outfile))
    else:
        outdir = args.outdir
    sampdf = load_sample_info(cp, args.sampleinfo)
    (rtprimer, site, brain, region) = guess_site(args.infile, sampdf)
    logging.info(f'guessed rtprimer={rtprimer} site={site} brain={brain} region={region}')
    logging.info(f'outdir={outdir} outfile={args.outfile}')
    outdf = process_ssifasta(cp, args.infile, outdir=outdir, site=site, datestr=None)
    outdf['site'] = site
    outdf['brain'] = brain
    outdf['region'] = region
    outdf.to_csv(args.outfile, sep='\t')
    logging.info(f'wrote output to {args.outfile}')
        
   
