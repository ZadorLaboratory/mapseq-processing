#!/usr/bin/env python
#
#  reads EXP.all.tsv table and calculates various stats.    
#

import argparse
import logging
import os
import sys
import traceback

from configparser import ConfigParser

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.utils import *
from mapseq.core import *


def calc_qcstats(config, infile, outfile):
    '''
    
    '''
    logging.debug(f'infile={infile} outfile={outfile} ')
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
        logging.debug(f'made outdir={outdir}')
    alldf = load_df(infile)
    logging.debug(f'loaded {infile} len={len(alldf)}')

    print('targets:')
    print(alldf[alldf.site == 'target']['type'].value_counts() )
    print('injection:')
    print(alldf[alldf.site == 'injection']['type'].value_counts() )
    


    


    
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
                        help='config file.')  
    
    parser.add_argument('-s','--sampleinfo', 
                        metavar='sampleinfo',
                        required=True,
                        default=None,
                        type=str, 
                        help='XLS sampleinfo file or sampleinfo.tsv. ')  
    
    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                   required=False,
                    default=None, 
                    type=str, 
                    help='Text file output.  "qcstats.txt" if not given')  

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default='./', 
                    type=str, 
                    help='outdir. input file base dir if not given.')     

    parser.add_argument('infile',
                        metavar='infile',
                        help='MXXX.all.tsv file')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    cp.read(args.config)
    cdict = {section: dict(cp[section]) for section in cp.sections()}
    
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infile={args.infile}')
    
    sampdf = load_sample_info(cp, args.sampleinfo)
    logging.debug(f'\n{sampdf}')
      
    outdir = None
    if args.outfile is not None:
        outdir = os.path.dirname(args.outfile)
        logging.debug(f'making outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)
        outfile = args.outfile
    else:
        outfile = './qcstats.txt'

    calc_qcstats(cp, args.infile, outfile = outfile)
    

        
        
        
      