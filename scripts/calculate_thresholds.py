#!/usr/bin/env python
#
#  creates read counts plot/plots to determine thresholds.   
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

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from cshlwork.utils import *

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import *
    
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

    parser.add_argument('-f','--fraction', 
                    metavar='fraction',
                    required=False,
                    default=None,
                    type=float, 
                    help='threshold that captures fraction of counts')
    
    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='PDF plot out file. "thresholds.tsv" if not given')  

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. input file base dir if not given.')     

    parser.add_argument('infiles',
                        metavar='infiles',
                        nargs ="+",
                        type=str,
                        help='BCXXX.44.counts.tsv files')
       
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
    
    sampdf = load_sample_info(cp, args.sampleinfo)
    logging.debug(f'\n{sampdf}')
      
    outdir = None
    if args.outdir is not None:
        outdir = os.path.abspath(args.outdir)
        logging.debug(f'making outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)
    else:
        afile = args.infiles[0]
        filepath = os.path.abspath(afile)    
        dirname = os.path.dirname(filepath)
        outdir = dirname

    #make_countsplots(cp, args.infiles)
    threshlist = calculate_thresholds_all(cp, sampdf, args.infiles, outfile=args.outfile, fraction=args.fraction)
    for item in threshlist:
        print(item)
      