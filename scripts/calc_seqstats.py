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


'''
'''

# 
# sh = StatsHandler(config, outdir=outdir, datestr=datestr)
# sh.add_value(f'/fastq/pair{pairshandled}','reads_total', num_handled)


def calc_seqstats(config, sampdf,  infile, outfile=None, outdir=None):
    '''
    calculate info about sequences, YY/RR spikeins.
    
    input. vs. 32 + 20 !!

                                     sequence                    tail
    0          TACTCGGAAGCGAGTGTGTCAAGGACGCTT  CCTGCTGATAGTGGTTAAAGTT
    1          GCTAACCATAGTCGGGTGCGTCAACAGCCG  CCGGTAGGAGTTGGGTCTCTCG
                        30nt                            22nt
    original regexes. 
    spikeinregex= CGTCAGTC$
    realregex = [TC][TC]$
    loneregex = [AG][AG]$
    
    modified. 
    spikeinregex = CGTCAG$ on sequence
    realregex =    ^[TC][TC] on tail. 
    loneregex =    ^[AG][AG] on tail. 
    
    '''
    logging.info(f'infile={infile} outfile={outfile} ')
    
    try:
        sh = get_default_stats()
    except:
        sh = StatsHandler(config, outdir=outdir, datestr=None)
    
    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
    else:
        outdir = './'

    if outfile is not None:
        odir = os.path.dirname(outfile)
        logging.debug(f'making outdir: {odir} ')
        os.makedirs(odir, exist_ok=True)
    else:
        outfile = f'{outdir}/seqstats.txt'    

    logging.debug(f'final infile={infile} outfile={outfile} outdir={outdir} ')
    alldf = load_df(infile)
    logging.info(f'loaded {infile} len={len(alldf)}')

    # Regexes
    sire = 'CGTCAG$'    
    #
    # what we really care about is mutated codes. not L2 or L1. 
    #
    # all spike-ins are   real = all_real - spike_ins
    #
    realre = '^[TC][TC]'
    lonere = '^[AG][AG]'
    
    logging.info(f'calculating L2s')
    # L2s (includes all reals and valid spikeins)   
    ltmap = alldf['tail'].str.contains(realre, regex=True) == True
    ltdf = alldf[ltmap]
    
    # Spike-in L2s
    logging.info(f'selecting spike-ins')
    simap = ltdf['sequence'].str.contains(sire, regex=True) == True
    sdf = ltdf[simap]
    
    # Real L2s. 
    rmap = ltdf['sequence'].str.contains(sire, regex=True) == False
    rdf = ltdf[rmap]
    
    # L1s
    logging.info(f'selecting L1s')    
    lmap = alldf['tail'].str.contains(lonere, regex=True) == True
    ldf = alldf[lmap]
    
    logging.info(f'calculating unmatched')
    remaindf = alldf[~lmap]
    remaindf = remaindf[~ltmap]
    
    
    n_total =  len(alldf)
    n_ltwo = len(ltdf)
    n_spikein = len(sdf) 
    n_real = len(rdf)
    n_lone = len(ldf)
    n_nomatch = len(remaindf)
    
    nomatch_proportion = ( n_nomatch + 1 ) /  ( n_total + 1 )
    lone_proportion = ( n_lone + 1 )  / ( n_total + 1 )

    sh.add_value(f'/seqstats','n_total', n_total)    
    sh.add_value(f'/seqstats','n_ltwo', n_lone )
    sh.add_value(f'/seqstats','n_spikein', n_spikein )
    sh.add_value(f'/seqstats','n_real', n_real )
    sh.add_value(f'/seqstats','n_lone', n_lone )
    sh.add_value(f'/seqstats','n_nomatch', n_nomatch )
    sh.add_value(f'/seqstats','nomatch_proportion', nomatch_proportion )
    sh.add_value(f'/seqstats','lone_proportion', lone_proportion  )
    

    
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
                    help='Text file output.  "seqstats.txt" if not given')  

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default='./', 
                    type=str, 
                    help='outdir. input file base dir if not given.')     

    parser.add_argument('infile',
                        metavar='infile',
                        help='MXXX.all.fulldf.tsv file from p process_fastq step. ')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    config = ConfigParser()
    config.read(args.config)
    cdict = {section: dict(config[section]) for section in config.sections()}
    
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infile={args.infile}')
    
    sampdf = load_sample_info(config, args.sampleinfo)
    logging.debug(f'\n{sampdf}')
    sampdf.to_csv('./sampinfo.tsv', sep='\t')
      
    outdir = None
    outfile = None
    if args.outfile is not None:
        outdir = os.path.dirname(args.outfile)
        logging.debug(f'making outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)
        outfile = args.outfile
    else:
        outfile = './seqstats.txt'

    if args.outdir is not None:
        os.makedirs(args.outdir, exist_ok=True)
        outdir = args.outdir


    sh = StatsHandler(config, outdir=outdir, datestr=None)
    calc_seqstats(config, sampdf, args.infile, outfile=outfile, outdir=outdir)
    

        
        
        
