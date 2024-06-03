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
templateswitching.m

    num_target_L2 = sum(sum(barcodematrix(:,target_L2)));          % total num of L2 molecules in L2 targets
    num_target_L1 = sum(sum(barcodematrixL1(:,target_L1)));        % total num of L1 molecules in L1 targets
    num_spikes_L1 = sum(num_spikein(target_L1));                   % total num of spike in molecules in L1 targets
    num_spikes_L2 = sum(num_spikein(target_L2));                   % total num of spike in molecules in L2 targets
    num_templateswitching = sum(sum(barcodematrix(:,target_L1)));  % num of L2 molecules detected in L1 targets
    c = num_templateswitching / (num_target_L2 * (num_spikes_L1 + num_target_L1));   % template switching coefficient

    ratio_ts = 0.5 * c * num_target_L2 + c * num_spikes_L2;         % ratio of template swtiching molecules in all L2 targets

'''

# 
# sh = StatsHandler(config, outdir=outdir, datestr=datestr)
# sh.add_value(f'/fastq/pair{pairshandled}','reads_total', num_handled)

def calc_template_switching(config, sampdf, alldf):
    '''
    calculate info about template switching rate/ratio
    '''
    logging.debug(f'calc_template_switching()')

    # total num of L2 molecules in L2 targets
    target_real = int( alldf.loc[ alldf['site'] == 'target' ].loc[ alldf['type'] == 'real' ]['umi_count'].sum() )
    
    # total num of spike in molecules in L2 targets
    target_spike = int( alldf.loc[ alldf['site'] == 'target' ].loc[ alldf['type'] == 'spike']['umi_count'].sum() )

    # total num of L1 molecules in L1 targets
    #lone_lone = alldf[alldf['site'] == 'target-lone'][alldf['type'] == 'lone'].umi_count.sum()
    lone_lone = int( alldf.loc[ alldf['site'] == 'target-lone' ].loc[ alldf['type'] == 'lone' ]['umi_count'].sum() )

    # total num of spike in molecules in L1 targets
    lone_spike = int( alldf.loc[ alldf['site'] == 'target-lone'].loc[ alldf['type'] == 'spike']['umi_count'].sum() )

    # num of L2 molecules detected in L1 targets
    lone_real = int( alldf.loc[ alldf['site'] == 'target-lone'].loc[ alldf['type'] == 'real']['umi_count'].sum() )
    logging.debug(f'target_real={target_real} target_spike={target_spike} lone_lone={lone_lone} lone_spike={lone_spike} lone_real={lone_real}')

    # c = num_templateswitching / (num_target_L2 * (num_spikes_L1 + num_target_L1));
    tsc = lone_real / (  target_real *   ( lone_spike + lone_real  ))
    logging.debug(f'tsc = {tsc}')
    
    #  0.5 * c * num_target_L2 + c * num_spikes_L2
    ratio_ts = ( 0.5 * tsc * target_real ) + ( tsc * target_spike) 
    logging.debug(f'ratio_ts = {ratio_ts}')

    try:
        sh = get_default_stats()
    except:
        sh = StatsHandler(config, outdir=outdir, datestr=None)
        
    sh.add_value(f'/qcstats/template_switching','target_real', target_real)
    sh.add_value(f'/qcstats/template_switching','target_spike', target_spike)
    sh.add_value(f'/qcstats/template_switching','lone_lone', lone_lone)
    sh.add_value(f'/qcstats/template_switching','lone_spike', lone_spike)
    sh.add_value(f'/qcstats/template_switching','lone_real', lone_real)
    sh.add_value(f'/qcstats/template_switching','ts_coeff', tsc)
    sh.add_value(f'/qcstats/template_switching','ratio', ratio_ts)



def calc_qcstats(config, sampdf,  infile, outfile=None, outdir=None):
    '''
    write stats to outfile. 
    may write other intermdiatefiles to outdir
    

    '''
    logging.debug(f'infile={infile} outfile={outfile} ')
    try:
        sh = get_default_stats()
    except:
        sh = StatsHandler(config, outdir=outdir, datestr=None)
    
    if outfile is not None:
        outdir = os.path.abspath( os.path.dirname(outfile))
        logging.debug(f'making outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)
    else:
        outfile = './qcstats.txt'

    if outdir is not None:
        os.makedirs(os.path.abspath( outdir) , exist_ok=True)

    logging.debug(f'final infile={infile} outfile={outfile} outdir={outdir} ')
    alldf = load_df(infile)
    logging.debug(f'loaded {infile} len={len(alldf)}')

    #print('targets:')
    #print(alldf[alldf.site == 'target']['type'].value_counts() )
    #print('injection:')
    #print(alldf[alldf.site == 'injection']['type'].value_counts() )
  
    calc_template_switching(config, sampdf, alldf)
    
  


    
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
        outdir = os.path.abspath( os.path.dirname(args.outfile))
        logging.debug(f'making outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)
        outfile = args.outfile
    else:
        outfile = './qcstats.txt'

    if args.outdir is not None:
        os.makedirs(args.outdir, exist_ok=True)
        outdir = args.outdir


    sh = StatsHandler(config, outdir=outdir, datestr=None, outfile='')
    calc_qcstats(config, sampdf, args.infile, outfile=outfile, outdir=outdir)
    

        
        
        
