#!/usr/bin/env python
#
# Create Python/Pandas-compatible readtable from MATLAB barseq barcode file. 
# Assumes format. First field is cellid, later fields are base codes, with lookup 
# below...
#
# 5_17_372270,1,2,4,2,2,2,1,2,2,4,1,1,1,4,2
# 5_17_380518,1,2,1,4,4,4,2,2,2,2,3,3,2,3,4
# 5_17_380548,1,2,3,1,1,3,2,1,3,1,3,2,2,3,4
#
# At minimum, needs label, type, brain columns.  
#   label = BC999 [ BC998, BC997 ]
#   type = injection
#   brain =  ( input required )
#

import argparse
import logging
import os
import random
import sys

from configparser import ConfigParser

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import *
from mapseq.barcode import *
from mapseq.utils import *
from mapseq.stats import *

BASE_MAP = { 1: 'G', 2 : 'T', 3: 'A', 4 :'C' }
BC_COL_START = 1
BC_COL_END = 15
BC_COLS = list( range(  BC_COL_START , BC_COL_END +1 ))

def read_barseq_csv( infile, cp=None):
    if cp is None:
        cp = get_default_config()
    
    logging.debug(f'BASE_MAP={BASE_MAP} BCOLS={BC_COLS}')
    
    df = pd.read_csv(infile, header=None)
    for cname in range( BC_COL_START , BC_COL_END +1):
        df[cname] = df[cname].map(BASE_MAP)
    
    newdf = pd.DataFrame()    
    newdf['cell_id'] = df[0]
    newdf['vbc_read'] = df[BC_COLS].astype(str).apply(''.join, axis=1 )
    return newdf

        
def expand_dataframe(df, cycles=5):
    '''
    Hackish, but a quick way to ensure each row of an original DF 
    has 2^<cycles> copies. 
    
    '''
    newdf = df.copy()
    for i in range(0,5):
        newdf = pd.concat( [ newdf, newdf])
    newdf.reset_index(inplace=True, drop=True)
    return newdf
  
   
def make_umi_list(listlength, strlength=12, alphabet=['A','C','G','T']):
    logging.debug(f'making str list len={listlength} from alphabet={alphabet}')
    umilist = []
    for i in range(0,listlength):
        umilist.append( ''.join(random.choices(alphabet, k=strlength)) )
    return umilist


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
    
    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=True,
                    type=str, 
                    help='Full dataset table TSV') 

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')

    parser.add_argument('-b','--brain', 
                        metavar='brain',
                        required=True,
                        type=str, 
                        help='brain id to assign.') 

    parser.add_argument('-r','--readtable', 
                        metavar='readtable',
                        required=False,
                        type=str, 
                        help='MAPseq readtable to merge with BARseq data.')

    parser.add_argument('-S','--ssi', 
                        metavar='ssi',
                        required=False,
                        default = 'BC999',
                        type=str, 
                        help='Fake rtprimer/barcode label [BC999]') 
   
    parser.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='BARseq encoded barcode CSV file.')
        
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
    logging.debug(f'outdir set to {outdir}')
    outdir = os.path.abspath(outdir)    
    os.makedirs(outdir, exist_ok=True)

    logging.info(f'handling {args.infile} to outdir {outdir}')    
    logging.debug(f'infile = {args.infile}')
        
    if args.datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    else:
        datestr = args.datestr
    sh = StatsHandler(outdir=outdir, datestr=datestr)

    logging.info(f'Opening {args.infile}')
    df = read_barseq_csv( args.infile, 
                          cp=cp)
    logging.info(f'Read CSV data. len={len(df)}')
    df = expand_dataframe(df, cycles=5)
    logging.debug(f'expanded DF len={len(df)}')
    
    df['type'] = 'injection'
    df['label'] = args.ssi
    df['brain'] = args.brain
    df['read_count'] = 999
    
    umilist = make_umi_list(len(df)) 
    umiser = pd.Series(umilist)
    df['umi'] = umiser

    if args.readtable is not None:
        logging.info(f'loading {args.readtable} for merging...')
        rdf = load_mapseq_df(args.readtable, fformat='readtable')
        df = pd.concat([rdf, df], ignore_index=True )
        df = fix_category_nans(df)
        df.fillna('',inplace=True)
        df['vbc_read_short'] = df['vbc_read'].str.slice(0,15)
        df = fix_mapseq_dtypes(df) 
    
    write_mapseq_df(df, outfile)    
    logging.info('Done make_barseq_readtable.')
        