#!/usr/bin/env python
#
# 
import argparse
import json
import logging
import os
import sys
import traceback

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import *
from mapseq.utils import *



def merge_mapseq_dataframes(infiles, outfile):
    '''

     
    '''
    mergelist = []
    for infile in infiles:
        logging.debug(f'reading {infile}')
        df = load_mapseq_df(infile)
        mergelist.append(df)
    logging.debug(f'created list of {len(mergelist)} dataframes.')
    log_objectinfo(mergelist, label='mergelist')

    outdf = pd.concat(mergelist, ignore_index =True, axis=0 )
    logging.info(f'created merged output datframe {outdf}, resetting index...')
    outdf.reset_index(inplace=True, drop=True)
    return outdf


 
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
       
    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=True,
                    type=str, 
                    help='merged stats json file.')  

    parser.add_argument('infiles' ,
                        metavar='infiles', 
                        type=str,
                        nargs='+',
                        default=None, 
                        help="One or more stats.<datestr<.json files.")
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    logging.debug(f'infile={args.infiles} outfile={args.outfile}')

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
        
    df = merge_mapseq_dataframes( args.infiles, 
                            outfile=outfile)
    
    write_mapseq_df(df, outfile)
    logging.info('Done aggregate_reads')