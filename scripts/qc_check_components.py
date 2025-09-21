#!/usr/bin/env python
import argparse
import itertools
import logging
import os
import sys
from configparser import ConfigParser
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.collapse import *
from mapseq.core import *
from mapseq.barcode import *
from mapseq.utils import *
from mapseq.stats import *


def assess_components(uniques, component_lists, min_seq_count = 2, top_x = 10):
    '''
    Take top 10 components and calculate fraction that exceed max_hamming. 
    
    '''
    logging.debug(f'assessing uniques/components')
    sh = get_default_stats()
    
    udf = uniques
    comps = component_lists
    logging.debug(f' {len(comps)} components ')
    logging.debug(f' {len(udf)} unique sequences ')
    
    sizemap = []
    for i, clist in enumerate(comps):
        csize = len(clist)
        sizemap.append( [ i, csize ])
    size_df = pd.DataFrame(sizemap, columns=['comp_idx','seq_count'])
    size_df.sort_values('seq_count', ascending=False, inplace=True)
    size_df.reset_index(inplace=True, drop=True)   
    logging.debug(f'stats: \n{size_df.seq_count.describe()}')

    size_df = size_df[size_df['seq_count'] >= min_seq_count ]

    # Assess top X
    for i in range(0, top_x):
        max_comp_idx = size_df.iloc[i].comp_idx
        max_comp = comps[max_comp_idx]
        
        slist = []
        for idx in max_comp:
            s = udf.iloc[idx].vbc_read
            slist.append(s)
        logging.debug(f'made list of len={len(slist)} component sequence list. checking...')
        ( hmax, n_pairs, n_exceed, max_ok)=  max_hamming(slist)
        sh.add_value(f'/collapse/top_{i}','n_seqs', len(slist) )
        sh.add_value(f'/collapse/top_{i}','max_hamming', hmax )
        sh.add_value(f'/collapse/top_{i}','n_pairs', n_pairs )
        sh.add_value(f'/collapse/top_{i}','n_pairs_exceed', n_exceed )
        good_prop = 1.0 - (n_exceed / n_pairs)
        sh.add_value(f'/collapse/top_{i}','pct_good', good_prop )        
    
    return size_df


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

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Component sets.') 

    parser.add_argument('-u', '--uniques',
                        metavar='uniques',
                        required=True,
                        type=str,
                        help='Dataframe TSV with vbc_read column and same index as components.txt')
    
    parser.add_argument('-C','--components',
                        metavar='components',
                        required=True,
                        type=str,
                        help='Text file with component sets. ')
        
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
    logging.debug(f'uniques={args.uniques} components={args.components}')
          
    # set outdir / outfile
    outdir = os.path.abspath('./')
    if args.outfile is not None:
        logging.debug(f'outdir not specified. outfile specified.')
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
    logging.info(f'handling {args.uniques} {args.components}  to outdir {args.outfile}')    


    if args.datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    else:
        datestr = args.datestr
         
    sh = StatsHandler(outdir=outdir, datestr=datestr)
 
    udf = load_mapseq_df( args.uniques, use_dask=False)
    logging.debug(f'loaded. len={len(udf)} dtypes =\n{udf.dtypes}') 
    
    # list of indices.
    components = parse_components(args.components)
           
    # Make final VBC/UMI based table (each row is a neuron)
    logging.debug(f'args={args}')
    results = assess_components(udf,
                                components)
    
    outfile = args.outfile
    logging.info('Done assessing components. ')
    
    