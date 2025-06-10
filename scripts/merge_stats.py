#!/usr/bin/env python
#
# merge stats.<datestring>.json files. 
#
import argparse
import json
import logging
import os
import sys
import traceback

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

def merge(a: dict, b: dict, path=[]):
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                merge(a[key], b[key], path + [str(key)])
            elif a[key] != b[key]:
                raise Exception('Conflict at ' + '.'.join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a

def merge_stats(infiles, outfile):
    '''
    combine multiple stats files into one unified (additive) JSON tree.
    needed because some processes put entries into the same top-level dict, 
    e.g. 'fastq' from both process_fastq_pairs, and filter_reads
     
    '''
    merged = None
    for infile in infiles:
        logging.debug(f'reading {infile}')
        with open(infile) as jfh:
            jo = json.load(jfh)
        if merged is None:
            logging.debug(f'merged is None. setting.')
            merged = jo
        else:
            logging.debug(f'merged exists. merging...')
            merged = merge(merged, jo)
    logging.debug(f'finished merging. writing to {outfile}')
    with open(outfile ,'w', encoding='utf-8') as jofh:
        json.dump(merged, jofh, ensure_ascii=False, indent=4)
    logging.debug('done.')

 
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
        
    merge_stats( args.infiles, 
                 outfile=outfile)
    