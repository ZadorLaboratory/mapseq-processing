#!/usr/bin/env python
#
#   processes per-barcode fasta files. 
#
#    awk "NR%2==0" BC${BCidx[$i]} | grep -v N | cut -b 1-44 | sort | uniq -c | sort -nr > Mseq204_YC_inj2processedBC${BCidx[$i]}.txt
#    #split output files into two files per index, one that is containing the read counts of each unique sequnce, the other the unique sequences themselves.
#    awk '{print $1}' Mseq204_YC_inj2processedBC${BCidx[$i]}.txt > Mseq204_YC_inj2BC${BCidx[$i]}_counts.txt
#    awk '{print $2}' Mseq204_YC_inj2processedBC${BCidx[$i]}.txt > Mseq204_YC_inj2_BC${BCidx[$i]}seq.txt
#
#   BC1.fasta  ->   BC1_processed.tsv
#
#    import resource, sys
#    resource.setrlimit(resource.RLIMIT_STACK, (2**29,-1))
#    sys.setrecursionlimit(10**6)
#


import argparse
import logging
import os
import sys
import traceback

from configparser import ConfigParser

import pandas as pd

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import load_sample_info, process_ssifasta, fix_columns_int, guess_site
from mapseq.utils import *
    
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

    parser.add_argument('-a','--aligner', 
                    metavar='aligner',
                    required=False,
                    default='bowtie', 
                    type=str, 
                    help='aligner tool  [bowtie | bowtie2]')

    parser.add_argument('-m','--max_mismatch', 
                        metavar='max_mismatch',
                        required=False,
                        default=None,
                        type=int, 
                        help='Max mismatch for aligner read collapse.')

    parser.add_argument('-r','--recursion', 
                        metavar='recursion',
                        required=False,
                        default=None,
                        type=int, 
                        help='Max recursion. Handle larger input to collapse() Default is ~5000.')

    parser.add_argument('-s','--sampleinfo', 
                        metavar='sampleinfo',
                        required=True,
                        default=None,
                        type=str, 
                        help='XLS sampleinfo file. ')

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='out file for all merged info')  

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
                        help='SSI fasta file(s)')
       
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
    
    logging.debug(f'recursionlimit = {sys.getrecursionlimit()}')
    if args.recursion is not None:
        rlimit = int(args.recursion)
        logging.info(f'set new recursionlimit={rlimit}')
        sys.setrecursionlimit(rlimit)
    
    outdir = None
    if args.outdir is not None:
        outdir = os.path.abspath(args.outdir)
        logging.debug(f'making missing outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)
    else:
        afile = args.infiles[0]
        filepath = os.path.abspath(afile)    
        dirname = os.path.dirname(filepath)
        outdir = dirname
    
    sampdf = load_sample_info(cp, args.sampleinfo)
    logging.debug(f'\n{sampdf}')
    sampdf.to_csv(f'{outdir}/sampleinfo.tsv', sep='\t')
    
    if args.aligner is not None:
        logging.info(f'setting aligner to {args.aligner}')
        cp.set('ssifasta', 'tool', args.aligner )  

    if args.max_mismatch is not None:
        tool = cp.get('ssifasta','tool')
        mm= str(args.max_mismatch)
        logging.info(f'setting max_mismatch to {mm} tool={tool}')
        cp.set(tool, 'max_mismatch', mm )
        #logging.debug(f"after set. max_mismatch={cp.get('bowtie', 'max_mismatch')} ")   

    cfilename = f'{outdir}/process_ssifasta.config.txt'
    
    write_config(cp, cfilename, timestamp=True)        

    outdflist = []    
    for infile in args.infiles:
        # rtprimer will be guessed from filename. "BC<rtprimer>.fasta"
        (rtprimer, site, brain, region) = guess_site(infile, sampdf)
        logging.info(f'guessed rtprimer={rtprimer} site={site} brain={brain} region={region}')
        try:
            outdf = process_ssifasta(cp, infile, outdir=outdir, site=site)
            logging.debug(f'initial outdf\n{outdf}')
            outdf['site'] = site
            outdf['brain'] = brain
            outdf['region'] = region
            logging.info(f'got outdf:\n{outdf}')
            outdflist.append(outdf)
    
        except Exception as ex:
            logging.warning(f'Possible problem with {infile}. ')
            logging.warning(traceback.format_exc(None))

    if len(outdflist) > 0:                     
        logging.debug(f'merging dfs in outdflist len={len(outdflist)}')
        outdf = merge_dfs(outdflist)
        if args.outfile is None:
            outfile = os.path.join(outdir, 'experiment.all.tsv')
        else:
            outfile = args.outfile
        outdf.to_csv(outfile, sep='\t')          
    