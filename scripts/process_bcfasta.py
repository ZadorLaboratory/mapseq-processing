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


import argparse
import logging
import os
import sys

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from configparser import ConfigParser

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import process_bcfasta
    
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

    parser.add_argument('-a','--aligntool', 
                    metavar='aligntool',
                    required=False,
                    default=None, 
                    type=str, 
                    help='aligner tool  [bowtie | bowtie2]')
    
    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='out file. stdout if not given.')  
    

    parser.add_argument('infile' ,
                        type=str,
                        help='barcode fasta file')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    cp.read(args.config)
    cdict = {section: dict(cp[section]) for section in cp.sections()}
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infiles={args.infile}')
    
    if args.aligntool is not None:
        logging.info(f'setting aligntool to {args.aligntool}')
        cp.set('bcfasta', 'tool', args.aligntool )  
    
    outdf = process_bcfasta(cp, args.infile)
    
    if args.outfile is None:
        print(outdf)
        #outfile = sys.stdout
    #else:
    #    outfile = args.outfile        
    #outdf.to_csv(outfile, sep='\t', index=False)    
    
    
    