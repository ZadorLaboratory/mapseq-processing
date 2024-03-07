#!/usr/bin/env python
#

import argparse
import logging
import os
import sys
import traceback

from configparser import ConfigParser

import pandas as pd
import scipy

from scipy.io import loadmat


gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import *
from mapseq.utils import *
from mapseq.stats import *


def dump_matlab(config, infile, outdir=None, expid=None, label=None):
    '''
    Brains = B1 B2 etc...
    barcodematrix  = 
    refbarcodes
    spikes
    B1 = real
    B1norm = filtered and normalized by spikes
    ''' 
    
    if outdir is None:
        outdir = "."
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logging.debug(f'made outdir={outdir}')

    logging.debug(f'dumping {infile} to {outdir}')
    annots = loadmat(infile)
    
    # pull sequences for index. 
    refbarcodes = annots['refbarcodes']
    refseqlist = []
    refseqlist
    for row in refbarcodes:
        s = ''.join([chr(item) for item in row])
        refseqlist.append(s)
    
    outfile = f'{outdir}/refbarcodes.txt'
    writelist(outfile, refseqlist)
    rbcdf = pd.DataFrame(refseqlist, columns=['sequence'])
    outfile = f'{outdir}/refbarcodes.fasta'
    write_fasta_from_df(rbcdf, './refbarcodes.fasta')
    
    barcodematrix = annots['barcodematrix']
    logging.debug(f'bcmatrix shape={barcodematrix.shape}')    
    bdf = pd.DataFrame(barcodematrix, index=refseqlist)
    outfile = f'{outdir}/BarcodeMatrix.tsv'
    bdf.to_csv(outfile, sep='\t')
    logging.debug(f'wrote {outdir}/{key}.norm.tsv')
    for brain_id in ['1','2']:
        key = f'B{brain_id}'
        bm = annots[f'{key}norm']
        seqlist = []
        for row in annots[f'{key}seq']:
            s = ''.join([chr(item) for item in row])
            seqlist.append(s)
        df = pd.DataFrame(bm, index=seqlist)
        outfile = f'{outdir}/{key}.norm.tsv'
        df.to_csv(outfile, sep='\t')
        logging.debug(f'wrote {outdir}/{key}.norm.tsv')
    



    
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

    parser.add_argument('-e','--expid', 
                    metavar='expid',
                    required=False,
                    default='MAPSeq Experiment',
                    type=str, 
                    help='Explicitly provided experiment id, e.g. M205')

    parser.add_argument('-l','--label', 
                    metavar='label',
                    required=False,
                    default='label', 
                    type=str, 
                    help='input column to use to label matrix columns | region [label] ')   

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. input file base dir if not given.')     

    parser.add_argument('-s','--sampleinfo', 
                        metavar='sampleinfo',
                        required=True,
                        default=None,
                        type=str, 
                        help='XLS sampleinfo file. ')

    parser.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='Single MAPseq BarcodeMatrix.mat file.')
       

    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    cp.read(args.config)
    cdict = format_config(cp)    
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infile={args.infile}')
       
    outdir = None
    if args.outdir is not None:
        outdir = os.path.abspath(args.outdir)
        logging.debug(f'making missing outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)
    else:
        afile = args.infile
        filepath = os.path.abspath(afile)    
        dirname = os.path.dirname(filepath)
        outdir = dirname
    cfilename = f'{outdir}/dump_matlab.config.txt'
    #write_config(cp, cfilename, timestamp=True)        
    
    #sampdf = load_sample_info(cp, args.sampleinfo)
    #logging.debug(f'\n{sampdf}')
        
    # create and handle 'real' 'spikein' and 'normalized' barcode matrices...
    dump_matlab(cp, args.infile, outdir, expid=args.expid, label=args.label )
    