#!/usr/bin/env python
#
#
# h5py example
#
#  infile = 'v73data.mat'
#  f = h5py.File(infile, "r") 
#  list( f.keys())
#  ['#refs#', '#subsystem#', 'chshift_20x']
#  f['chshift_20x'].shape
# (4,1)
#  r1 = f['chshift_20x'][0][0]
#  h1 = f[r1]
#  d1 = h1[:] 
#  type(d1)
#  #numpy.ndarray
#  d1.shape
# (1,6) 
#
#  converting mapseq barcodematrix files via h5py
#
#  f = h5py.File(infile, "r")
#  bcmatrix = np.array(f['barcodematrix']) 
#  df = pd.DataFrame(bcmatrix.T) 
#  df.to_csv('bcmatrix.tsv', sep='\t')
#
#  refbarcodes = np.array(f['refbarcodes'])  
#  rdf = pd.DataFrame( refbarcodes.T )   # ascii integers, 1 per column
#  cdf = pd.DataFrame()
#  for col in rdf.columns:
#     cdf[col] = rdf[col].apply(chr)
#
#  def apply_join(row):
#      return ''.join(row)
#
#  joined = cdf.apply(apply_join, axis=1)
#  joined = joined.astype('string')


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


def map_columns(df, sampdf):
    '''
    for all entries in sampldf.matrixcolumn make mapping to sampldf.rtprimer BC<rtprimer> string. 
    apply map to dataframe columns
    return dataframe with column names in proper order. 
    '''
    mdf = sampdf[['rtprimer','matrixcolumn']].copy()
    for col in mdf.columns:
        mdf[col] = pd.to_numeric(mdf[col], errors='coerce')
    mdf.dropna(inplace = True)
    mdf.matrixcolumn = mdf.matrixcolumn.astype(int)
    mdf = mdf.astype('string')
    col_names = [ str(x) for x in list(df.columns)]
    
    nmap = {}
    for row in mdf.iterrows():
        nmap[ str( row[1]['matrixcolumn'] ) ] =  f"BC{row[1]['rtprimer']}"
    
    new_colnames = []
    for cn in col_names:
        try:
            new_colnames.append( str(nmap[str(cn)] ) )
        except KeyError:
            print(f'error {cn} type={type(cn)}')
            new_colnames.append(cn)
    df.columns = new_colnames
    return df


def get_brain_columns(sampdf, brain_id):
    '''
    
    '''
    logging.debug(f'sampdf len={len(sampdf)} ')     
    bcols = sampdf[sampdf.brain == brain_id]['rtprimer']
    logging.debug(f'{len(bcols)} cols for brain {brain_id} ')
    bcols = list(bcols)
    bcols = [ f'BC{x}' for x in bcols ]
    bcols.sort()
    logging.debug(f'brain columns: {bcols}')
    return bcols


def dump_matlab(config, infile, sampdf=None, outdir=None, expid=None, label=None):
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
    logging.debug(f'refseqlist dataframe={rbcdf}')
    outfile = f'{outdir}/refbarcodes.fasta'
    write_fasta_from_df(rbcdf, './refbarcodes.fasta')
    
    barcodematrix = annots['barcodematrix']
    logging.debug(f'bcmatrix shape={barcodematrix.shape}')    
    bdf = pd.DataFrame(barcodematrix, index=refseqlist)
    logging.debug(f'bcmatrix columns={list(bdf.columns)} converting to string.')
    bdf.columns = list( [str( x + 1)  for x in list(bdf.columns)] )
    
    if sampdf is not None:
        bdf = map_columns(bdf, sampdf)
    outfile = f'{outdir}/BarcodeMatrix.tsv'
    bdf.to_csv(outfile, sep='\t')
    logging.debug(f'wrote {outdir}/BarcodeMatrix.tsv')
    for brain_id in ['1','2','3']:
        key = f'B{brain_id}'
        try:
            bm = annots[f'{key}norm']
            seqlist = []
            for row in annots[f'{key}seq']:
                s = ''.join([chr(item) for item in row])
                seqlist.append(s)
            df = pd.DataFrame(bm, index=seqlist)
            df.columns = list( [str( x + 1)  for x in list(df.columns)] )
            if sampdf is not None:
                df = map_columns(df, sampdf)
                bcols = get_brain_columns(sampdf, brain_id)
                logging.debug(f'df cols = {df.columns}')
                df = df[ bcols ]
            
            outfile = f'{outdir}/{key}.norm.tsv'
            df.to_csv(outfile, sep='\t')
            logging.debug(f'wrote {outdir}/{key}.norm.tsv')
        except KeyError:
            logging.debug(f'No key {key}')
   
   
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
                        required=False,
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
    
    if args.sampleinfo is not None:
        sampdf = load_sample_info(args.sampleinfo, cp=cp)
        logging.debug(f'\n{sampdf}')
    else:
        sampdf = None
        
    # create and handle 'real' 'spikein' and 'normalized' barcode matrices...
    dump_matlab(cp, args.infile, sampdf, outdir, expid=args.expid, label=args.label )
    