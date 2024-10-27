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
from scipy import sparse

import numpy as np 
import h5py


gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)


NTMAP=['N','G','A','T','C']

def map_char(x):
    return NTMAP[x]

def log_objectinfo(obj, label):
    '''
    print info about object to the debug log 
    label should be the variable name for the object
    '''
    ot = type(obj)
    size_gb = ((( sys.getsizeof(obj) ) / 1024 ) / 1024 ) / 1024
    logging.debug(f"variable= '{label}' type={ot} size={size_gb:.4f} GB.")
 

def dump_barseq_matlab(infile, outdir=None):
    '''
    filt_neurons ->  
        ['all_bc', 'all_bc_count', 'angle', 'axonalbc_count', 
        'axonalbc_id', 'axonalbc_pos', 'axonalbc_slice', 'clustid', 
        'clustname', 'depth', 'dom_bc', 'dom_bc_count', 'duplicate_bc', 
        'expmat', 'genes', 'has_proj', 'id', 'is_barcoded', 'pos', 'pos40x', 
        'proj_area', 'proj_brain', 'projmat', 'projmat_norm', 'repeat_bc', 
        'slice', 'soma_bc', 'soma_bc_score', 'soma_bc_sig', 'subclass']
    
        proj_area    list of anatomical regions
        proj_brain
    
    
        soma_bc   
        soma_bc_score
        soma_bc_sig
        

    Q: what is the encoding of soma_bc (1=?  2=? 3=? 4=?)  GATC


    ''' 
    # N included to make integers map naturally by indexing

    
    DUMP_ITEMS=['soma_bc']
    CHAR_COLS = ['soma_bc']
    
    if outdir is None:
        outdir = os.path.abspath("./")
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logging.debug(f'made outdir={outdir}')
    logging.debug(f'dumping {infile} to {outdir}')

    hdfo = h5py.File(infile)

    bdf = extract_bc_df(hdfo)
    of = os.path.join(outdir, 'soma_bc_data.tsv')
    logging.info(f'writing DF to {of}')
    bdf.to_csv(of, sep='\t')

    edf = extract_exp_df(hdfo)
    of = os.path.join(outdir, 'exp_data.tsv')
    logging.info(f'writing DF to {of}')
    edf.to_csv(of, sep='\t')
    
    logging.info(f'files written.')
    

def extract_exp_df(hdfo):
    '''
    https://stackoverflow.com/questions/8403144/loading-matlab-sparse-matrix-saved-with-v7-3-hdf5-into-python-and-operating-o 
    
    extract and label expression matrix. 
    /filt_neurons/expmat/data
                        /ir
                        /jc
    SPARSE
    bc x genes x double 
    
    list(ex.attrs)
        ['H5PATH', 'MATLAB_class', 'MATLAB_sparse']
     
     dereferencing objects:
     https://stackoverflow.com/questions/35490204/how-to-dereference-hdf5-references-in-python
https://stackoverflow.com/questions/59659715/how-to-access-cell-variables-stored-in-mat-file-using-h5py-module 
                        
    '''
    # get genes
    fn = hdfo['filt_neurons']
    
    genes = np.array(fn['genes']).transpose()
    gdf = pd.DataFrame(genes, columns=['gene_sym','gene_seq'])
   
    
    # make expression matrix
    ex = hdfo['filt_neurons/expmat']
    # csc_matrix    
    logging.debug(f'loading as csc_matrix.')
    A = sparse.csc_matrix((ex["data"], ex["ir"], ex["jc"]))
    logging.debug(f'done.')
    log_objectinfo(A, 'csc_matrix')
    
    # coo_matrix
    # logging.debug(f'loading as coo_matrix.')
    # data = np.asarray(ex["data"])
    # ir = np.asarray(ex["ir"])
    # jc = np.asarray(ex["jc"])    
    # B = sparse.coo_matrix(data, (ir, jc))    
    # log_objectinfo(B, 'coo_matrix')
    
    logging.debug('creating sparse dataframe..')
    edf = pd.DataFrame.sparse.from_spmatrix(A)
    logging.debug('done')    
    log_objectinfo(edf, 'csc dataframe')    

    return edf

def extract_bc_df(hdfo):
    '''
    Extract all barcode-oriented data. Columns:
    soma_bc       barcode sequence
    pos_x            
    pos_y
    ...
    
    remember to subtract 1 from any index values (matlab indexes from 1 
    while numpy/python indexes from 0
    
    '''
    fn = hdfo['filt_neurons']
    soma_df = pd.DataFrame()
    
    # arrays of chars/one-character strings
    # soma_bc   soma barcodes
    soma_bc = np.array(fn['soma_bc']).astype(int, copy=False).transpose()
    vfunc = np.vectorize(map_char)
    soma_bc_char = vfunc(soma_bc)
    slist = []
    for row in soma_bc_char:
        slist.append(''.join(row))
    soma_bc_ser = pd.Series(slist)
    soma_df['soma_bc'] = soma_bc_ser
    
    ############### two-columns ###############
    # position.  
    pos = np.array(fn['pos']).transpose()
    pdf = pd.DataFrame(pos, columns=['pos_x','pos_y'])
    soma_df['pos_x'] = pdf['pos_x']
    soma_df['pos_y'] = pdf['pos_y']

    # position 40x.  
    pos40x = np.array(fn['pos40x']).transpose()
    pdf = pd.DataFrame(pos40x, columns=['pos40x_x','pos40x_y'])
    soma_df['pos40x_x'] = pdf['pos40x_x']
    soma_df['pos40x_y'] = pdf['pos40x_y']

    dep = np.array(fn['depth']).transpose()
    pdf = pd.DataFrame(dep, columns=['depth_x','depth_y'])
    soma_df['depth_x'] = pdf['depth_x']
    soma_df['depth_y'] = pdf['depth_y']
    
    ############### single-column ###############
    # id
    id_arr = np.array(fn['id'])
    id_ser = pd.Series(id_arr.ravel())
    soma_df['id'] = id_ser.astype(int)

    # angle 
    angle_arr = np.array(fn['angle'])
    angle_ser = pd.Series(angle_arr.ravel())
    soma_df['angle'] = angle_ser.astype(int)
    
    # slice 
    slice_arr = np.array(fn['slice'])
    slice_ser = pd.Series(slice_arr.ravel())
    soma_df['slice'] = slice_ser.astype(int)    
    
    
    logging.info(f'got soma df len={len(soma_df)}')
    return soma_df

def oldfunc():    
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
    
    parser.add_argument('-e','--expid', 
                    metavar='expid',
                    required=False,
                    default='MAPSeq Experiment',
                    type=str, 
                    help='Explicitly provided experiment id, e.g. M205')

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. working dir if not given.')     

    parser.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='Single BARseq filtered neurons.mat file.')
       

    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    logging.debug(f'infile={args.infile}')
       
    outdir = os.path.expanduser('./')
    if args.outdir is not None:
        outdir = os.path.abspath(args.outdir)
        
    logging.debug(f'making missing outdir: {outdir} ')
    os.makedirs(outdir, exist_ok=True)

   
    # create and handle 'real' 'spikein' and 'normalized' barcode matrices...
    dump_barseq_matlab( args.infile, outdir )
    