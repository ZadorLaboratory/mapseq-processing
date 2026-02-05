import copy
import glob
import gzip
import itertools
import json
import logging
import math
import os
import random
import subprocess
import sys
import traceback

import datetime as dt

from configparser import ConfigParser
from collections import defaultdict

import pandas as pd
import numpy as np
from natsort import natsorted, index_natsorted, order_by_index

import scipy

import dask
import dask.dataframe as dd
from dask.dataframe import from_pandas
from dask.distributed import Client

from mapseq.utils import *
from mapseq.bowtie import *
from mapseq.barcode import *
from mapseq.stats import *
from mapseq.plotting import *
from mapseq.collapse import *

#
# STANDARD FORMATS AND DATATYPES
#

FFORMATS = ['reads','aggregated','filtered','readtable','collapsed','vbctable']

STR_COLS = {
    'reads'      : ['sequence'],
    'aggregated' : ['sequence'],
    'filtered'   : ['vbc_read', 'spikeseq', 'libtag', 'umi',  'ssi'],
    'readtable'  : ['vbc_read', 'libtag', 'umi' ],
    'collapsed'  : ['vbc_read', 'libtag', 'umi' ],
    'vbctable'   : ['vbc_read'],      
    }

INT_COLS = {
    'reads'      : [],
    'aggregated' : ['read_count'],
    'filtered'   : ['read_count'],
    'readtable'  : ['read_count'],
    'vbctable'   : ['read_count','umi_count'],
    }

CAT_COLS = {
    'reads'      : ['source'],
    'aggregated' : ['source'],
    'filtered'   : ['source'],
    'readtable'  : ['label','site','type','brain','region','source','ourtube','rtprimer'],
    'vbctable'   : ['label','site','type','brain','region','ourtube'],        
    }

FMT_DTYPES = {      'read_count'    : 'int64',
                    'umi_count'     : 'int64',
                    'sequence'      : 'string[pyarrow]',
                    'vbc_read'      : 'string[pyarrow]',
                    'vbc_read_col'  : 'string[pyarrow]',
                    'vbc_read_short': 'string[pyarrow]',        
                    'spikeseq'      : 'string[pyarrow]',    
                    'libtag'        : 'string[pyarrow]',    
                    'umi'           : 'string[pyarrow]',    
                    'ssi'           : 'string[pyarrow]',
                    'cell_id'       : 'string[pyarrow]',
                    'source'        : 'category',
                    'label'         : 'category',
                    'rtprimer'      : 'category',
                    'type'          : 'category',
                    'site'          : 'category',
                    'brain'         : 'category',
                    'region'        : 'category',
                    'ourtube'       : 'category',

    }

#                 Sample Information metadata site labels
# LABEL                        SOURCE                MEANING
# EXPERIMENTAL
# target                    user-provided        treated sample
# target-negative           user-provided        treated but expected to be low
#       (Both handled exactly the same and included in matrices) 
# injection                 user-provided         treated sample
#
# CONTROLS
# target-negative-control   user-provided       untreated sample
# target-wt-control         core added          untreated biological sample
# target-water-control      core added          empty sample
# injection-water-control   core added          empty sample

#CONTROL_SITES=['target-negative-control', 
#               'target-wt-control',
#               'target-water-control',
#               'injection-water-control',
#               'injection-wt-control']

CONTROL_SITES=['target-wt-control',
               'target-water-control',
               'injection-water-control',
               'injection-wt-control',
               'target-lone']

#
#             MAPSEQ UTILITY 
#

def get_default_config():
    dc = os.path.expanduser('~/git/mapseq-processing/etc/mapseq.conf')
    cp = ConfigParser()
    cp.read(dc)
    return cp

    
def load_mapseq_df( infile, fformat=None, use_dask=False, use_arrow=True):
    '''
    Abstracted loading code for all MAPseq pipeline dataframe formats. 
    Uses dtypes above FMT_DTYPES
    
    Fixes NaNs in categoricals. 

    '''    
    logging.info(f"loading {infile} format='{fformat}' use_dask={use_dask} use_arrow={use_arrow}")
    if fformat is None:
        logging.info(f'fformat is None. Will not apply datatype checks/correction.')
    elif fformat not in FFORMATS :
        logging.warning(f'fformat {fformat} not in {FFORMATS} . Datatypes may be incorrect.')
    
    ftype = None
    df = None
    
    if infile.endswith('.tsv'):
         ftype = 'tsv'
    elif infile.endswith('.parquet'):
       ftype = 'parquet'
    else:
        logging.error('input file must have extension .tsv or .parquet')
        sys.exit(1)
    logging.debug(f'input filetype={ftype}')
    
    start = dt.datetime.now()
    
    if ftype == 'tsv':
        if use_dask:
            logging.debug(f'loading via Dask ftype=tsv use_dask=True with dtype dict...')
            df = dd.read_csv(infile, sep='\t', dtype=FMT_DTYPES )  
        else:
            logging.debug(f'loading via Pandas ftype=tsv use_dask=False with dtype dict...')
            df = pd.read_csv(infile, sep='\t', index_col=0, dtype=FMT_DTYPES )
     
        logging.debug(f'after: dtypes=\n{df.dtypes}')        
    
    elif ftype == 'parquet':
        if use_dask:
            if use_arrow:
                logging.debug(f'loading via Dask ftype=parquet use_dask=True filesystem=arrow')
                df = dd.read_parquet(infile, filesystem='arrow')
            else:
                logging.debug(f'loading via Dask ftype=parquet use_dask=True filesystem=fsspec (default)')
                df = dd.read_parquet(infile)  
        else:
            logging.debug(f'loading via Pandas ftype=parquet use_dask=False use_arrow={use_arrow}')
            df = pd.read_parquet(infile)
    logging.debug(f'fixing dataframe types')
    df = fix_mapseq_df_types(df, 
                             fformat=fformat, 
                             use_arrow=use_arrow) 
    
    end = dt.datetime.now()
    delta_seconds = (dt.datetime.now() - start).seconds    
    log_transferinfo(infile, delta_seconds)
    log_objectinfo(df, 'loaded_df')
    return df       


def write_mapseq_df(df, outfile, outformats=['tsv','parquet'] ):
    '''
    Wrapper to time writeout performance. 
    Assumes outfile is TSV string, but will adjust as needed. 
    
    '''
    log_objectinfo(df, 'write_df')
    outfile = os.path.abspath(outfile)
    dirname, base, ext = split_path(outfile)
    logging.debug(f'dirname={dirname} base={base} ext={ext}')
    
    if 'tsv' in outformats:
        of = os.path.join(dirname, f'{base}.tsv')
        logging.info(f'Saving len={len(df)} as TSV to {of}...')
        start = dt.datetime.now()
        df.to_csv(of, sep='\t')
        
        end = dt.datetime.now()
        delta_seconds = (dt.datetime.now() - start).seconds
        log_transferinfo(of, delta_seconds)
        logging.info(f'Done writing {of}')
    
    if 'parquet' in outformats:
        of = os.path.join(dirname, f'{base}.parquet')
        logging.info(f'Saving len={len(df)} as Parquet to {of}...')
        start = dt.datetime.now()
        df.to_parquet(of)        
        end = dt.datetime.now()
        delta_seconds = (dt.datetime.now() - start).seconds
        log_transferinfo(of, delta_seconds)
        logging.info(f'Done writing {of}')     
    logging.info(f'Done writing DF to desired format(s).')

    
def fix_mapseq_df_types(df, 
                        fformat='reads', 
                        use_arrow=True):
    '''
    confirm that columns are proper types. 
    use string[pyarrow] for strings, else string
    
    Fix NaNs. 
    
    '''
    logging.info(f'old dataframe dtypes=\n{df.dtypes}')
    
    if fformat in FFORMATS:
        if use_arrow:
            tt = 'string[pyarrow]'
        else:
            tt = 'string[python]'
            
        for scol in STR_COLS[fformat]:
            try:
                dt = df[scol].dtype 
                if dt != tt:
                    logging.debug(f'converting col={scol} from {dt} to {tt} ...')
                    df[scol] = df[scol].astype(tt) 
            except KeyError:
                logging.warning(f'column {scol} not found. Vital for {fformat}')
            
        for icol in INT_COLS[fformat]:
            dt = df[icol].dtype
            if dt != 'int64':
                logging.debug(f"converting col={icol} from {dt} to 'int64' ...")
                df[icol] = df[icol].astype('int64')                
            
        for ccol in CAT_COLS[fformat]:
            try:
                dt = df[ccol].dtype
                if dt != 'category':
                    logging.debug(f"converting col={ccol} from {dt} to 'category' ...")
                    df[ccol] = df[ccol].astype('category')
                    try:
                        df[ccol] = df[ccol].cat.add_categories([''])
                    except ValueError:
                        logging.debug(f"category ''  already exists." )                
                else:
                    try:
                        df[ccol] = df[ccol].cat.add_categories([''])
                    except ValueError:
                        logging.debug(f"category ''  already exists." )
                #df[ccol].fillna('', inplace=True)                        
                df.fillna( {ccol: ''}, inplace=True)
            except KeyError:
                logging.warning(f'column {ccol} not found. Vital for {fformat}?')        
            
    else:
        logging.warning('unrecognized mapseq format. return original')
    logging.info(f'new dataframe dtypes=\n{df.dtypes}')
    has_nan = df.isna().any().any()
    if has_nan:
        logging.warning(f'Dataframe has NaN values, columns={df.isna().any()}')
    else:
        logging.info(f'No NaN values in dataframe.')
    
    return df


def load_mapseq_matrix_df( infile, use_dask = False ):
    '''
    convenience method to align with load_mapseq_df()
    just ensures validity, data type. 
    
    '''
    df = None
    if infile.endswith('.tsv'):
        if use_dask:
            df = dd.read_csv(infile, sep='\t')  
        else:
            df = pd.read_csv(infile, sep='\t', index_col=0)

    elif infile.endswith('.parquet'):
        if use_dask:
            df = dd.read_parquet(infile)  
        else:
            df = pd.read_parquet(infile)
    else:
        logging.error('input file must have relevant extension .tsv or .parquet')
        sys.exit(1)
    
    logging.debug(f'naive: dtypes=\n{df.dtypes}')    
    return df    



#
#          SAMPLEINFO/METADATA UTILITY
#
def load_sample_info(file_name, 
                     sheet_name='Sample information', 
                     cp=None):
    #
    # Parses Excel spreadsheet to get orderly sample metadata, saves as ./sampleinfo.tsv.     
    # OR Reads in sampleinfo.tsv
    # Assumes various properties of spreadsheet that need to stay static. 
    #
    #   ['Tube # by user', 'Our Tube #', 'Sample names provided by user',
    #   'Site information', 'RT primers for MAPseq', 'Brain ', 'Column#']
    #
    # If brain is not given, or is empty, all are set to 'brain1'. 
    # If region is not given, or is empty, all are set to <rtprimer>
    # 
    
    if cp is None:
        cp = get_default_config()
    
    # Mappings for excel columns. 
    sheet_to_sample = {
            'Tube # by user'                  : 'usertube', 
            'Our Tube #'                      : 'ourtube', 
            'Sample names provided by user'   : 'samplename', 
            'Site information'                : 'siteinfo',
            'Spike-in Ratio'                  : 'si_ratio',
            'Read Count Minimum'              : 'min_reads',
            'RT primers for MAPseq'           : 'rtprimer',
            'Brain'                           : 'brain',
            'Region'                          : 'region',
            'Matrix Column'                   : 'matrixcolumn',
        }
    
    sample_columns = ['usertube', 'ourtube', 'samplename', 'siteinfo', 'si_ratio', 'rtprimer', 'brain', 'region', 'matrixcolumn', 'min_reads'] 

    if file_name.endswith('.xlsx'):
        edf = pd.read_excel(file_name, sheet_name=sheet_name, header=1, dtype=str)        
        sdf = pd.DataFrame()
        
        for ecol in edf.columns:
            ecol_stp = ecol.strip()    
            try:
                # map using stripped column name, retrieve using actual excel column name
                # which may have trailing spaces...
                scol = sheet_to_sample[ecol_stp]
                logging.debug(f'found mapping {ecol} -> {scol}')
                cser = edf[ecol]
                sdf[scol] = cser
            
            except KeyError:
                logging.debug(f'no mapping for {ecol} continuing...')
            
            except Exception as ex:
                logging.error(f'error while handling {ecol} ')
                logging.info(traceback.format_exc(None))
        
        # Only keep rows with rtprimer info. 
        sdf = sdf[sdf['rtprimer'].isna() == False]
        sdf.fillna('', inplace=True)

        # Ensure required columns are present, fill in with reasonable values.         
        for scol in sample_columns:
            try:
                ser = sdf[scol]
            except KeyError as ke:
                logging.info(f'no column {scol}, required. Creating...')
                if scol == 'samplename':
                    sdf[scol] = sdf['ourtube']
                elif scol == 'region':
                    sdf[scol] = sdf['rtprimer']
                elif scol == 'si_ratio':
                    sdf[scol] = 1.0
                elif scol == 'min_reads':
                    logging.info('No per-sample min_reads column. setting from config...')
                    sdf = set_min_reads(sdf, cp=cp)
        logging.info(f'loaded DF from Excel {file_name}')

        # Fix missing values with defaults. 
        mmap = sdf['min_reads'] == ''
        sdf.loc[mmap, 'min_reads'] = 1
        
        rmap = sdf['si_ratio'] == ''
        sdf.loc[rmap, 'si_ratio'] = 1.0
        
    elif file_name.endswith('.tsv'):
        sdf = pd.read_csv(file_name, sep='\t', index_col=0, keep_default_na=False, dtype =str, comment="#")
        sdf = sdf.astype('str', copy=False)    
    else:
        logging.error(f'file {file_name} neither .xlsx or .tsv')
        sdf = None

    # fix datatypes for important columns
    sdf['min_reads'] = sdf['min_reads'].astype(int)
    sdf['si_ratio'] = sdf['si_ratio'].astype(float)
    
    logging.debug(f'created reduced sample info df:\n{sdf}')
    return sdf


def set_min_reads(df, cp):
    '''
    apply target/injection min_reads from config if not provided
    explicitly in sampleinfo. 
    '''
    inj_min_reads = int(cp.get('vbctable','inj_min_reads'))
    target_min_reads = int(cp.get('vbctable','target_min_reads'))
    logging.info(f'setting min_reads: target_min_reads={target_min_reads} inj_min_reads={inj_min_reads}')
    df['min_reads'] = 1
    tmask = df['siteinfo'].str.startswith('target')
    imask = df['siteinfo'].str.startswith('injection')
    df.loc[tmask,'min_reads'] = target_min_reads
    df.loc[imask, 'min_reads'] = inj_min_reads
    return df


def get_rtlist(sampledf):
    '''
    handle numeric rtprimer convert to BC<int>
    otherwise, return literal list of labels. 
    
    '''
    rtlist = list(sampledf['rtprimer'].dropna())
    nrtlist = []
    for x in rtlist:
        try:
            y = int(float(x))
            nrtlist.append(y)
        except ValueError:
            logging.debug(f'ignoring bad int() input.')
        
    if len(nrtlist) > 1: 
        nrtlist = [f'BC{x}' for x in nrtlist]
    else:
        nrtlist = list(sampledf['rtprimer'].dropna())
    logging.debug(f'got rtlist = {nrtlist}')
    return nrtlist





############################################################
#
#                       Pipeline process phases
#
############################################################
def process_fastq_pairs(infilelist, 
                        outdir,                         
                        force=True,
                        min_reads=None,
                        cp = None ):
    '''
    Top-level command entry point to simplify command line script. 
    
    Originally used to combine operations, but these must be separate with large datasets. 
    
    Pull out relevant leading sections of read1 and read2, join them  
    df = process_fastq_pairs_pd_join( infilelist, 
                                      outdir, 
                                      force=args.force, 
                                      cp=cp) 
    '''
    if cp is None:
        cp = get_default_config()
    
    logging.info('Handling FASTQ parsing...')
    df = process_fastq_pairs_pd_chunked( infilelist,
                                      outdir, 
                                      force=force, 
                                      cp=cp)    
    logging.info(f'done with FASTQ parsing.')
    return df


def process_fastq_pairs_pd_chunked( infilelist, 
                                    outdir,                         
                                    force=False, 
                                    cp = None):
    '''
    only parse out read lines to pandas, then join with pandas. 
    also AGGREGATE by identical read to minimize data. read_count column is counts of full sequence. 
    handle files with chunking so that strings can be squished to pyarrow. 

    Output:
        sequence    read_count
    
    df = pd.read_csv(fh, header=None, skiprows = lambda x: x % 4 != 1, dtype="string[pyarrow]")
    df_iter = pd.read_csv(fh, header=None, skiprows = lambda x: x % 4 != 1, dtype="string[pyarrow]", chunksize=1000000)
    
    https://pythonspeed.com/articles/chunking-pandas/
    
    DASK??
    https://saturncloud.io/blog/how-to-efficiently-read-large-csv-files-in-python-pandas/

    '''
    r1s = int(cp.get('fastq','r1start'))
    r1e = int(cp.get('fastq','r1end'))
    r2s = int(cp.get('fastq','r2start'))
    r2e = int(cp.get('fastq','r2end'))
    
    chunksize = int(cp.get('fastq','chunksize'))
    source_regex = cp.get('fastq','source_regex')
    logging.info(f'chunksize={chunksize} lines.')
    logging.debug(f'read1[{r1s}:{r1e}] + read2[{r2s}:{r2e}]')
    df = None
    sh = get_default_stats()
    pairnum = 1
    chunknum = 1
    old_len = 0

    for (read1file, read2file) in infilelist:
        source_label = parse_sourcefile(read1file, source_regex)
        logging.info(f'handling {read1file}, {read2file} source_label={source_label}')
        start = dt.datetime.now()
        fh1 = get_fh(read1file)
        dfi1 = pd.read_csv(fh1, 
                           header=None, 
                           skiprows = lambda x: x % 4 != 1, 
                           dtype="string[pyarrow]", 
                           chunksize=chunksize) 
        fh2 = get_fh(read2file)
        dfi2 = pd.read_csv(fh2, 
                           header=None, 
                           skiprows = lambda x: x % 4 != 1, 
                           dtype="string[pyarrow]", 
                           chunksize=chunksize)
        
        for chunk1 in dfi1:
            chunk2 = dfi2.get_chunk()
            logging.debug(f'chunk1={chunk1} chunk2={chunk2}')
            if df is None:
                logging.debug(f'making new read sequence DF...')
                df = pd.DataFrame(columns=['sequence'])
                logging.info(f'got chunk len={len(chunk1)} slicing...')
                df['sequence'] = chunk1[0].str.slice(r1s, r1e) + chunk2[0].str.slice(r2s, r2e) 
                logging.info(f'chunk df type={df.dtypes}')           
                df['source'] = source_label
            else:
                logging.debug(f'making additional read sequence DF...')
                ndf = pd.DataFrame(columns=['sequence'], dtype="string[pyarrow]")
                logging.info(f'got chunk len={len(chunk1)} slicing...')
                ndf['sequence'] = chunk1[0].str.slice(r1s, r1e) + chunk2[0].str.slice(r2s, r2e) 
                logging.info(f'chunk df type={ndf.dtypes}')                
                ndf['source'] = source_label                
                df = pd.concat([df, ndf], copy=False, ignore_index=True)
            logging.debug(f'handled chunk number {chunknum}')
            chunknum += 1

        logging.debug(f'handled pair number {pairnum}')
        
        # Measure file read speed.
        end = dt.datetime.now()
        delta_seconds = (dt.datetime.now() - start).seconds
        log_transferinfo( [read1file, read2file] , delta_seconds)
        pair_len = len(df) - old_len
        old_len = len(df)
        sh.add_value('/fastq',f'pair{pairnum}_len', pair_len )
        pairnum += 1


    logging.debug(f'dtypes =\n{df.dtypes}')
    logging.info('Finished processing all input.')
    sh.add_value('/fastq','reads_handled', len(df) )
    sh.add_value('/fastq','pairs_handled', pairnum )
    return df          


def process_fastq_grouped(   infilelist, 
                             outdir,
                             sampdf=None,                         
                             force=False, 
                             cp = None,
                             ):
    '''
    parse infiles by pairs. 
    get extra field(s)
    optionally filter by extra fields values. 
    Useful for Novaseq raw data, with 1 sample per FASTQ pair. 
        
    '''
    r1s = int(cp.get('fastq','r1start'))
    r1e = int(cp.get('fastq','r1end'))
    r2s = int(cp.get('fastq','r2start'))
    r2e = int(cp.get('fastq','r2end'))
    
    chunksize = int(cp.get('fastq','chunksize', fallback=50000))
    source_regex = cp.get('fastq','source_regex', fallback='(.+?_S.+?_L.+?)_')
    ourtube_regex = cp.get('fastq','ourtube_regex', fallback='.+?-(.+?)_')

    logging.info(f' source_regex = {source_regex} ourtube_regex = {ourtube_regex}')

    write_each = cp.getboolean('fastq','write_each', fallback=False)
    filter_by_non_dominant = cp.getboolean('fastq','filter_by_non_dominant',fallback=False )
    filter_column = cp.get('fastq','filter_column', fallback='ssi')
    drop_filter_column = cp.getboolean('fastq','drop_filter_column', fallback=False)
    filter_by_sampleinfo_ssi = cp.getboolean('fastq','filter_by_sampleinfo_ssi', fallback=False)
    if sampdf is None:
        filter_by_sampleinfo_ssi = False
        logging.warning('No sampleinfo provided. Filename-based SSI filtering not possible.')
    logging.info(f' write_each = {write_each} filter_by_non_dominant={filter_by_non_dominant} filter_column={filter_column} drop_filter_column={drop_filter_column} filter_by_sampleinfo_ssi={filter_by_sampleinfo_ssi} chunksize={chunksize} lines.')
    logging.debug(f'read1[{r1s}:{r1e}] + read2[{r2s}:{r2e}]')
    
    outdf = None
    sh = get_default_stats()
    pairnum = 1
    chunknum = 1

    for (read1file, read2file) in infilelist:
        df = None
        source_label = parse_sourcefile(read1file, source_regex, 1 )
        if filter_by_sampleinfo_ssi:
            source_tube = parse_sourcefile(read1file, ourtube_regex, 1 )
            logging.info(f'source_tube={source_tube}')
        logging.info(f'handling {read1file}, {read2file} source_label={source_label} ')
        start = dt.datetime.now()
        fh1 = get_fh(read1file)
        dfi1 = pd.read_csv(fh1, 
                           header=None, 
                           skiprows = lambda x: x % 4 != 1, 
                           dtype="string[pyarrow]", 
                           chunksize=chunksize) 
        fh2 = get_fh(read2file)
        dfi2 = pd.read_csv(fh2, 
                           header=None, 
                           skiprows = lambda x: x % 4 != 1, 
                           dtype="string[pyarrow]", 
                           chunksize=chunksize)
        
        for chunk1 in dfi1:
            chunk2 = dfi2.get_chunk()
            logging.debug(f'chunk1={chunk1} chunk2={chunk2}')
            if df is None:
                logging.debug(f'making new read sequence DF...')
                df = pd.DataFrame(columns=['sequence'])
                logging.info(f'got chunk len={len(chunk1)} slicing...')
                df['sequence'] = chunk1[0].str.slice(r1s, r1e) + chunk2[0].str.slice(r2s, r2e) 
                logging.info(f'chunk df type={df.dtypes}')           
                df['source'] = source_label
            else:
                logging.debug(f'making additional read sequence DF...')
                ndf = pd.DataFrame(columns=['sequence'], dtype="string[pyarrow]")
                logging.info(f'got chunk len={len(chunk1)} slicing...')
                ndf['sequence'] = chunk1[0].str.slice(r1s, r1e) + chunk2[0].str.slice(r2s, r2e) 
                logging.info(f'chunk df type={ndf.dtypes}')                
                ndf['source'] = source_label                
                df = pd.concat([df, ndf], copy=False, ignore_index=True)
            logging.debug(f'handled chunk number {chunknum}')
            chunknum += 1
        # We now have all input from 2 input files.
        
        # split out requested fields. 
        logging.info(f'getting additional field(s)...')
        fieldlist = add_split_fields(df, 'sequence', cp, 'fastq')
        
        # save to file.
        if write_each: 
            logging.info(f'write_each={write_each}  Saving to readfile...')
            of = os.path.join(outdir, f'{source_label}.reads.tsv')
            write_mapseq_df(df, of)
        else:
            logging.debug(f'write_each={write_each} skipping.')
        
        # gather extra field stats. 
        for field in fieldlist:
            try:
                vc = df[field].value_counts()
                n_dom = vc.iloc[0]
                logging.debug(f'dominant value count for {field} = {n_dom} ') 
                pct_dom = n_dom / len(df) * 100
                spct = f'{pct_dom:.2f}'
                sh.add_value('/fastq', f'pair{pairnum}_{source_label}_{field}_percent', spct )
            except Exception as e:
                logging.warning(f'issue dealing with source={source_label} field={field}')
        logging.debug(f'handled pair number {pairnum}')
        
        # Measure file read speed.
        end = dt.datetime.now()
        delta_seconds = (dt.datetime.now() - start).seconds
        log_transferinfo( [read1file, read2file] , delta_seconds)
        sh.add_value('/fastq',f'pair{pairnum}_len', len(df) )
        pairnum += 1

        # Optionally handle per-file non-dominant filtering.
        if filter_by_non_dominant:
            filter_column = cp.get('fastq','filter_column')
            logging.debug(f'filtering by dominant value in source.')
            df = filter_non_dominant( df, 
                                      filter_column=filter_column, 
                                      drop_filter_column=drop_filter_column)
        else:
            logging.debug(f'no filtering by non-dominant column value.')

        if filter_by_sampleinfo_ssi:
            logging.debug(f'filtering by ssi. source_label={source_label} source_tube={source_tube}')
            filter_column = cp.get('fastq','filter_column')
            df = filter_by_sampleinfo(df,
                                      sampdf = sampdf,
                                      ourtube = source_tube,
                                      cp=cp
                                      )
        else:
            logging.debug(f'no filtering by sampleinfo SSI.') 

        logging.debug(f'continuing creation of full outdf...')
        outdf = pd.concat([outdf, df], copy=False, ignore_index=True)

    logging.debug(f'dtypes =\n{outdf.dtypes}')
    logging.info('Finished processing all input.')
    sh.add_value('/fastq','reads_handled', len(outdf) )
    sh.add_value('/fastq','pairs_handled', pairnum )
    outdf.reset_index(inplace=True, drop=True)
    return outdf          


def filter_non_dominant(df,
                        filter_column = 'ssi',
                        drop_filter_column = True,
                        ):
    '''
    filter rows of df that do not have the dominant value of <dom_column>
    '''
    init_len = len(df)
    logging.debug(f'handling df len={init_len}')
    vcser = df[filter_column].value_counts()
    dom_val = vcser.index[0]
    logging.debug(f'dominant value={dom_val}')
    df = df[ df[filter_column]  ==  dom_val ]
    logging.debug(f'filtered len={len(df)}')
    if drop_filter_column:
        df = df.drop(filter_column, axis=1)
    return df 


def filter_by_sampleinfo(df,
                         sampdf,
                         ourtube,
                         cp=None
                         ):
    '''
    fairly hard-coded, assumes 'ssi' column. 
    assign rtprimer from ssi sequence
    choose rtprimer form sampdf given ourtube value. 
    only retain rows with the correct rtprimer.

    '''
    if cp is None:
        cp = get_default_config()
    logging.debug(f'len(df)={len(df)} ourtube={ourtube} type(ourtube)={type(ourtube)} ')
    bcfile = os.path.expanduser( cp.get('barcodes','ssifile') )    
    logging.debug(f'bcfile={bcfile}')

    drop_sampleinfo_columns = cp.getboolean('fastq','drop_sampleinfo_columns')

    # Map label, rtprimer to SSIs    
    logging.debug(f'getting rt labels...')
    labels = get_rtlist(sampdf)
    logging.debug(f'rtlabels={labels}')
    bcdict = get_barcode_dict(bcfile, labels)
    rtdict = get_rtprimer_dict(bcfile, labels)
    rtbcdict = get_rtbc_dict(bcfile, labels)
    logging.debug(f'got {len(bcdict)} barcodes with labels and primer number.')    
   
    logging.info('filling in rtprimer number by SSI sequence...')    
    df['rtprimer'] = df['ssi'].map(rtdict)
    df['rtprimer'] = df['rtprimer'].astype('category')
    
    # Normalize case for both filename and sampleinfo
    ourtube = ourtube.lower()
    sampdf['ourtube'] = sampdf['ourtube'].str.lower()
    
    sr = sampdf.loc[ sampdf['ourtube'] == ourtube]
    logging.debug(f'ourtube={ourtube} row={sr}')
    sv = sr['rtprimer']
    logging.debug(f'rtprimer ser = {sv}')
    if len(sv) > 0:
        rtprimer_val = sv.iloc[0]
        logging.info(f'rtprimer value = {rtprimer_val}')

        df = df[df['rtprimer'] == rtprimer_val]
        df.reset_index(inplace=True, drop=True)
    else:
        logging.warning(f'no ')
    
    if drop_sampleinfo_columns:
        logging.debug('dropping ssi, rtprimer')
        df = df.drop(['ssi','rtprimer'], axis=1)   
    logging.debug(f'filtered len(df)={len(df)} ')
    return df
    

def add_split_fields(df, column, cp, section):
    '''
    given a dataframe, column, and a configparser w/ section
    pull out all items in section with _st and _end, and create
    new columns in <df> from <column>. 
    
    '''
    plist = cp.items(section)
    fdict = defaultdict(dict)
    for (k,v) in plist:
        if k.endswith('_st'):
           fname =  k[:-3]
           fd = fdict[fname]
           fd['start'] = v
           
        elif k.endswith('_end'):
           fname =  k[:-4]
           fd = fdict[fname]
           fd['end'] = v
    
    for field in list(fdict.keys()):
        logging.debug(f'handling field={field}')
        try:
            fd = fdict[field]
            logging.debug(f'info for field={field} = {fdict[field]}')
            fstart = int( fdict[field]['start'])
            fend = int( fdict[field]['end'] )
            logging.debug(f'splitting field={field} start = {fstart} end = {fend}')
            df[field] = df[column].str.slice(fstart, fend).astype('string[pyarrow]')
            logging.debug(f'successfully split field={field}')
        
        except Exception as ex:
            logging.warning(f'problem with field={field} bad start/end? Check config?')
    
    return list(fdict.keys())




#
#        AGGREGATE  READS 
#

def aggregate_reads(df, 
                    column=['sequence','source'],
                    outdir=None, 
                    min_reads=None, 
                    use_dask=None, 
                    dask_temp = None, 
                    cp=None,
                    ):
    '''
    Top level entry point for aggregate. Chooses underlying algorithm/system. 
    '''
    if cp is None:
        cp = get_default_config()
    
    if use_dask is None:
        use_dask = cp.getboolean('fastq', 'use_dask')
    
    if use_dask:
        if dask_temp is None:
            dask_temp = os.path.abspath(os.path.expanduser( cp.get('fastq','dask_temp')))
    
    if min_reads is None:
        min_reads = int(cp.get('fastq','min_reads'))
    else:
        min_reads = int(min_reads)
    
    if outdir is None:
        outdir = os.path.abspath('./')
    logging.info(f'aggregate_reads: use_dask={use_dask} dask_temp={dask_temp} min_reads={min_reads}')
    
    if use_dask:
        df = aggregate_reads_dd(df, 
                                column, 
                                outdir=outdir, 
                                min_reads=min_reads, 
                                dask_temp=dask_temp, 
                                cp=cp )
    else:
        df = aggregate_reads_pd(df, 
                                column, 
                                outdir=outdir, 
                                min_reads=min_reads )
    return df
    
    
    
def aggregate_reads_pd(df, 
                       column=['sequence','source'], 
                       outdir=None, 
                       min_reads=1):
    initlen = len(df)
    logging.debug(f'aggregating read counts DF len={len(df)} column={column}')
    try:
        # some test files may only have 'sequence' column.
        vcs = df.value_counts( column )
    except KeyError:
        vcs = df.value_counts( column[0])
    ndf = pd.DataFrame( vcs )
    ndf.reset_index(inplace=True, drop=False)
    ndf.rename({'count':'read_count'}, inplace=True, axis=1)
    logging.debug(f'DF len={len(ndf)}')
    
    if min_reads > 1:
        logging.info(f'Dropping reads with less than {min_reads} read_count.')
        logging.debug(f'Length before read_count threshold={len(ndf)}')
        ndf = ndf[ndf['read_count'] >= min_reads]
        ndf.reset_index(inplace=True, drop=True)
        logging.info(f'Length after read_count threshold={len(ndf)}')    
    else:
        logging.info(f'min_reads = {min_reads} skipping initial read count thresholding.')  
    logging.info(f'final output DF len={len(ndf)}')
    
    return ndf

def aggregate_reads_dd(seqdf, 
                       column=['sequence','source'], 
                       outdir=None, 
                       min_reads=1, 
                       chunksize=50000000,
                       dask_temp=None, 
                       cp=None):
    '''
    ASSUMES INPUT IS DASK DATAFRAME
    retain other columns and keep first value

    Suggestions for partitions. 
    https://stackoverflow.com/questions/44657631/strategy-for-partitioning-dask-dataframes-efficiently

    Limiting resource usage (CPU, memory)
    https://stackoverflow.com/questions/59866151/limit-dask-cpu-and-memory-usage-single-node
    
    '''
    if cp is None:
        cp = get_default_config()
    
    if dask_temp is None:
        dask_temp = os.path.abspath( os.path.expanduser(cp.get('fastq','dask_temp')))
    logging.debug(f'setting dask temp. dask_temp={dask_temp}')
    
    if not os.path.isdir(dask_temp):
        os.makedirs(dask_temp, exist_ok=True)        
    dask.config.set(temporary_directory=dask_temp)
        
    initlen = len(seqdf)
    logging.debug(f'aggregating read counts for col={column} len={len(seqdf)}')
    #ndf = seqdf[column].value_counts().compute()
    ndf = seqdf.value_counts(column).compute()
    ndf = ndf.reset_index()
    
    ndf.rename({'count':'read_count'}, inplace=True, axis=1)
    logging.info(f'computed counts. new DF len={len(ndf)}')
    
    # get back all other columns from original set. e.g. 'source'
    logging.info(f'merging to recover other columns from original DF')
    ndf = dd.from_pandas(ndf)
    #ndf = pd.merge(ndf, seqdf.drop_duplicates(subset=column,keep='first'),on=column, how='left')  
    result = ndf.merge(seqdf.drop_duplicates(subset=column, keep='first'), on=column, how='left')  
    ndf = result.compute()
    ndf.drop('Unnamed: 0', inplace=True, axis=1 )
    logging.info(f'got merged DF=\n{ndf}')
  
    if min_reads > 1:
        logging.info(f'Dropping reads with less than {min_reads} read_count.')
        logging.debug(f'Length before read_count threshold={len(ndf)}')
        ndf = ndf[ndf['read_count'] >= min_reads]
        ndf.reset_index(inplace=True, drop=True)
        logging.info(f'Length after read_count threshold={len(ndf)}')    
    else:
        logging.info(f'min_reads = {min_reads} skipping initial read count thresholding.')  
    logging.info(f'final output DF len={len(ndf)}')    
    return ndf


def sequence_value_counts(x):
    '''
    Function to be handed to Dask client for scalable calculation. 
    
    '''
    return x['sequence'].value_counts()



#
#    SPLIT MAPSEQ COLUMNS, FILTER             
#

def split_fields(df, column='sequence', drop=False, cp=None):
    '''
    E.g.
        vbc_read_st = 0
        vbc_read_end = 30
        spikeseq_st=24
        spikeseq_end = 32
        libtag_st=30
        libtag_end=32
        rtag_st=32
        rtag_end=37
        umi_st = 37
        umi_end = 49
        ssi_st = 49
        ssi_end = 57
    

    '''    
    logging.info(f'pulling out all fields defined in config.')
    if cp is None:
        cp = get_default_config()

    if not cp.has_section('split'):
        logging.error('This function requires config with [split] section...')
        sys.exit(1)

    plist = cp.items('split')
    fdict = defaultdict(dict)
    for (k,v) in plist:
        if k.endswith('_st'):
           fname =  k[:-3]
           fd = fdict[fname]
           fd['start'] = v
           
        elif k.endswith('_end'):
           fname =  k[:-4]
           fd = fdict[fname]
           fd['end'] = v
    
    for field in list(fdict.keys()):
        logging.debug(f'handling field={field}')
        try:
            fd = fdict[field]
            logging.debug(f'info for field={field} = {fdict[field]}')
            fstart = int( fdict[field]['start'])
            fend = int( fdict[field]['end'] )
            logging.debug(f'splitting field={field} start = {fstart} end = {fend}')
            df[field] = df[column].str.slice(fstart, fend).astype('string[pyarrow]')
            logging.debug(f'successfully split field={field}')
        
        except Exception as ex:
            logging.warning(f'problem with field={field} bad start/end? Check config?')
            
    sh = get_default_stats()
    
    if drop:
        logging.info(f'dropping {column} column to slim.')
        df.drop( column, axis=1, inplace=True)
    else:
        logging.info(f'drop is false. keeping {column} column.')
    logging.info(f'df done. len={len(df)} returning...')
    return df



def filter_reads_pd(df,
                    max_repeats=None,
                    max_n_bases=None,
                    min_reads=None, 
                    column='sequence',
                    remove=True,
                    cp=None):
    '''
    Used by filter_split.
    column should be nucleotide string. 

    '''
    if cp is None:
        cp = get_default_config()
    if min_reads is None:
        min_reads = int(cp.get('fastq','min_reads'))
    if max_repeats is None:
        max_repeats = int(cp.get('fastq','max_repeats')) 
    if max_n_bases is None:
        max_n_bases = int(cp.get('fastq','max_n_bases'))    

    logging.info(f'max_repeats = {max_repeats} max_n_bases={max_n_bases} min_reads={min_reads}')
    num_initial = str(len(df))
    
    sh = get_default_stats()
    
    # Filter by read_count 
    if min_reads > 1:
        logging.debug(f'Dropping reads with less than {min_reads} read_count.')
        logging.info(f'Length before read_count threshold={len(df)}')
        df = df[df['read_count'] >= min_reads]
        df.reset_index(inplace=True, drop=True)
        logging.info(f'Length after read_count threshold={len(df)}')    
    else:
        logging.info(f'min_reads = {min_reads} skipping initial read count thresholding.')  
    
    # Find Ns
    df['nb_valid'] = ~df[column].str.contains('N')
    num_has_n = str( len(df) - df['nb_valid'].sum() )
    
    # Find homopolymer runs
    df['nr_valid'] = True
    for nt in ['A','C','G','T']:
        rnt =  nt * (max_repeats + 1)
        logging.info(f'checking for {rnt}')
        rmap = ~df[column].str.contains(rnt)
        n_bad = len(df) - rmap.sum()
        logging.info(f'found {n_bad} max_repeat sequences for {nt}') 
        sh.add_value('/fastq',f'n_repeat_{nt}', str(n_bad) )
        df['nr_valid'] = df['nr_valid'] & rmap
        
    logging.debug(f"found {len(df) - df['nr_valid'].sum() }/{len(df)} bad sequences.")
    num_has_repeats = str( len(df) - df['nr_valid'].sum() )
    
    if remove:
        validmap = df['nr_valid'] & df['nb_valid']
        df = df[validmap]
        df.drop(['nb_valid','nr_valid'], inplace=True, axis=1)
        df.reset_index(drop=True, inplace=True)
        logging.info(f'remove=True, new len={len(df)}')
    logging.debug(f'types = {df.dtypes}')
    sh.add_value('/fastq_filter','num_initial', num_initial )
    sh.add_value('/fastq_filter','max_repeats', max_repeats )
    sh.add_value('/fastq_filter','num_has_repeats', num_has_repeats )
    sh.add_value('/fastq_filter','num_has_n', num_has_n )
    sh.add_value('/fastq_filter','num_kept', str(len(df)) )
    # changes made to inbound df, but return anyway
    return df


def filter_homopolymers(df, max_repeats=7, column='sequence', remove=True):
    '''
    find/remove Ns
    Find and remove homopolymer runs
    
    '''    
    # Find Ns
    df['nb_valid'] = ~df[column].str.contains('N')
    num_has_n = str( len(df) - df['nb_valid'].sum() )
    
    # Find homopolymer runs
    df['nr_valid'] = True
    for nt in ['A','C','G','T']:
        rnt =  nt * (max_repeats + 1)
        logging.info(f'checking for {rnt}')
        rmap = ~df[column].str.contains(rnt)
        n_bad = len(df) - rmap.sum()
        logging.info(f'found {n_bad} max_repeat sequences for {nt}') 
        df['nr_valid'] = df['nr_valid'] & rmap
        
    logging.debug(f"found {len(df) - df['nr_valid'].sum() }/{len(df)} bad sequences.")
    num_has_repeats = str( len(df) - df['nr_valid'].sum() )
    
    if remove:
        validmap = df['nr_valid'] & df['nb_valid']
        df = df[validmap]
        df.drop(['nb_valid','nr_valid'], inplace=True, axis=1)
        df.reset_index(drop=True, inplace=True)
        logging.info(f'remove=True, new len={len(df)}')
    logging.debug(f'types = {df.dtypes}')
    return df



def filter_fields(df,
                    drop=None,
                    cp=None):
    '''
    Used by filter_split.
    Remove reads that have bad fields. 
    
    Config values, e.g. 
        [readfilter]
        rtag_seq = GTACT
        rrtag_seq = CACGA

    will force 'rtag' column to contain literal string value GTACT or be dropped. 

    '''
    if cp is None:
        cp = get_default_config()

    if not cp.has_section('readfilter'):
        logging.error('This function requires config with [readfilter] section...')
        sys.exit(1)

    if drop is None:
        drop_mismatch = cp.getboolean('readfilter','drop_mismatch', fallback=False)

    sh = get_default_stats()

    plist = cp.items('readfilter')
    fdict = defaultdict(dict)
    
    for (k,v) in plist:
        logging.debug(f'handling {k}')
        if k.endswith('_seq'):
           fname =  k[:-4]
           fd = fdict[fname]
           fd['seq'] = v
           
        elif k.endswith('_regex'):
           fname =  k[:-6]
           fd = fdict[fname]
           fd['regex'] = v

    logging.debug(f'dict = {fdict}')
    num_initial = len(df)

    df['valid'] = True        
    df_cols = list( df.columns )
    for column in list( fdict.keys()):
        if column in df_cols:
            logging.info(f'handling field {column}')
            rdict = fdict[column]
            logging.debug(f'req dict for {column} is {rdict}')
            for reqtype in list( rdict.keys() ):
                logging.debug(f'handling requirement {reqtype}')
                if reqtype == 'seq':
                    rstr = rdict['seq']
                    logging.debug(f"checking {column} must match {rstr}  ")
                    vmap = df[column] == rstr
                    badmap = ~(vmap)
                    n_bad = len(df) - vmap.sum()
                    logging.info(f'found {n_bad} for {column}') 
                    sh.add_value('/field_filter',f'n_bad_{column}', str(n_bad) )
                    df['valid'] = df['valid'] & vmap
                if reqtype == 'regex':
                    logging.warning(f'regex filtering not implemented.')
        else:
            logging.warning(f'column {column} not in dataframe. Ignoring...')

    validmap = df['valid']
    n_bad = len(df) - validmap.sum()
    n_good = len(df) - n_bad
    df.drop(['valid'], inplace=True, axis=1)
    if drop_mismatch:
        df = df[validmap]
        df.reset_index(drop=True, inplace=True)
        logging.info(f'Removed bad, new len={len(df)}')
    else:
        logging.info('Not removing mismatches.')
    
    sh.add_value('/field_filter','num_initial', num_initial )
    sh.add_value('/field_filter','num_good', str(n_good)  )
    sh.add_value('/field_filter','num_bad',  str(n_bad )  )
    sh.add_value('/field_filter','num_kept', str(len(df)) )
    pct = ( n_bad / num_initial ) * 100
    spct = f'{pct:.2f}'
    sh.add_value('/field_filter','percent_bad', spct  ) 
       
    pct = ( n_good / num_initial ) * 100
    spct = f'{pct:.2f}'
    sh.add_value('/field_filter','percent_good', spct  )

    return df

#
#            READTABLE,  READ FREQUENCY PLOTS
#
def process_make_readtable_pd(df,
                              sampdf,
                              bcfile=None, 
                              outdir=None, 
                              cp = None):
    '''
    take split/aligned read file and produce a produce 
    fully tagged and read-oriented data table, with 
    validated/consistent contents
    
    -- convert SSI column to BCXXX tag. 
    -- add optional region label
    -- classify by site-type
    -- set brain label
    -- aggregate and sum read counts 
    -- remove redudant info (SSI, libtag, spikeseg)

    '''        
    logging.info(f'inbound df len={len(df)} columns={list( df.columns )}')
    sh = get_default_stats()
    
    n_initial = len(df)
    if outdir is None:
        outdir = os.path.abspath('./')
    outdir = os.path.abspath(outdir)    
    os.makedirs(outdir, exist_ok=True)
    
    if cp is None:
        cp = get_default_config()

    project_id = cp.get('project','project_id')
    spikeseq = cp.get('readtable','spikeseq')
    realregex = cp.get('readtable', 'realregex' )
    loneregex = cp.get('readtable', 'loneregex' )    

    # establish logic for using libtag, if defined. 
    # optionally keep/remove L1s
    use_libtag = cp.getboolean('readtable','use_libtag', fallback=True)
    filter_by_libtag = cp.getboolean('readtable','filter_by_libtag', fallback=True)
    use_lone = cp.getboolean('readtable','use_lone', fallback=True)
    include_lone = cp.getboolean('vbcfilter','include_lone', fallback=False )
    if not use_libtag:
        filter_by_libtag = False
    
    if bcfile is None:
        bcfile = os.path.expanduser( cp.get('barcodes','ssifile') )    
    
    logging.debug(f'spikeseq={spikeseq} realregex={realregex} loneregex={loneregex} bcfile={bcfile} use_lones={use_libtag}')

    # Map label, rtprimer to SSIs    
    logging.debug(f'getting rt labels...')
    labels = get_rtlist(sampdf)
    logging.debug(f'rtlabels={labels}')
    bcdict = get_barcode_dict(bcfile, labels)
    rtdict = get_rtprimer_dict(bcfile, labels)
    rtbcdict = get_rtbc_dict(bcfile, labels)
    logging.debug(f'got {len(bcdict)} barcodes with labels and primer number.')    
   
    logging.info('filling in rtprimer number by SSI sequence...')    
    df['rtprimer'] = df['ssi'].map(rtdict)
    df['rtprimer'] = df['rtprimer'].astype('category')
    
    logging.info('filling in label by rtprimer...')
    df['label'] = df['rtprimer'].map(rtbcdict)
    df['label'] = df['label'].astype('category')

    # identify and remove unmatched SSI sequences? 
    logging.debug('removing unmatched SSI')
    badssidf = df[ df.isna().any(axis=1) ]
    n_badssi = len(badssidf)
    df.drop(badssidf.index, inplace=True)
    df.reset_index(inplace=True, drop=True)    

    # We don't need the ssi sequence once label/rtprimer are set and we've
    # removed and written the bad ones.  
    logging.debug('dropping SSI sequence')
    df.drop('ssi', inplace=True, axis=1)

    logging.debug('writing out bad_ssi.')
    of = os.path.join( outdir, 'bad_ssi.tsv')
    badssidf.reset_index(inplace=True, drop=True)
    start = dt.datetime.now()    
    badssidf.to_csv(of, sep='\t')
    end = dt.datetime.now()
    delta_seconds = (dt.datetime.now() - start).seconds
    log_transferinfo(of, delta_seconds)    
    logging.info(f'Wrote bad_ssi DF len={len(badssidf)} to {of}')
    badssidf = None   

    # set site from rtprimer which should now be correct.      
    sdf = sampdf[['rtprimer','siteinfo']]
    sdf = sdf[sdf['siteinfo'] != '']
    smap = dict(zip(sdf['rtprimer'], sdf['siteinfo']))
    df['site'] = df['rtprimer'].map(smap)


    if use_libtag:
        # calculate and save libtag abundance.
        of = os.path.join( outdir, 'libtag_counts.tsv') 
        lcdf = df['libtag'].value_counts().reset_index()
        lcdf.to_csv(of, sep='\t')
        logging.info(f'Saved libtag counts to {of}')
        
        logging.info(f'use_libtag = True. Using libtag to determine L2 type.')
        logging.info('identifying L2 (reals) by libtag...')
        rmap = df['libtag'].str.match(realregex)
        df.loc[rmap, 'type'] = 'real'

        # must be set after identifying reals. all spikes are reals, but
        # not all reals are spikes. 
        logging.info(f'Identifying spikeins by spikeseq={spikeseq}')
        smap = df['spikeseq'] == spikeseq
        df.loc[smap, 'type'] = 'spike'

        if use_lone:        
            logging.info('identifying L1s by libtag...')
            lmap = df['libtag'].str.match(loneregex)
            df.loc[lmap, 'type'] = 'lone'

            # find and remove (at least) known template-switch rows from dataframe.
            # template switch type is valid L1/lone (from libtag) but is a valid target (from SSI)
            # OR target-lone site but NOT a valid L1 type from libtag.   
            nonlone = df[ ~df['site'].str.startswith('target-lone') ]
            nonlone_tsdf = nonlone[ nonlone['type'] == 'lone' ]
            
            lone = df[ df['site'].str.startswith('target-lone') ]
            lone_tsdf = lone[ lone['type'] != 'lone' ]

            tsdf = pd.concat( [nonlone_tsdf, lone_tsdf])
            # tsdf = nonlone[ ((nonlone['type'] == 'lone') & ( nonlone['site'].str.startswith('target'))) ]

            of = os.path.join(outdir, f'{project_id}.template_switch.tsv') 
            if len(tsdf) > 0:
                logging.info(f'Writing template switch DF len={len(tsdf)} Writing to {of}')
                tsdf.to_csv(of, sep='\t')
            n_tswitch = len(tsdf)
            
            # remove known template switch from readtable
            df.drop(nonlone_tsdf.index, inplace=True)
            df.drop(lone_tsdf.index, inplace=True)
            df.reset_index(drop=True, inplace=True)
            sh.add_value('/readtable', 'n_tswitch', str(n_tswitch) )
            
            # remove and save epected lones.
            lones = df[ df['type'] == 'lone']
            if not include_lone:
                df = df [df['type']  != 'lone']
            of = os.path.join(outdir, f'{project_id}.valid_lones.tsv') 
            lones.to_csv(of, sep='\t')
            n_lones = len(lones) 
            sh.add_value('/readtable', 'n_lones', str(n_lones) )

        if filter_by_libtag:
            logging.debug('Filtering by bad libtag/type')  
            badtypedf = df[ df.isna().any(axis=1) ]
            n_badtype = len(badtypedf)
            df.drop(badtypedf.index, inplace=True)
            df.reset_index(inplace=True, drop=True)    
    
            of = os.path.join( outdir, 'bad_type.tsv')
            badtypedf.reset_index(inplace=True, drop=True)
            badtypedf.to_csv(of, sep='\t')
            logging.info(f'Wrote bad_type DF len={len(badtypedf)} to {of}')
            badtypedf = None   
            sh.add_value('/readtable', 'n_badtype', str(n_badtype) )

        else:
                logging.info('No filtering by libtag. Setting unidentified to real.')
                namap = df['type'].isna()
                df.loc[namap, 'type'] = 'real'            

    else:
        logging.info('use_libtag = False. Ignoring libtag. Just using spike sequence. ')
        logging.info(f'Identifying spikeins by spikeseq={spikeseq}')
        smap = df['spikeseq'] == spikeseq
        df.loc[smap, 'type'] = 'spike'

        logging.info('ignoring L1/L2 libtag. All non-spikes are real.')        
        nsmap = df['type'] != 'spike'
        df.loc[nsmap, 'type'] = 'real'

    logging.debug('Dropping redundant spikeseq field.')
    #df.drop(['spikeseq','libtag'], inplace=True, axis=1)
    df.drop(['spikeseq'], inplace=True, axis=1)
    
    # NOT REQUIRED VALUES, so replace NaNs 
    # set brain
    bdf = sampdf[['rtprimer','brain']]
    bdf = bdf[bdf['brain'] != '']
    bmap = dict(zip(bdf['rtprimer'],bdf['brain']))
    df['brain'] = df['rtprimer'].map(bmap)
    
    try:
        df['brain'] = df['brain'].cat.add_categories([''])
    except:
        logging.debug('brain not categorical?')
           
    df.fillna({'brain': ''}, inplace=True)
    bdf = None
    
    # set region
    rdf = sampdf[['rtprimer','region']]
    rdf = rdf[rdf['region'] != '']
    rmap = dict(zip(rdf['rtprimer'],rdf['region']))
    df['region'] = df['rtprimer'].map(rmap)
    
    try:
        df['region'] = df['region'].cat.add_categories([''])
    except:
        logging.debug('region not categorical?')
    df.fillna({'region': ''}, inplace=True)
    rdf = None

    # set ourtube
    tdf = sampdf[['rtprimer','ourtube']]
    tdf = tdf[tdf['ourtube'] != '']
    tmap = dict(zip(tdf['rtprimer'],tdf['ourtube']))
    df['ourtube'] = df['rtprimer'].map(tmap)
    try:
        df['ourtube'] = df['ourtube'].cat.add_categories([''])
    except:
        logging.debug('ourtube not categorical?')
    df.fillna({'ourtube': ''}, inplace=True)
    tdf = None
        
    sh.add_value('/readtable','n_full_sequences', str(len(df)) )    
    n_final = len(df)

    sh.add_value('/readtable', 'n_initial', str(n_initial) ) 
    sh.add_value('/readtable', 'n_badssi', str(n_badssi) )
    sh.add_value('/readtable', 'n_final', str(n_final) )     
    
    df = fix_mapseq_df_types(df, fformat='readtable')
    return df


def threshold_by_sample(df, sampdf):
    '''
    min_reads by sample from sampleinfo
    
    '''
    logging.debug(f'thresholding per sample: initial size={len(df)}')
    outdf = None
    sampmap = list( zip(sampdf['rtprimer'], sampdf['min_reads']))
    for (rtprimer, min_reads) in sampmap:
        min_reads = int(min_reads)
        logging.debug(f'thresholding rtprimer {rtprimer} to min_reads={min_reads}')
        rtdf = df[df['rtprimer'] == rtprimer]
        rtdf = rtdf[rtdf['read_count'] >= min_reads ]
        if outdf is None:
            outdf = rtdf
        else:
            outdf = pd.concat([outdf, rtdf], ignore_index=True)
    outdf.reset_index(inplace=True, drop=True)
    logging.debug(f'thresholding per sample: final size={len(outdf)}')    
    return outdf


#
#    ALIGN, COLLAPSE BY VBC             
#
def align_collapse(df,
                      column='vbc_read',
                      pcolumn='read_count',
                      gcolumn=None,  
                      max_mismatch=None,
                      max_recursion=None, 
                      outdir=None, 
                      datestr=None,
                      force=False,
                      min_reads = None,
                      drop_pcolumn = True, 
                      cp=None):
    '''
    Entry point for CLI align_collapse. 
    Determine grouped vs. non-grouped. 
    
    @arg column         Column to align and collapse
    @arg pcolumn        Column to decide which sequence to assign to all. 
    @arg gcolumn        Column to group processing within
    @arg max_mismatch   Max Hamming distance between collapsed sequences 
    @arg max_recursion  Max recursion to set for Python interpreter. 
    @arg outdir         Base outdir for intermediate output. 
    @arg datestr        Optional date string for logging
    @arg force          Recalculate intermediate steps, otherwise pick up.
    @arg min_reads      Optional thresholding before collapse 
    @arg cp             ConfigParser object with [collapse] section.
    
    '''
    sh = get_default_stats()
    
    collapse_lib = cp.get('collapse','library', fallback = 'networkx')

    if max_mismatch is None:
        max_mismatch = cp.get('collapse','max_mismatch')
    
    
    if gcolumn is not None:
        logging.info(f'Grouped align_collapse. Group column = {gcolumn}')
        sh.add_value(f'/collapse','n_initial', str(len(df)) )
        sh.add_value(f'/collapse','collapse_mode', 'grouped' )
        sh.add_value(f'/collapse','collapse_library', collapse_lib )
        sh.add_value(f'/collapse','cli_max_mismatch', str(max_mismatch) )
        
        logging.info(f'Grouped align_collapse. Saving rows with no group.')
        ngdf = df[ ( df[gcolumn].isna() ) | ( df[gcolumn] == '' )  ]
        sh.add_value(f'/collapse','n_nogroup', str(len(ngdf)) )
        
        df = df.drop(ngdf.index)
        df.reset_index(inplace=True, drop=True)
        sh.add_value(f'/collapse','n_withgroup', str(len(df)) )
        
        if collapse_lib == 'networkx':
            logging.debug(f'networkx collapse...')
            if max_mismatch is None:
                max_mismatch = cp.get('collapse','max_distance')
            df = align_collapse_nx_grouped( df, 
                                        column = column,
                                        pcolumn = pcolumn,
                                        gcolumn = gcolumn,
                                        max_distance=max_mismatch,
                                        max_recursion=max_recursion, 
                                        outdir=outdir, 
                                        datestr=datestr,
                                        force=force,
                                        min_reads = min_reads,
                                        cp=cp )            
        elif collapse_lib == 'custom':
            logging.debug(f'custom collapse...')
            df = align_collapse_pd_grouped( df, 
                                        column = column,
                                        pcolumn = pcolumn,
                                        gcolumn = gcolumn,
                                        max_mismatch=max_mismatch,
                                        max_recursion=max_recursion, 
                                        outdir=outdir, 
                                        datestr=datestr,
                                        force=force,
                                        min_reads = min_reads,
                                        cp=cp )
        logging.info(f'merging back saved nogroup rows. ')
        ngdf[f'{column}_col'] = ngdf[column]
        df = pd.concat( [ df , ngdf ], ignore_index=True )
        df.reset_index(inplace=True, drop=True)
        sh.add_value(f'/collapse','n_remerged', str(len(df)) )
    
    else:
        logging.info(f'Global align_collapse.')
        sh.add_value(f'/collapse','collapse_mode', 'global' )
        sh.add_value(f'/collapse','collapse_library', collapse_lib )
        sh.add_value(f'/collapse','cli_max_mismatch', str(max_mismatch) )
        df = align_collapse_pd( df = df,
                                column = column,
                                pcolumn = pcolumn,
                                max_mismatch=max_mismatch,
                                max_recursion=max_recursion, 
                                outdir=outdir, 
                                datestr=datestr,
                                force=force, 
                                min_reads = min_reads,
                                cp=cp ) 
    logging.info(f'Got DF len={len(df)} Fixing dtypes...')
    try:
        df = fix_mapseq_df_types(df, fformat='readtable')
    except KeyError:
        logging.warning(f'cannot fix types. May not be readtable format. ')

    logging.info('Renaming collapsed column(s)')
    df.rename(columns={f'{column}' : f'{column}_orig'}, inplace=True)
    df.rename(columns={f'{column}_col' : f'{column}'}, inplace=True)

    logging.info(f'Returning DF len={len(df)} dtypes={df.dtypes}')
    return df
                            

#
#            VBCTABLE
#
def process_make_vbctable_pd(df,
                          outdir=None,
                          inj_min_reads = 2,
                          target_min_reads = 2,
                          gcolumn = 'vbc_read',
                          sampdf= None, 
                          cp = None):
    '''   
    
    primary columns:  label, umi, type  
    derived columns:  brain, region, rtprimer, label, ssi
    -- remove nomatches
    -- threshold read_count by site type
    -- collapse by VBC sequence, calculating UMI count
    
    Drops source, and possibly other fields since aggregation would make them ambiguous.
    
    This should  be prepared for input it make VBC matrices    
    '''
    if cp is None:
        cp = get_default_config()

    project_id = cp.get('project','project_id')
    use_lone = cp.getboolean('readtable','use_lone')
    filter_by_libag = cp.getboolean('readtable','filter_by_libtag')
    
    logging.debug(f'inbound df len={len(df)} columns={df.columns}')
    log_objectinfo(df, 'input-df')
    ndf = df
    
    ndf.replace('nomatch', np.nan, inplace=True)
    logging.info(f'DF before removing nomatch: {len(ndf)}')
    log_objectinfo(ndf, 'new-df')
    ndf.dropna(inplace=True, axis=0, ignore_index=True)
    logging.info(f'DF after removing nomatch/NaN: {len(ndf)}')
    log_objectinfo(ndf, 'new-df')

    # create spike-in QC tables before thresholding. 
    of = os.path.join( outdir, f'{project_id}.spikeqc.prethresh.xlsx' ) 
    make_spikein_qctable( ndf,
                          outfile = of,
                          cp = cp ) 

    # perform optional per-sample read thresholding.
    if sampdf is not None:
        logging.debug('sampleinfo DF provided....')
        sampdf['min_reads'] = sampdf['min_reads'].astype(int)
        if int( sampdf['min_reads'].max()) > 1:
            logging.info('performing per-sample read count thresholding...')
            ndf = threshold_by_sample(ndf, sampdf)
            ndf.reset_index(drop=True, inplace=True)
    else:
        logging.debug('sampleinfo DF not provided.')
    
    # threshold by read counts, by siteinfo, before starting...
    logging.info(f'thresholding by read count. inj_min_reads={inj_min_reads} target_min_reads={target_min_reads} len={len(ndf)}') 
    tdf = ndf[ndf['site'].str.startswith('target')]
    tdf = tdf[tdf['read_count'] >= int(target_min_reads)]
    idf = ndf[ndf['site'].str.startswith('injection')]
    idf = idf[idf['read_count'] >= int(inj_min_reads)]
    thdf = pd.concat([tdf, idf])
    thdf.reset_index(drop=True, inplace=True)

    # create spike-in QC tables after thresholding. 
    of = os.path.join( outdir, f'{project_id}.spikeqc.postthresh.xlsx' ) 
    make_spikein_qctable( thdf,
                          outfile = of,
                          cp = cp )
    
    logging.info(f'DF after threshold inj_min_reads={inj_min_reads} target_min_reads={target_min_reads}: {len(thdf)}')
    log_objectinfo(thdf, 'threshold-df')    
    
    agg_params = {  'umi'       : 'nunique',
                    'read_count': 'sum', 
                    'brain'     : 'first',
                    'region'    : 'first',
                    'site'      : 'first',
                    'ourtube'   : 'first',    
                }
   
    udf = thdf.groupby([ gcolumn, 'label', 'type'], observed=True ).agg( agg_params ).reset_index()
    udf.rename( {'umi':'umi_count'}, axis=1, inplace=True)
    logging.info(f'DF after umi/label collapse: {len(udf)}')

    # output controls by SSI/site, save to TSV
    controls = udf[ udf['site'].isin( CONTROL_SITES ) ]
    controls.reset_index(inplace=True, drop=True)
    controls.to_csv(f'{outdir}/{project_id}.controls.tsv', sep='\t')

    udf.reset_index(inplace=True, drop=True)  
    sh = get_default_stats()
    sh.add_value('/vbctable','n_vbcs', len(udf) )        
    log_objectinfo(udf, 'umi-df')
    
    udf = fix_mapseq_df_types(udf, fformat='vbctable')
    return udf


#
#        VBCFILTER
#
def process_filter_vbctable(df, 
                               inj_min_umi = None,
                               target_min_umi = None,
                               target_min_umi_absolute = None,
                               outdir = None,
                               cp=None):
    '''
    Take unfiltered vbctable. 
    apply all thresholds/filters for input to matrices. 
    
    Per brain...
    -- <type>_min_umi against reals (but not spikes).
    
    threshold logic. 
    inj_min_umi                VBC UMI must exceed to be kept.
    target_min_umi             if ANY target area exceeds, keep all of that VBC targets. 
    target_min_umi_absolute    hard threshold cutoff
    
    How to handle controls when include_controls=True?
    Which controls to return for reporting?
    
    '''
    if cp is None:
        cp = get_default_config()

    if outdir is None:
        outdir = os.path.abspath('./')
        
    if inj_min_umi is None:
        inj_min_umi = int(cp.get('vbcfilter','inj_min_umi'))
    if target_min_umi is None:
        target_min_umi = int(cp.get('vbcfilter','target_min_umi'))   
    if target_min_umi_absolute is None:
        target_min_umi_absolute = int(cp.get('vbcfilter','target_min_umi_absolute'))

    project_id = cp.get('project','project_id')
    require_injection = cp.getboolean('vbcfilter','require_injection')
    include_injection = cp.getboolean('vbcfilter','include_injection')
    include_controls = cp.getboolean('vbcfilter','include_controls')
    include_lone = cp.getboolean('vbcfilter','include_lone', fallback=False )
    use_target_negative=cp.getboolean('vbcfilter','use_target_negative')
    use_target_water_control=cp.getboolean('vbcfilter','use_target_water_control')

    logging.debug(f'filtering VBC table.  n_rows={len(df)} outdir={outdir} inj_min_umi={inj_min_umi} target_min_umi={target_min_umi} ')

    df['brain'] = df['brain'].astype('string')
    bidlist = list(df['brain'].dropna().unique())
    bidlist = [ x for x in bidlist if len(x) > 0 ]
    bidlist.sort()
    logging.debug(f'handling brain list: {bidlist}')
    
    sh = get_default_stats()

    # optionally keep/remove for inclusion in each brain matrix.
    controls = df[ df['site'].isin( CONTROL_SITES ) ]

    # optionally keep/remove L1s if not removed prior to now. 
    # include_lone = cp.getboolean('vbcfilter','include_lone', fallback=False )
    if not include_lone:
        logging.info('include_lone is False. Removing L1s...')
        controls = controls[ controls['site'] != 'target-lone']

    # save for reference
    if len(controls) > 0:
        controls.reset_index(inplace=True, drop=True)
        controls.to_csv(f'{outdir}/{project_id}.controls.tsv', sep='\t')

    norm_dict = {}

    bdflist = []

    for brain_id in bidlist:
        valid = True
        ndf = None 
        logging.debug(f'handling brain_id={brain_id}')
        bdf = df[df['brain'] == brain_id ]
        
        if include_controls:
            logging.info(f'include_controls=True, adding controls to brain.')
            logging.debug(f'before merging: len(controls) = {len(controls)} len(bdf)={len(bdf)}')
            bcdf = controls.copy()
            bcdf['brain'] = brain_id
            bdf = pd.concat( [bdf, bcdf], ignore_index = True)
            logging.debug(f'after merging: len(bdf) = {len(bdf)}') 
        bdf.reset_index(inplace=True, drop=True)
        initial_len = len(bdf)  
        logging.info(f'[{brain_id}] full initial len={len(bdf)}')
        
        # defaults for using controls/target-negatives
        max_negative = 1
        max_water_control = 1
        
        logging.debug(f'[{brain_id}]: inj_min_umi={inj_min_umi} target_min_umi={target_min_umi}')
        logging.debug(f'[{brain_id}]:use_target_negative={use_target_negative} use_target_water_control={use_target_water_control}')
        sh = get_default_stats() 

        # remove spikes and save them 
        # Spikes do NOT get thresholded by UMI, or restricted by injection/target presence 
        spikes = bdf[ bdf['type'] == 'spike']
        reals = bdf[ bdf ['type'] != 'spike']
    
        # Separate targets and injections for specific UMI thresholding. 
        targets = reals[ reals['site'].str.startswith('target') ]
        injection = reals[ reals['site'].str.startswith('injection') ]    
        logging.info(f'[{brain_id}]: inputs raw={initial_len} real_inj={len(injection)} real_targets={len(targets)} spikes={len(spikes)} controls={len(controls)} ')
    
        # threshold injection(s)
        if inj_min_umi > 1:
            before = len(injection)            
            injection = injection[ injection['umi_count'] >= inj_min_umi ]
            injection.reset_index(inplace=True, drop=True)
            logging.debug(f'[{brain_id}]: filtering by inj_min_umi={inj_min_umi} before={before} after={len(injection)}')
        else:
            logging.debug(f'[{brain_id}]:  inj_min_umi={inj_min_umi} no filtering needed.')
        logging.info(f'[{brain_id}]: after inj_min_umi real_inj={len(injection)} real_targets={len(targets)} spikes={len(spikes)} controls={len(controls)} ')
    
        # threshold target(s)
        # target_min_umi             if any target exceeds, include for all targets
        # target_min_umi_absolute    unconditionally threshold target*
        # inj_min_umi                unconditionally threshold injection
    
        # apply absolute if relevant    
        if target_min_umi_absolute > 1:
            before = len(targets)
            targets = targets[targets['umi_count'] >= target_min_umi_absolute ]
            targets.reset_index(drop=True, inplace=True)
            after = len(targets)
            logging.debug(f'[{brain_id}]: filtering by target_min_umi_absolute={target_min_umi_absolute} before={before} after={after}')
        #
        # handle target thresholding for non-absolute case.     
        # threshold by target_min_umi or threshold by target-negative
        # if use_target_negative is true, but no target negative site 
        # defined, use target_min_umi and throw warning. 
        if use_target_negative:
            logging.info(f'[{brain_id}]: use_target_negative is {use_target_negative}')
            max_negative = calc_min_umi_threshold(targets, 'target-negative', cp)
            logging.debug(f'[{brain_id}]: target-negative UMI count = {max_negative}')
    
        if use_target_water_control:
            logging.info(f'[{brain_id}]: use_target_water_control is {use_target_water_control}')
            max_water_control = calc_min_umi_threshold(targets, 'target-water-control',cp)
            logging.debug(f'[{brain_id}]: target_water_control UMI count = {max_water_control}')
    
        target_min_umi = max([target_min_umi, max_negative, max_water_control ])
        logging.debug(f'[{brain_id}]: min_target UMI count after all constraints = {target_min_umi}')   
    
        # threshold by target_min_umi but if any target exceeds, keep all        
        if target_min_umi > 1:
            before = len(targets)
            targets = filter_targets_min_umi_any(targets, min_umi = target_min_umi)
            after = len(targets)
            logging.info(f'[{brain_id}]: before={before} -> {after} filtering for VBC any target >= {target_min_umi}')
            
        # get injection-filtered real target table, and target-filtered real injection table
        # in case either is needed. 
        (ftargets, finjection ) = filter_non_injection(targets, injection)            
        logging.debug(f'[{brain_id}]: {len(ftargets)} real target VBCs after injection filtering.')
        logging.debug(f'[{brain_id}]: {len(finjection)} real injection VBCs after target filtering.')
    
        if require_injection:
            logging.debug(f'[{brain_id}]: require_injection={require_injection} inj_min_umi={inj_min_umi}')
            targets = ftargets
            if len(targets) < 1:
                logging.warning(f'[{brain_id}]: brain {brain_id} did not have any targets passing injection filtering.')
        else:
            logging.debug(f'[{brain_id}]: require_injection={require_injection} proceeding...')
        
        if len(targets) > 1 : 
            if include_injection:
                    logging.debug(f'[{brain_id}]: include_injection={include_injection} Merging {len(targets)} targets + {len(injection)} injection VBCs')
                    ndf = pd.concat( [ targets, finjection], ignore_index=True )
                    ndf.reset_index(inplace=True, drop=True) 
            else:
                logging.debug(f'[{brain_id}]: include_injection={include_injection} excluding injection VBCs from table.')
                ndf = targets
        else:
            logging.debug(f'[{brain_id}]: include_injection={include_injection} but no targets passed. Output empty targets.')
            ndf = targets
        
        # re-create df with un-thresholded spike-ins
        # ndf may be empty
        if len(ndf) > 1:               
            ndf = pd.concat([spikes, ndf ], ignore_index=True)   
            ndf.reset_index(inplace=True, drop=True)
            logging.debug(f'[{brain_id}]: brain {brain_id} len={len(bdf)} output DF:\n{ndf}')
            logging.info(f'[{brain_id}] final len={len(ndf)} appending to bdflist.')
            bdflist.append(ndf)
        else:
            logging.info(f'[{brain_id}]: No targets passed filtering. Not including.')
    
    # merge all brains into one dataframe...
    logging.debug(f'sizes={[ len(x) for x in bdflist ]} ')

    if len(bdflist) > 0:    
        df = pd.concat(bdflist, ignore_index = True)
        df.reset_index(inplace=True, drop=True)
        logging.info(f'All brains. Final merged filtered DF len={len(df)}')
        df = fix_mapseq_df_types(df, fformat='vbctable')
    else:
        logging.warning('No data passed filtering. Creating empty dataframe.')
        df = pd.DataFrame(columns=df.columns)
    return ( df, controls )


def filter_targets_min_umi_any(targets, min_umi ):
    '''
    retain all VBC rows for which *any* target label umi_count exceeds target_min_umi
    remove all others. 
    '''
    logging.debug(f'inbound targets len={len(targets)} min_umi={min_umi}')
    vmax = targets.groupby('vbc_read').agg({'umi_count': 'max' }).copy()
    vmax.reset_index(inplace=True, drop=False)
    init_len = len(vmax)
    vmax = vmax[vmax['umi_count'] >= min_umi ]
    vmax.drop('umi_count', inplace=True, axis=1)
    vmax.reset_index(inplace=True, drop=True)
    logging.debug(f'unique VBC by max: before={init_len} after={len(vmax)}')

    result = vmax.merge(targets, on='vbc_read', how='left')
    #result.drop('umi_count_x', inplace=True, axis=1)
    #result.rename( {'umi_count_y':'umi_count'}, inplace=True, axis=1)
    result.reset_index(inplace=True, drop=True)
    
    # check again..
    vmax = result.groupby('vbc_read').agg({'umi_count': 'max' })
    minmax = vmax.min()
    logging.debug(f'result len={len(result)} new smallest max per VBC = {minmax}')
    return result
   


def filter_non_injection(rtdf, ridf, write_out=False):
    '''
    rtdf and ridf should already be filtered by brain, type, and minimum umi_count
    remove rows from rtdf that do not also have VBC in ridf with the same sequience
    Does an inner join() on the dataframes, keyed on sequence. 
    Keeps values and columns from first argument (rtdf)
    Also filters injection df by presence in target, for possible inclusion.
    
    '''
    logging.info(f'filtering non-injection. ')
    logging.debug(f'before threshold inj df len={len(ridf)}')
    
    # get target VBCs that are in injection
    mtdf = merge_and_filter(rtdf, ridf)
    mtdf = fix_mapseq_df_types(mtdf, fformat='vbctable')

    # get injection VBCs that are in at least one target, similarly 
    midf = merge_and_filter(ridf, rtdf)
    midf = fix_mapseq_df_types(midf, fformat='vbctable')
    return ( mtdf, midf)


def merge_and_filter(adf, bdf, on='vbc_read', indicator=True, how='outer'):
    '''
    Return filtered ADF consisting of only entries in ADF that also exist in BDF as joined on 
    the 'on' column.  
    
    '''
    #suffixes=('_a','_b')
    # get targets in injection merge using vbc_read as common field. 
    #mdf = rtdf.merge(ridf, on='vbc_read', indicator=True, how='outer', suffixes=('_t','_i'))
    mdf = adf.merge(bdf, on='vbc_read', indicator=True, how='outer', suffixes=('_a','_b'))
    
    # only keep target VBCs that are in injection
    mdf = mdf[mdf['_merge'] == 'both']    
    incol = mdf.columns
    outcol = []
    selcol =[]
    for c in incol:
        if not c.endswith('_b'):
            # _t or join column. 
            selcol.append(c)
            outcol.append(c.replace('_a',''))
    mdf = mdf[selcol]
    mdf.columns = outcol
    logging.debug(f'before drop_duplicates. len={len(mdf)}')
    mdf.drop_duplicates(inplace=True)
    mdf.drop('_merge', inplace=True, axis=1)
    mdf.reset_index(drop=True, inplace=True)
    logging.debug(f'after drop_duplicates. len={len(mdf)} columns={mdf.columns}')
    logging.debug(f'created merged/joined DF w/ common sequence items.  df=\n{mdf}')
    return mdf


def filter_all_lt(df, key_col='sequence', val_col='umi_count', threshold=5):
    '''
    filter dataframe where *all* are less than a threshold. 
    
    takes df, groups by key_col and sums val_col. 
    keeps *all* rows for keys where *any* val_col > threshold for that key.  
    '''
    logging.debug(f'inbound df len={len(df)}')
    maxdf = df.groupby(key_col)[val_col].max()
    logging.debug(f'maxdf: unique  key_col len={len(maxdf)}')
    fmdf = maxdf[maxdf > threshold]
    logging.debug(f'fmdf: unique keys that pass threshold len={len(fmdf)}')
    keeplist = list(fmdf.index)
    logging.debug(f'keep list len={len(keeplist)}')
    outdf = df[df[key_col].isin(keeplist)]
    outdf.reset_index(inplace=True, drop=True)
    logging.debug(f'outdf len={len(outdf)}')
    return outdf


def calc_min_umi_threshold(df, site='target-negative', cp=None):
    '''
    retrieve maximum umi_count value for provided site type. 
    returns max max value...

    '''
    if cp is None:
        cp = get_default_config()
    
    logging.debug(f'getting UMI threshold type={site}')
    countlist = []
    min_threshold = 0
    tndf = df[ df['site'] == site]
    lablist = list(tndf['label'].dropna().unique())
    for label in lablist:
        ldf = tndf[tndf['label'] == label]
        if len(ldf) > 0:
            maxumi = ldf['umi_count'].max()
            logging.debug(f'max umi_count for label {label} = {maxumi}')
            countlist.append( maxumi )
    if len(countlist) > 0:
        min_threshold = max(countlist)
        logging.debug(f'calculated UMI threshold={min_threshold} type={site} ')    
    return min_threshold

#
#    MATRICES
#
def process_make_matrices(df,
                          sampdf=None, 
                          outdir=None,
                          exp_id=None,  
                          label_column='label',
                          cp = None):
    '''
    
    Simplified matrix creation. 
    NOTE: Assumes ALL input is valid, thresholded, and to be included.
    
    -- per brain, pivot real VBCs (value=umi counts)
    -- create real, real normalized by spike-in, and scaled normalized.   
    -- use label_column to pivot on, making it the y-axis, x-axis is vbc sequence.  
    -- calculate and export spike-in counts. 
    
    '''
    if cp is None:
        cp = get_default_config()
    if outdir is None:
        outdir = os.path.abspath('./')
    if exp_id is None:
        exp_id = cp.get('project','project_id')

    spikein_ratio_map = None
    if sampdf is not None:
        # get si_ratio map for BCXXX columns
        sirmap = get_si_ratio_map(sampdf, column=label_column)
        logging.info(f'spikein ratio map exists. Checking for any != 1.0...')
        for k in sirmap.keys():
            v = sirmap[k]
            if v > 1.0:
                spikein_ratio_map = sirmap
                logging.debug(f'found non-1.0 spike-in ratio. setting for use.')
                break
    logging.debug(f'spikein_ratio_map = {spikein_ratio_map}')

    logging.debug(f'running exp_id={exp_id} ')

    df['brain'] = df['brain'].astype('string')
    bidlist = list(df['brain'].dropna().unique())
    bidlist = [ x for x in bidlist if len(x) > 0 ]
    bidlist.sort()
    logging.debug(f'handling brain list: {bidlist}')
    
    sh = get_default_stats()

    norm_dict = {}
    spike_dict = {}
    real_dict = {}

    for brain_id in bidlist:
        valid = True 
        logging.debug(f'handling brain_id={brain_id}')
        bdf = df[df['brain'] == brain_id]                               
        bdf.to_csv(f'{outdir}/{brain_id}.vbcdata.tsv', sep='\t')
        bdf.to_parquet(f'{outdir}/{brain_id}.vbcdata.parquet')

        # separate reals and spikes
        reals = bdf[bdf['type'] == 'real']       
        ridf = reals[reals['site'].str.startswith('injection')]
        rtdf = reals[reals['site'].str.startswith('target')]

        spikes = bdf[bdf['type'] == 'spike']
        stdf = spikes[spikes['site'].str.startswith('target')]
        sidf = spikes[spikes['site'].str.startswith('injection')]
        
        target_columns = list( rtdf[label_column].unique())
        injection_columns = list( ridf[label_column].unique())

        # reals        
        rbcmdf = reals.pivot(index='vbc_read', columns=label_column, values='umi_count')
        scol = natsorted(list(rbcmdf.columns))
        rbcmdf = rbcmdf[scol]
        rbcmdf.fillna(value=0, inplace=True)
        logging.debug(f'brain={brain_id} raw real barcode matrix len={len(rbcmdf)}')            
        if not len(rbcmdf) > 0:  
            valid = False
            logging.warning(f'brain={brain_id} not valid.')
        
        # spikes    
        sbcmdf = spikes.pivot(index='vbc_read', columns=label_column, values='umi_count')
        spcol = natsorted(list(sbcmdf.columns))
        sbcmdf = sbcmdf[spcol]
        sbcmdf.fillna(value=0, inplace=True)    
        logging.debug(f'brain={brain_id} spike barcode matrix len={len(sbcmdf)}')
        if not len(sbcmdf) > 0:  
            valid = False
            logging.warning(f'brain={brain_id} not valid.')

        if valid:
            (rbcmdf, sbcmdf) = sync_columns(rbcmdf, sbcmdf)
            nbcmdf = normalize_weight_grouped(rbcmdf, sbcmdf, columns = [target_columns, injection_columns])
            
            if spikein_ratio_map is not None:
                logging.info(f'spikein ratio map exists. Using.')
                n_sicols = 0
                for col in nbcmdf.columns:
                    sir = spikein_ratio_map[col]
                    if sir != 1.0:
                        logging.debug(f'applying spikein ratio {sir} to column {col}')
                        nbcmdf[col] = nbcmdf[col] * sir
                        n_sicols += 1
                logging.info(f'brain={brain_id} spikein ratio(s) applied to {n_sicols} columns.')

            # fix column order AFTER potential spikein adjustments...
            nbcmdf = nbcmdf[ natsorted( list(nbcmdf.columns) ) ]
         
            logging.debug(f'nbcmdf.describe()=\n{nbcmdf.describe()}')
            
            norm_dict[brain_id] = nbcmdf
            spike_dict[brain_id] = sbcmdf
            real_dict[brain_id] = rbcmdf
            
            sh.add_value(f'/matrices/brain_{brain_id}','n_vbcs_real', len(rbcmdf) )
            sh.add_value(f'/matrices/brain_{brain_id}','n_vbcs_spike', len(sbcmdf) )
            
            # real matrix
            rbcmdf.to_csv(f'{outdir}/{brain_id}.rbcm.tsv', sep='\t')
            # spike-in matrix
            sbcmdf.to_csv(f'{outdir}/{brain_id}.sbcm.tsv', sep='\t')
            # real normalized by spike-ins.    
            nbcmdf.to_csv(f'{outdir}/{brain_id}.nbcm.tsv', sep='\t')          
            logging.info(f'done with brain={brain_id}')
        else:
            logging.warning(f'brain={brain_id} skipped.')    
        
    logging.info(f'got dict of {len(norm_dict)} normalized barcode matrices. returning.')
    return (real_dict, spike_dict, norm_dict)


def sync_columns(df1, df2, fillval=0.0):
    '''
    If two DFs don't have same columns, add empty columns to make the same. 
    Note columns may be missing in either DF, so number of columns is not enough to compare...
    '''
    x = set(df1.columns)
    y = set(df2.columns)
    w = x.difference(y) # columns to be added to df2
    z = y.difference(x) # columns to be added to df1    

    # only process if there is a problem...
    if len(w) > 0 or len(z) > 0:
        logging.debug('mismatched matrix columns, fixing..')
        for c in w:
            df2[c] = fillval
        for c in z:
            df1[c] = fillval
        
        scol = natsorted(list(df1.columns))
        df1 = df1[scol]
        scol = natsorted(list(df2.columns))
        df2 = df2[scol]    
        logging.debug(f'df1 col={df1.columns}')
        logging.debug(f'df2 col={df2.columns}')    
    else:
        logging.debug('matrix columns match.')
    
    return (df1, df2)


def get_si_ratio_map(sampdf, column = 'label'):
    '''
    create BCXXX -> si_ratio dict, with correct data types.
    handle possible other column as header
    
    '''
    # make label column as it isn't native. 
    sampdf['label'] = 'BC' + sampdf['rtprimer']
    sirdict = {}
    for i, ser in sampdf.iterrows():
        colval = ser[column]
        try:
            si_ratio = float(ser['si_ratio'])
        except ValueError:
            logging.warning(f"{ser['si_ratio']} not convertible to float. Setting to 1.0 ")
            si_ratio = 1.0
        sirdict[colval] = si_ratio
    return sirdict 


def normalize_weight_grouped(df, weightdf, columns=None, sampdf=None):
    '''
    Weight values in df by weightdf, by column groups
    
    df dataframe to be scaled, weightdf is spikeins
    columns is list of column groupings. weighting will done within the groups.     
    e.g.   [  ['BC3', 'BC9', 'BC1', 'BC2'], ['BC89'] ] where BC89 is the injection column.
      
    Assumes matrix index is sequence.
    Assumes matrices have same columns  
    If column numbers are mis-matched, will create empty column
    If columns is none, use/weight all columns
    Will only return columns somewhere in columns. 

    '''
    logging.debug(f'normalizing df=\n{df}\nby weightdf=\n{weightdf} columns={columns}')
    
    # sanity checks, fixes. 
    if len(df.columns) != len(weightdf.columns):
        logging.error(f'mismatched matrix columns df:{len(df.columns)} weightdf: {len(weightdf.columns)} !!')
    
    if columns is None:
        columns = [ list( df.columns ) ]    
    
    outdf = None
    
    for i, group in enumerate( columns ):
        logging.debug(f'handling group {i} of {len(group)} columns.')
        # we may get empty categories, depending on experiment. 
        if len(group) > 0:
            try: 
                gwdf = weightdf[group]
                gdf = df[group]
                sumlist = []
                for col in group:
                    sum = gwdf[col].sum()                       
                    sumlist.append(sum)
                sum_array = np.array(sumlist)
                maxidx = np.argmax(sum_array)
                maxval = sum_array[maxidx]  
                maxcol = gwdf.columns[maxidx]
                logging.debug(f'largest spike sum for {maxcol} sum()={maxval}')
                factor_array =  maxval / sum_array
                factor_array = np.nan_to_num(factor_array, posinf=0.0)
                logging.debug(f'factor array= {list(factor_array)}')
            
                max_list = []
                sum_list = []
                for col in gdf.columns:
                    max_list.append(gdf[col].max())
                    sum_list.append(gdf[col].sum())
                logging.debug(f'real max_list={max_list}')
                logging.debug(f'real sum_list={sum_list}')
                
                normdf = gdf.copy()
                for i, col in enumerate(normdf.columns):
                    logging.debug(f'handling column {col} idx {i} * factor={factor_array[i]}')
                    normdf[col] = (normdf[col] * factor_array[i] ) 
            
                max_list = []
                sum_list = []
                for col in normdf.columns:
                    max_list.append(normdf[col].max())
                    sum_list.append(normdf[col].sum())
                logging.debug(f'norm max_list={max_list}')
                logging.debug(f'norm sum_list={sum_list}')
                if outdf is None:
                    logging.debug(f'outdf is None, setting to new.')
                    outdf = normdf
                else:
                    logging.debug(f'outdf not None, before: len={len(outdf)}. cols={len(outdf.columns)}')
                    outdf = pd.concat([ outdf, normdf], axis=1)
                    logging.debug(f'after: len={len(outdf)}. cols={len(outdf.columns)}')
            except Exception as e:
                logging.warning(f'got exception processing a group.')    
        else:
            logging.warning('column group is empty, ignoring...')
       
    logging.info(f'normalizing done.len={len(outdf)}. cols={len(outdf.columns)} ')
    return outdf



#
#         REPORTS / QC
#
def make_read_report_xlsx(df,
                     outdir=None,
                     step='collapsed',                      
                     cp=None):
    '''
    creates excel sheet of unique VBC counts and UMI sums for 
    SSIs and types. 
    
    produced for readtable and collapsed. 
    collapsed numbers will be less than readtable since there 
    are fewer unique VBCs after collapse. Difference reflects components. 
     
    '''
    TYPE=['real','spike']
    
    if cp is None:
        cp = get_default_config()
    
    if outdir is None:
        outdir = os.path.abspath('./')
    project_id = cp.get('project','project_id')
    
    outfile = os.path.join(outdir, f'{project_id}.{step}.readreport.xlsx')    
    
    logging.info(f'creating unique VBC/read XLSX report: {outfile} ')
    
    vdf = df.groupby(by=['label','type'],observed=True).agg( {'vbc_read':'nunique','umi':'nunique'} )
    vdf.reset_index(inplace=True, drop=False)
    vdf.sort_values(by='label', inplace=True, key=lambda x: np.argsort(index_natsorted( vdf['label'])))
    vdf.reset_index(inplace=True, drop=True)

    sdf = df[df['type'] == 'spike']
    sdf = sdf.groupby(by=['label','type'],observed=True).agg( {'vbc_read':'nunique','umi':'nunique'} )
    sdf.reset_index(inplace=True, drop=False)
    sdf.sort_values(by='label', inplace=True, key=lambda x: np.argsort(index_natsorted( sdf['label'])))
    sdf.reset_index(inplace=True, drop=True)    

    rdf = df[df['type'] == 'real']
    rdf = rdf.groupby(by=['label','type'],observed=True).agg( {'vbc_read':'nunique','umi':'nunique'} )
    rdf.reset_index(inplace=True, drop=False)
    rdf.sort_values(by='label', inplace=True, key=lambda x: np.argsort(index_natsorted( rdf['label'])))
    rdf.reset_index(inplace=True, drop=True)

    rcdf = df.groupby(by=['label','type'],observed=True).agg( {'read_count':'sum'} )
    rcdf.reset_index(inplace=True,drop=False)
    rcdf.sort_values(by='label', inplace=True, key=lambda x: np.argsort(index_natsorted( rcdf['label'])))
    rcdf.reset_index(inplace=True, drop=True)    

    with pd.ExcelWriter(outfile) as writer:
        vdf.to_excel(writer, sheet_name='All Unique VBC')
        sdf.to_excel(writer, sheet_name='Spike Unique VBC')
        rdf.to_excel(writer, sheet_name='Real Unique VBC')
        rcdf.to_excel(writer, sheet_name='Read Count Sum')

    of = os.path.join(outdir, f'{project_id}.{step}.uvbc.tsv')
    vdf.to_csv(of, sep='\t')
    of = os.path.join(outdir, f'{project_id}.{step}.spike.tsv')
    sdf.to_csv(of, sep='\t')
    of = os.path.join(outdir, f'{project_id}.{step}.real.tsv')
    rdf.to_csv(of, sep='\t')
    of = os.path.join(outdir, f'{project_id}.{step}.rcount.tsv')
    rcdf.to_csv(of, sep='\t')
    
    logging.info(f'Wrote XLSX report: {outfile} ')


def make_rtag_counts(df,
                     outdir,
                     cp=None):
    '''
    Save summary of rtag/rrtag value counts. 
    '''
    if cp is None:
        cp = get_default_config()
    
    try:
        rtag = df['rtag']
        rvc = rtag.value_counts()
        rdf = pd.DataFrame()
        rdf['rseq_count'] = rvc
        rdf.reset_index(inplace=True)
        tot_count = rdf['rseq_count'].sum()
        rdf['rcum_prop'] = rdf['rseq_count'].cumsum() / tot_count
    
        of = os.path.join(outdir, 'rtag_counts.tsv')
        write_mapseq_df(rdf, of, outformats=['tsv'])
        logging.debug(f'wrote rtag info to {of}')
        
        rrtag = df['rrtag']
        rrvc = rrtag.value_counts()
        rrdf = pd.DataFrame()
        rrdf['rrseq_count'] = rrvc
        rrdf.reset_index(inplace=True)    
        tot_count = rrdf['rrseq_count'].sum()
        rrdf['rrcum_prop'] = rrdf['rrseq_count'].cumsum() / tot_count
        of = os.path.join(outdir, 'rrtag_counts.tsv')
        write_mapseq_df(rrdf, of, outformats=['tsv'])
        logging.debug(f'wrote rrtag info to {of}')
    except Exception as ex:
        logging.warning(f'DF does not have rtag column?')
        logging.warning(traceback.format_exc(None))

def make_spikein_qctable(df,
                         outfile = None, 
                         cp=None
                         ):
    '''
    gather QC info about readtable data.
         
    ''' 
    logging.debug(f'QC table.')
    if cp is None:
        cp = get_default_config()

    project_id = cp.get('project','project_id')
        
    if outfile is None:
        outdir = os.path.abspath('./')
        outfile = os.path.join( outdir, f'{project_id}.spike.qctable.xlsx')
        logging.debug(f'outdir not provided, set to {outdir}')
    else:
        dirpath, base, ext = split_path( os.path.abspath(outfile))
        outdir = dirpath
        logging.debug(f'outdir = {outdir} ')
            
    sh = get_default_stats()
    xlout = outfile
    
    logging.debug(f'writing to {xlout}')
    with pd.ExcelWriter( xlout) as writer:
        srdf = df[df['type'] == 'spike']
        srdf['vbcumi'] = srdf['vbc_read'] + srdf['umi']
        sinfo = srdf.groupby(['label'], observed=True).agg( {'vbc_read':'nunique',
                                                             'umi':'nunique', 
                                                             'vbcumi':'nunique'}).reset_index()
        sinfo = sinfo.rename(columns={ 'vbc_read':  'uniq_vbc_read',
                                        'umi'    :  'uniq_umi',
                                        'vbcumi' :  'total_umi' })
        # Fix label order
        sinfo.sort_values(by='label', inplace=True, key=lambda x: np.argsort( index_natsorted( sinfo['label'])))
        sinfo.reset_index(inplace=True, drop=True)
        sinfo.to_excel(writer)
    logging.info(f'wrote out to {xlout}')        


def make_controls_umireport_xlsx(df, 
                                 outdir=None, 
                                 cp=None, 
                                 cols=['site','type'], 
                                 vals = ['vbc_read','label','umi_count','read_count'], 
                                 sort_by='umi_count',
                                 step='vbctable',
                            ):
    '''
    write subsetted data from vbctable dataframe for QC checks. 

    vbc_read                  label  type  umi_count  read_count brain region  site
0 AAAAAAACAGCTAAAGAATCCTTGTTCACC  BC207  lone          3          25               target-water-control
1 AAAAAAACATTCACGCGTATGGCCTGAGGG  BC207  lone          1          20               target-water-control
2 AAAAAAACCGGCCTTGTACTTGGTTCTCTT  BC203  real         13         509     6                       target
  
    '''
    if cp is None:
        cp = get_default_config()
        
    if outdir is None:
        outdir = os.path.abspath('./')
        logging.debug(f'outdir not provided, set to {outdir}')
    else:
        outdir = os.path.abspath(outdir)
        logging.debug(f'outdir = {outdir} ')
    
    project_id = cp.get('project','project_id')
    sh = get_default_stats()
       
    # for all sites, for all types, create
    # <site>.<type>.tsv  ->  label  umi_count read_count 
    #
    labeldict = {}
    for c in cols:
        clist = list(df[c].dropna().unique())
        clist = [ x for x in clist if len(x) > 0 ]        
        clist.sort()
        labeldict[c] = clist
    
    logging.debug(f'labeldict={labeldict}')
    alist = []
    for c in cols:
        alist.append(labeldict[c])
    logging.debug(f'arglist = {alist}')
    combinations = list(itertools.product(*alist))
    
    try: 
        xlout = os.path.join(outdir, f'{project_id}.{step}.controls.umireport.xlsx')
        logging.debug(f'writing to {xlout}')
        with pd.ExcelWriter( xlout) as writer:
            for tup in combinations:
                logging.debug(f'handling columnset {tup}')
                tsvname = []
                ndf = df.copy()
                for i, col_name in enumerate(cols):
                    col_val = tup[i]
                    tsvname.append(col_val)
                    logging.debug(f'filtering {col_name} by value={col_val}')
                    ndf = ndf[ ndf[col_name] == col_val ]
                fname = '.'.join(tsvname)
                outfile = os.path.join(outdir, f'{fname}.tsv')
                if len(ndf) > 1:
                    logging.debug(f'writing to {outfile}')
                    ndf = ndf[vals]
                    ndf.sort_values(by=[sort_by], ascending=False, inplace=True)
                    ndf.reset_index(inplace=True, drop=True)
                    ndf.to_csv(outfile, sep='\t')
                    sh.add_value('/{step}',f'n_{fname}_vbcs', len(ndf) )
                    # If it is a real control, add to control report. 
                    for s in CONTROL_SITES:
                        if s in fname:
                            logging.debug(f"writing control '{fname}' len={len(ndf)} to {xlout}")
                            ndf.to_excel(writer, sheet_name=fname)
                        else:
                            pass
                else:
                    logging.info(f'no entries for {fname}')
             
            logging.debug(f'done with all tuples in {combinations}')
    except Exception as ex:
        logging.warning(f'exception = {ex}')    


def make_vbctable_parameter_report_xlsx(df, 
                                        outdir=None, 
                                        cp=None, 
                                        params=None
                                        ):
    '''
        Create single XLSX with table of number of VBCs per SSI at the various 
        thresholds for injection, target UMI counts. 
            For reference:
            VBCs per brain at different inj UMI thresholds:
            
                   10      5
            B13   1124    1244
            B14   2232    2445
            B16   2029    2129
            B18   608     708
            B19   1295    1450
            B21   3175    3629
            B22   3986    4673
            B23   1322    1547
            B24   3570    3954
            
    PARAMS are  <inj_min_umi> , <target_min_umi>         
    '''
    DEFAULT_PARAMS = [ (5,  3), 
                       (10, 3),
                       (10, 5),
                       (20, 5),
                       (30, 5),
                       (30, 10)]
      
    if cp is None:
        cp = get_default_config()
    if params is None:
        params = eval(  cp.get( 'vbctable', 'test_params', fallback='[ (5,3) ,(10,3), (10,5), (20,5), (30,5),(30,10)]')  ) 
    if outdir is None:
        outdir = os.path.abspath('./')
    else:
        outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    
    project_id = cp.get('project','project_id')        
    outfile = os.path.join(outdir, f'{project_id}.vbc_parameters.xlsx')    
    logging.info(f'testing range of UMI parameters (inj,target) to {outfile} ')
    logging.debug(f'params={params} ')

    outdf = None
    
    # Explicitly set include_injection for prosective matrix measurement. 
    # 
    testcp = copy.deepcopy(cp)
    testcp.set('vbcfilter', 'include_injection', 'True')
    
    for i, (inj_min_umi, target_min_umi) in enumerate(params):
        logging.debug(f'testing parameters: inj_min_umi = {inj_min_umi} target_min_umi = {target_min_umi} ')
        colname = f'inj:{inj_min_umi}, tar:{target_min_umi}'
        (fdf, controls) = process_filter_vbctable(df, 
                                      inj_min_umi=inj_min_umi, 
                                      target_min_umi = target_min_umi, 
                                      target_min_umi_absolute=1, 
                                      outdir = outdir, 
                                      cp=testcp)
        logging.debug(f'got filtered vbctable:\n{fdf}')
        fdf = fdf[fdf['type'] == 'real']
        # deal with duplicated controls originally without a brain ID. 
        xdf = fdf.groupby('label').agg({'vbc_read':'nunique'})
        xdf.reset_index(inplace=True, drop=False)
        xdf.sort_values(by='label', inplace=True, key=lambda x: np.argsort(index_natsorted( xdf['label'])))
        xdf.reset_index(inplace=True, drop=True)
        xdf.rename( {'vbc_read': colname}, inplace=True, axis=1)
        xdf.index = xdf['label']
        xdf.drop('label', inplace=True, axis=1)
        
        if outdf is None:
            outdf = xdf
        else:
            outdf = pd.concat([ outdf, xdf], axis=1)   

    with pd.ExcelWriter(outfile) as writer:
        outdf.to_excel(writer, sheet_name='Unique VBCs by Parameters')
        of = os.path.join(outdir, f'{project_id}.vbc_parameters.tsv')
        outdf.to_csv(of, sep='\t')
    logging.info(f'Wrote XLSX report: {outfile} ')

    return outdf



############################################################
#
#        Assessment/QC/Validation/Simulation
#
############################################################
def qc_make_readmatrix( df, sampdf=None, outdir='./', cp=None):
    '''
    make matrix of SSI vs. FASTQ  
    
    '''
    logging.debug(f'inbound df len={len(df)} columns={df.columns}')
    sh = get_default_stats()

    if cp is None:
        cp = get_default_config()
    bcfile = os.path.expanduser( cp.get('barcodes','ssifile') ) 
    project_id = cp.get('project','project_id')

    # Map label, rtprimer to SSIs    
    logging.debug(f'getting rt labels...')
    labels = None
    if sampdf is not None:
        labels = get_rtlist(sampdf)
    n_reads = len(df)
    sh.add_value('/qcreads/sources', 'n_reads', str(n_reads ) )
    logging.debug(f'rtlabels={labels}')
    bcdict = get_barcode_dict(bcfile, labels)
    n_ssis =  len(bcdict)
    sh.add_value('/qcreads/sources', 'n_ssis', str(n_ssis) )
    rtdict = get_rtprimer_dict(bcfile, labels)
    rtbcdict = get_rtbc_dict(bcfile, labels)
    logging.debug(f'got {len(bcdict)} barcodes with labels and primer number.')    

    logging.info('filling in rtprimer number by SSI sequence...')    
    df['rtprimer'] = df['ssi'].map(rtdict)
    df['rtprimer'] = df['rtprimer'].astype('category')
    
    logging.info('filling in label by rtprimer...')
    df['label'] = df['rtprimer'].map(rtbcdict)
    df['label'] = df['label'].astype('category')
    
    gdf = df.dropna()
    gdf.reset_index(inplace=True, drop=True)
    n_readswssi =  len(gdf)
    sh.add_value('/qcreads/sources', 'n_readswssi', str(n_readswssi) )
    good_pct = ( n_readswssi / n_reads ) * 100
    good_spct = f'{good_pct:.2f}'
    sh.add_value('/qcreads/sources', 'n_reads_valid_ssi', str(good_spct) )
    
    gadf = gdf.groupby(by=['label','source'], observed=True).agg( {'source':'count'})
    gadf.columns = ['count']
    gadf.reset_index(inplace=True,drop=False)
    
    sldf = gadf.pivot(index='source',columns='label',values='count')
    sldf.fillna(0.0, inplace=True)
    scol = natsorted(list(sldf.columns))
    sldf = sldf[scol]
    nsindex = natsorted(sldf.index)
    sldf = sldf.loc[nsindex]
    
    # Calc percentage in non-dominant, by SSI
    byssi = []
    for col in list(sldf.columns):
        bcc = sldf[col].sort_values(ascending=False)
        n_dom = bcc.iloc[0]
        f_dom = bcc.index[0]
        n_total = bcc.sum()
        dom_pct = ( n_dom / n_total ) * 100
        sdom_pct = f'{dom_pct:.2f}'
        sh.add_value(f'/qcreads/sources/ssi/{col}', 'dom_pct', str(sdom_pct) )
        sh.add_value(f'/qcreads/sources/ssi/{col}', 'dom_file', str(f_dom) )
        byssi.append( { 'ssi'     : col , 
                          'dom_pct': dom_pct , 
                         'dom_file':  str(f_dom)} )
    byssidf = pd.DataFrame(byssi)
    
    # Calc percentage in non-dominant, by source file.
    byfile = []
    for ridx in list(sldf.index):
        scc = sldf.loc[ridx].sort_values(ascending=False)
        n_dom = scc.iloc[0]
        s_dom = scc.index[0]
        n_total = scc.sum()
        dom_pct = ( n_dom / n_total ) * 100
        sdom_pct = f'{dom_pct:.2f}'
        sh.add_value(f'/qcreads/sources/ssi/{ridx}', 'dom_pct', str(sdom_pct) )
        sh.add_value(f'/qcreads/sources/ssi/{ridx}', 'dom_ssi', str(s_dom) )
        byfile.append( {'file': str(ridx), 
                          'dom_pct': dom_pct , 
                          'dom_ssi':  str(s_dom)} )
    byfiledf = pd.DataFrame(byfile)
            
    outfile = os.path.join(outdir, f'{project_id}.read_sources.xlsx')
    with pd.ExcelWriter(outfile) as writer:
        sldf.to_excel(writer, sheet_name='Source SSI matrix')
        byssidf.to_excel(writer,sheet_name='Dominant by SSI' )
        byfiledf.to_excel(writer,sheet_name='Dominant by file' )
    return sldf


class RandomSet():
    '''
    Allows random sampling without cast set to list. 
    No removals!! Add only. 
    '''
    
    def __init__(self):
        self.data_map = {}
        self.data = []

    def add(self, element):
        if element in self.data_map:
            return False

        self.data_map[element] = len(self.data)
        self.data.append(element)

    def update(self, element_list):
        for elem in element_list:
            self.add(elem)
        
    def discard(self, element):
        if not element in self.data_map:
            return False
        
        last_element = self.data[-1]
        rem_idx = self.data_map[element]
        self.data_map[last_element] = rem_idx
        self.data[rem_idx] = last_element
        self.data[-1] = element
        
        self.data.pop()
        self.data_map.pop(element)
        return True
    

    def getRandom(self):
        return random.choice(self.data)

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        return str(self.data)


def generate_mutation_cycle( sequence_list,
                             n_copies=5,
                             max_mismatch=3,
                             alphabet='CAGT'
                             ):
    '''
    Simulation of virus replication cycle. 
    
    '''
    logging.debug(f'mutation_cycle: input: {len(sequence_list)} sequences. max_mismatch={max_mismatch}')
    output_list = sequence_list.copy()
    sequences_handled = 0
    LOG_INTERVAL =  50
    for seq in sequence_list:
        for i in range(0, n_copies):
            mseq = mutate_sequence(seq, alphabet, max_mismatch)
            output_list.append(mseq)
        sequences_handled += 1
        if i % LOG_INTERVAL == 0:
            logging.debug(f'handled {i} sequences...')
    logging.debug(f'mutation_cycle: output: {len(output_list)} sequences.')
    return output_list

def generate_mutated_list( sequence_list, 
                           n_copies=5, 
                           n_cycles=5, 
                           max_mismatch=3, 
                           alphabet='CAGT'):
    '''
    Simple list version of mutation cycle. 
    Each cycle takes output of entire previous cycle. 
    Uses list and allows duplicates. 
    Mimics dominance of parent sequence.  
    
    '''
    logging.debug(f'args sequence_list n={len(sequence_list)} n_cycles={n_cycles} max_mismatch={max_mismatch}')
    output_list = sequence_list.copy()
    for i in range(0, n_cycles):
        logging.debug(f'cycle {i} ')
        output_list = generate_mutation_cycle( output_list, n_copies, max_mismatch, alphabet)
    logging.info(f'generated mutation list. len={len(output_list)} n_cycles={n_cycles} max_mismatch={max_mismatch} ')
    return output_list


def generate_mutated_df(   sequence_df, 
                           n_copies=5, 
                           n_cycles=5, 
                           max_mismatch=3, 
                           alphabet='CAGT'):
    '''
    Each cycle takes output of entire previous cycle. 
    Uses list and allows duplicates. 
    Mimics dominance of parent sequence.
    Labels every mutated sequence with the parent_idx of its parent.   
    
    '''
    logging.debug(f'args sequence_df n={len(sequence_df)} n_cycles={n_cycles} max_mismatch={max_mismatch}')
    for i in range(0, n_cycles):
        logging.debug(f'cycle {i} ')
        sequence_df = mutate_sequence_df(  sequence_df,
                                           n_copies, 
                                           max_mismatch, 
                                           alphabet )        
    logging.info(f'generated mutation list. len={len(sequence_df)} n_cycles={n_cycles} max_mismatch={max_mismatch} ')
    return sequence_df


def mutate_sequence_df( sequence_df,
                        n_copies, 
                        max_mismatch, 
                        alphabet ):
    '''
    Performs a single cycle of mutation on a sequence DF.
    Each sequence is mutated <n_copies> times and added to output. 
    Keep copy of all other columns in output, and include all contents of sequence_df.
    
    '''
    logging.debug(f'args sequence_df n={len(sequence_df)} n_copies={n_copies} max_mismatch={max_mismatch}')
    other_columns = list( sequence_df.columns )
    other_columns.remove('sequence')
    logging.debug(f'will retain columns {other_columns}')
    concat_list = []
    concat_list.append(sequence_df)
    for i in range(0, len(sequence_df)):
        mlist = []
        row =  sequence_df.iloc[i]
        for j in range(0, n_copies):
            mseq = mutate_sequence(row['sequence'])
            rlist = list( sequence_df.iloc[i])
            mlist.append(mseq)
        copy_df = pd.DataFrame(columns = sequence_df.columns)
        mser = pd.Series(mlist)
        copy_df['sequence'] = mser
        for c in other_columns:
            copy_df[c] = row[c]
        #logging.debug(f'copy_df=\n{copy_df}')
        concat_list.append(copy_df)
    #logging.debug(f'made concat list = {concat_list}')
    output_df = pd.concat(concat_list, copy=False, ignore_index=True )
    logging.debug(f'completed cycle. new n={len(output_df)}')
    return output_df


def mutate_sequence(sequence, alphabet='CAGT', max_mismatch=3):
    '''
    Take sequence, perform random mutation from alphabet to random position(s) 
    between 1 and max_mismatch times, with equal probability of each.  
    Return mutated sequence. 
    
        
    '''            
    alph_ba = bytearray(alphabet, encoding='UTF-8')
    n_edits = random.randint(1, max_mismatch)
    n_chars = len(sequence)
    sba = bytearray(sequence, encoding='UTF-8')
    for i in range(0, n_edits):
        m_idx = random.randint(0, n_chars - 1)
        sba[m_idx] = random.choice(alph_ba)
        #logging.debug(f'mutating index {m_idx} to {sba[m_idx]} ')
    return sba.decode('UTF-8')


def generate_simulated_vbcreads(n_parents=100):
    '''
    Simulate viral replication with mutation. 
    '''
    (parent_list, parent_df) = generate_random(n_sequences=n_parents)
    logging.debug(f'got parent list, df len={len(parent_df)}')

    mutated_list = generate_mutated_list(parent_list)    
    logging.debug(f'got flattened mutated list len={len(mutated_list)}')
    
    df = pd.DataFrame()
    sser = pd.Series(mutated_list)
    df['vbc_read'] = sser
    return df

def generate_simulated_vbcreads_df(n_parents=100):
    '''
    Simulate viral replication with mutation.
    Include parent_idx column for all sequences, including mutated. 
    '''
    (parent_list, parent_df) = generate_random(n_sequences=n_parents)
    parent_df['parent_idx'] = parent_df.index
    logging.debug(f'got parent df, df len={len(parent_df)}')

    mutated_df = generate_mutated_df(parent_df)    
    logging.debug(f'got flattened mutated list len={len(parent_df)}')
  

      
def generate_random(n_sequences=1000, n_bases=30, alphabet='CAGT'):
    '''
    Generate random vbc_reads/sequence dataframe.  
    Aggregate to remove (unlikely) duplicates. 
    Filter by homopolymers

    vbc_read   count

    '''
    logging.debug(f'args n_bases={n_bases} n_sequences={n_sequences}')
    ALPH=['C','A','G','T']
    LOG_INTERVAL=1000000
    
    sequence_list = []
    for i in range(0, n_sequences):
        seq = ''.join(random.choices(ALPH, k=n_bases))
        sequence_list.append(seq)
        if i % LOG_INTERVAL == 0:
            logging.debug(f'created {i} sequences...')
    ss = pd.Series(sequence_list, dtype='string[pyarrow]')
    df = pd.DataFrame(ss, columns=['sequence'])

    vcs = df.value_counts( ['sequence'] )
    vcs = vcs[ vcs == 1 ]
    df = vcs.reset_index()
    df.columns = ['sequence','read_count']
    df = filter_homopolymers(df, column='sequence')
    df['parent_idx'] = df.index
    seqlist = list(df['sequence'])
    return ( seqlist, df )




