import glob
import gzip
import json
import logging
import math
import os
import subprocess
import sys
import traceback

import datetime as dt

from configparser import ConfigParser
from collections import defaultdict

import pandas as pd
import numpy as np
from natsort import natsorted
import matplotlib.pyplot as plt
import seaborn as sns
import scipy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from dask.dataframe import from_pandas
import dask
import dask.dataframe as dd

from mapseq.utils import *
from mapseq.bowtie import *
from mapseq.barcode import *
from mapseq.stats import *


def get_default_config():
    dc = os.path.expanduser('~/git/mapseq-processing/etc/mapseq.conf')
    cp = ConfigParser()
    cp.read(dc)
    return cp


def get_rtlist(sampledf):
    rtlist = list(sampledf['rtprimer'].dropna())
    nrtlist = []
    for x in rtlist:
        try:
            y = int(float(x))
            nrtlist.append(y)
        except ValueError:
            logging.debug(f'ignoring bad int() input.')
    
    #rtlist = [int(x) for x in rtlist]
    
    nrtlist = [f'BC{x}' for x in nrtlist]
    return nrtlist


def package_pairfiles(infiles):
    '''
    pack up input list of elements into list of paired tuples. 
    ['a','b','c','d'] -> [('a','b'),('c','d')]

    '''
    if len(infiles) %2 != 0:
        logging.error(f'number of elements must be multiple of two!')
        
    infilelist = []
    a = None
    b = None
    for i,v in enumerate(infiles):
        if i % 2 == 0:
            a = v
        else:
            b = v
            t = (a,b)
            logging.info(f'input pair of readfiles: r1={a} r2={b}')
            infilelist.append(t)
    return infilelist


def guess_site(infile, sampdf):
    '''
    will look at filename and try to guess rt primer number, then 
    look for siteinfo in sampledf
    
    NOTE: assumes BC<rtprimer>.fasta or SSI<rtprimer>.fasta and <rtprimer> identifiers
    consist of digits. 
    
    '''
    logging.info(f'guessing site/brain/region for FASTA file {infile}')
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    head = base.split('.')[0]
    rtprimer_num = ''.join(i for i in head if i.isdigit())
    rtprimer_num = str(rtprimer_num)
    logging.debug(f'base={base} head={head} guessing rtprimer={rtprimer_num} sampdf=\n{sampdf}')
    logging.debug(f'sampdf.dtypes = {sampdf.dtypes}')
    df = sampdf[sampdf['rtprimer'] == rtprimer_num]
    df.reset_index(inplace=True, drop=True)
    logging.debug(f'rtprimer specific df = \n{df}')
    site = None
    brain = None
    region = None
    if len(df)> 0:
        try:
            site = str(df['siteinfo'][0])
            logging.debug(f'got site={site} successfully')
        except Exception as e:
            logging.warning(f'unable to get siteinfo for {infile}')
            site = 'target'  # default to target. 
        
        try:        
            brain = str(df['brain'][0])
            logging.debug(f'got brain={brain} successfully')
        except Exception as e:
            logging.warning(f'unable to get brain info for {infile}')
            brain = '0'
            
        try:        
            region = str(df['region'][0])
            logging.debug(f'got region={region} successfully')
        except Exception as e:
            logging.warning(f'unable to get region info for {infile}')
            region = str(rtprimer_num) # default to SSI label      
            
        logging.debug(f'got site={site} brain={brain} region={region} for rtprimer={rtprimer_num}') 

    logging.debug(f'got site={site} for rtprimer guessed from {infile}')
    return (rtprimer_num, site, brain, region )
    



def read_threshold_all(cp, alldf):
    '''
    we did not threshold BC by BC when processing fastqs, so we can threshold by proxy
    using a ratio of UMIs to read_count behind them. 
    '''
    target_ratio = float( cp.get('ssifasta','target_read_ratio'))
    injection_ratio = float( cp.get('ssifasta','injection_read_ratio'))
    
    alldf['read_ratio'] = alldf['read_count'] / alldf['umi_count']
        
    injdf = alldf[ alldf['site'].str.startswith('injection') ]
    tardf = alldf[ alldf['site'].str.startswith('target') ]
    injdf = injdf[ injdf['read_ratio'] > injection_ratio ]
    tardf = tardf[ tardf['read_ratio'] > target_ratio ]
    df_merged = pd.concat([injdf, tardf], ignore_index = True, sort=False)
    logging.info(f'target_ratio = {target_ratio} injection_ratio = {injection_ratio} before len={len(alldf)} after len={len(df_merged)}')
    return df_merged
       
    


def edges_from_btdf(btdf):
    readlist = btdf.name_read.values.tolist()
    alignlist = btdf.name_align.values.tolist()  
    edgelist = [ list(t) for t in zip(readlist, alignlist)]
    return edgelist

def get_components(edgelist, integers=True):
    '''
    returns strongly connected components from list of edges via Tarjan's algorithm. 
    assumes labels are integers (for later use as indices in dataframes. 
    '''
    complist = []
    logging.debug(f'getting connected components from edgelist len={len(edgelist)}')
    if len(edgelist) < 100:
        logging.debug(f'{edgelist}')
    for g in tarjan(from_edges(edgelist)):
        #logging.debug(f'g={g}')
        complist.append(g)
    logging.debug(f'{len(complist)} components.')
    if len(edgelist) < 100:
        logging.debug(f'{complist}')
    if integers:
        outlist = []
        for g in complist:
            outlist.append( [int(x) for x in g])
        complist = outlist    
    return complist

#
#  Tarjan's algorithm, same as Matlab graphconncomp()   
#  https://rosettacode.org/wiki/Tarjan#Python:_As_function
#
def from_edges(edges):        
    class Node:
        def __init__(self):
            # root is one of:
            #   None: not yet visited
            #   -1: already processed
            #   non-negative integer: what Wikipedia pseudo code calls 'lowlink'
            self.root = None
            self.succ = []

    nodes = defaultdict(Node)
    for v,w in edges:
        nodes[v].succ.append(nodes[w])

    for i,v in nodes.items(): # name the nodes for final output
        v.id = i

    return nodes.values()
    
    
def tarjan(V):
    '''
    May get recursion limit errors if input is large. 
    https://stackoverflow.com/questions/5061582/setting-stacksize-in-a-python-script/16248113#16248113
    
    import resource, sys
    resource.setrlimit(resource.RLIMIT_STACK, (2**29,-1))
    sys.setrecursionlimit(10**6)
    
    Same algorithm as used in MATLAB pipeline. 
    https://www.mathworks.com/help/matlab/ref/graph.conncomp.html 
    
    
    '''
    def strongconnect(v, S): 
        v.root = pos = len(S)
        S.append(v)
        for w in v.succ:
            if w.root is None:  # not yet visited
                yield from strongconnect(w, S)

            if w.root >= 0:  # still on stack
                v.root = min(v.root, w.root)
        if v.root == pos:  # v is the root, return everything above
            res, S[pos:] = S[pos:], []
            for w in res:
                w.root = -1
            yield [r.id for r in res]
    
    for v in V:
        if v.root is None:
            yield from strongconnect(v, [])


def make_read_counts_df(config, seqdf, label=None):
    '''
    input dataframe with 'sequence' column
    make counts column for identical sequences.  
    optionally assign a label to set in new column
    
    '''
    logging.debug(f'seqdf=\n{seqdf}')
    ser = seqdf['sequence'].value_counts()
    df = pd.DataFrame(columns=['sequence','read_count'])
    df['sequence'] = ser.index
    df['read_count'] = ser.values
    logging.debug(f'counts df = \n{df}')
    if label is not None:
        df['label'] = label
    return df

def make_umi_counts_df(config, seqdf, label=None):
    '''
    input dataframe with 'sequence' column
    make counts column for identical sequences.  
    optionally assign a label to set in new column
    combine values from collapsed read_counts 
    sequence   read_count  label     
    AAA        23           A
    AAA        5            A
    BBB        7            A
    BBB        6            A
    CCC        5            A
    
    ->
    AAA    read_count  umi_count label
    BBB    13           2          A
    CCC    28           2          A
    CCC    5            1          A
    
    https://stackoverflow.com/questions/73874908/value-counts-then-sum-of-a-different-column
    https://stackoverflow.com/questions/46431243/how-to-get-groupby-sum-of-multiple-columns
     
    '''
    logging.debug(f'seqdf=\n{seqdf}')
    # keep original label
    label = seqdf['label'].unique()[0]
    vdf = seqdf.drop('label',axis=1).groupby('sequence').agg({'read_count':'sum','umi_count':'sum'}).reset_index()
    #vdf = seqdf.groupby('sequence')['read_counts'].agg(count='count',sum='sum').reset_index()
    #vdf.rename(columns={'count': 'umi_counts','sum': 'read_counts'}, inplace=True)
    vdf.sort_values(by=['umi_count'], ascending=False, inplace=True)
    vdf.reset_index(inplace=True, drop=True)
    vdf['label'] = label
    logging.debug(f'counts df = \n{vdf}')
    return vdf


def make_fasta_df(config, infile, ignore_n=True):
    '''
    input fasta 
    ignore 'N' sequences.
    '''   
    slist = []
    rcs = SeqIO.parse(infile, "fasta")
    handled = 0
    for sr in rcs:
        s = sr.seq
        if ('N' in sr) and ignore_n :
            pass
        else:
            slist.append(str(s))
        handled += 1    
    logging.debug(f"kept {len(slist)} sequences out of {handled}")    
    df = pd.DataFrame(slist, columns=['sequence'] )
    return df

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
    #df['umi_count'] = df['umi_count'].astype(int)
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



def calc_min_threshold(config, df, site='target-negative'):
    '''
    retrieve maximum umi_count value for provided site type. 
    returns max max value...

    '''
    logging.debug(f'getting UMI threshold type={site}')
    countlist = []
    min_threshold = 0
    df['umi_count'] = df['umi_count'].astype(int)
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



def get_read_count_threshold(config, cdf, site=None):
    '''
    site = ['target-X','injection-X'] where X is '', lone, control, negative  
        Will use relevant threshold. If None, will use default threshold
    
    target_threshold=100
    target_ctrl_threshold=1000
    inj_threshold=2
    inj_ctrl_threshold=2
    
    '''
    logging.debug(f'getting info for site={site}')
    count_threshold=1
    if site is None:
        count_threshold = int(config.get('ssifasta', 'default_threshold'))
    else:
        count_threshold = int(config.get('ssifasta', f'{site}_threshold'))
    logging.debug(f'count threshold for {site} = {count_threshold}')
    return count_threshold


def calc_freq_threshold(df, fraction=0.9, column = 'read_count'):
    '''
    sorts column of input column
    calculates index of point at which <fraction> of data points are less 
    returns column value at that point + 1 
    '''
    ser = df[column].copy()
    ser.sort_values(ascending = False, inplace=True)
    ser.reset_index(drop=True, inplace=True)
    idx = int(len(ser) * fraction)
    t = int( ser.iloc[idx] + 1 )    
    return t
        

def threshold_read_counts(config, df, threshold=1):
    '''
    
    '''
    logging.debug(f'threshold counts threshold={threshold}')
    threshold= int(threshold)   
    df = df[df['read_count'] >= threshold].copy()
    return df

def load_sample_info(config, file_name, sheet_name='Sample information'):
    #
    # Parses Excel spreadsheet to get orderly sample metadata, saves as sampleinfo.tsv.     
    # OR Reads in sampleinfo.tsv
    # Assumes various properties of spreadsheet that need to stay static. 
    #
    #   ['Tube # by user', 'Our Tube #', 'Sample names provided by user',
    #   'Site information', 'RT primers for MAPseq', 'Brain ', 'Column#']
    #
    # If brain is not given, or is empty, all are set to 'brain1'. 
    # If region is not given, or is empty, all are set to <rtprimer>
    # 
    
    # Mappings for excel columns. 
    sheet_to_sample = {
            'Tube # by user'                  : 'usertube', 
            'Our Tube #'                      : 'ourtube', 
            'Sample names provided by user'   : 'samplename', 
            'Site information'                : 'siteinfo',
            'RT primers for MAPseq'           : 'rtprimer',
            'Brain'                           : 'brain',
            'Region'                          : 'region',
            'Matrix Column'                   : 'matrixcolumn',
        }
    
    sample_columns = ['usertube', 'ourtube', 'samplename', 'siteinfo', 'rtprimer', 'brain', 'region', 'matrixcolumn'] 
    int_sample_col = ['usertube', 'ourtube', 'rtprimer', 'region', 'matrixcolumn']     # brain is sometimes not a number. 
    str_sample_col = ['usertube', 'ourtube', 'samplename', 'siteinfo', 'rtprimer', 'brain' ,'region']

    if file_name.endswith('.xlsx'):
        edf = pd.read_excel(file_name, sheet_name=sheet_name, header=1)        
        sdf = pd.DataFrame()
        
        for ecol in edf.columns:
            ecol_stp = ecol.strip()    
            try:
                # map using stripped column name, retrieve using actual excel column name
                # which may have trailing spaces...
                scol = sheet_to_sample[ecol_stp]
                logging.debug(f'found mapping {ecol} -> {scol}')
                cser = edf[ecol]
                logging.debug(f'column for {scol}:\n{cser}')
                sdf[scol] = cser
            
            except KeyError:
                logging.debug(f'no mapping for {ecol} continuing...')
            
            except Exception as ex:
                logging.error(f'error while handling {ecol} ')
                logging.info(traceback.format_exc(None))
                
        logging.debug(f"rtprimer = {sdf['rtprimer']} dtype={sdf['rtprimer'].dtype}")
        
        # Only keep rows with rtprimer info. 
        sdf = sdf[sdf['rtprimer'].isna() == False]

        # Fix brain column
        sdf.loc[ sdf.brain.isna(), 'brain'] = 0
        try:
            sdf.brain = sdf.brain.astype('int')
        except ValueError:
            pass
        sdf.brain = sdf.brain.astype('string')
        sdf.loc[ sdf.brain == '0', 'brain'] = ''

        # fix rtprimer column, if it was read as float string (e.g '127.0' )
        #sdf['rtprimer'] =  sdf['rtprimer'].astype(float).astype(int).astype(str)
        
        for scol in sample_columns:
            try:
                ser = sdf[scol]
            except KeyError as ke:
                logging.info(f'no column {scol}, required. Creating...')
                if scol == 'samplename':
                    sdf[scol] = sdf['ourtube']
                elif scol == 'region':
                    sdf[scol] = sdf['rtprimer']
        
        # fix empty rows. 
        sdf.brain = sdf.brain.astype(str)
        #sdf.brain[sdf.brain == 'brain-nan'] = ''
        sdf.loc[sdf.brain == 'brain-nan','brain'] = ''
                    
        sdf.replace(r'^s*$', float('NaN'), regex = True, inplace=True)
        sdf.dropna(how='all', axis=0, inplace=True)        
        sdf = fix_columns_int(sdf, columns=int_sample_col)
        sdf = fix_columns_str(sdf, columns=str_sample_col)

    elif file_name.endswith('.tsv'):
        sdf = pd.read_csv(file_name, sep='\t', index_col=0, keep_default_na=False, dtype =str, comment="#")
        #df.fillna(value='', inplace=True)
        sdf = sdf.astype('str', copy=False)    
        sdf = fix_columns_int(sdf, columns=int_sample_col)
    else:
        logging.error(f'file {file_name} neither .xlsx or .tsv')
        sdf = None
        
    logging.debug(f'created reduced sample info df:\n{sdf}')
    return sdf

def load_seqtsv_pd(infile):
    '''
    handle reads output of process_fastq_pairs.
    52nt string sequence 
    '''
    logging.debug(f'loading reads TSV from {infile}')
    STR_COLUMNS = [ 'vbc_read','spikeseq', 'ssi','umi','libtag']
    INT_COLUMNS = ['read_count']
    df = pd.read_csv(infile, sep='\t', index_col=0)
    logging.debug(f'dtypes={df.dtypes}')
    for col in STR_COLUMNS:
        logging.debug(f'converting {col} to string[pyarrow]')
        try:
            df[col] = df[col].astype('string[pyarrow]')
        except KeyError:
            logging.warning(f'no {col} column. continue...')
        
    for col in INT_COLUMNS:
        logging.debug(f'converting {col} to category')
        try:
            df[col] = df[col].astype('uint32')    
        except KeyError:
            logging.warning(f'no {col} column. continue...')
            
    log_objectinfo(df, 'reads-df')
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
    
    
def load_mapseq_df( infile, fformat='reads', use_dask=False, chunksize=50000000):
    '''
    Abstracted loading code for all MAPseq pipeline dataframe formats. 
    
    '''
    FFORMATS = ['reads','aggregated','filtered','collapsed','readtable','vbctable']
    STR_COLS = {
        'reads'      : ['sequence'],
        'aggregated' : ['sequence'],
        'filtered'   : ['vbc_read', 'spikeseq', 'libtag', 'umi',  'ssi'],
        'collapsed'   : ['vbc_read_col','spikeseq', 'libtag', 'umi',  'ssi'],
        'readtable'  : ['vbc_read_col','spikeseq', 'libtag', 'umi',  'ssi'],
        'vbctable'   : ['vbc_read_col','spikeseq', 'libtag', 'umi',  'ssi'],      
        }
    
    INT_COLS = {
        'reads'      : [],
        'aggregated' : ['read_count'],
        'filtered'   : ['read_count'],
        'collapsed'   : ['read_count'],
        'readtable'  : ['read_count'],
        'vbctable'   : ['read_count','umi_count'],
        }

    CAT_COLS = {
        'reads'      : [],
        'aggregated' : [],
        'filtered'   : [],
        'collapsed'   : [],
        'readtable'  : ['label','site','type','brain','region'],
        'vbctable'   : ['label','site','type','brain','region'],        
        }
    
    logging.info(f'loading {infile} as format {fformat} use_dask={use_dask} chunksize={chunksize}')
    if fformat not in FFORMATS :
        logging.warning(f'fformat {fformat} not in {FFORMATS} . Datatypes may be incorrect.')
    
    ftype = None
    df = None
    
    if infile.endswith('.tsv'):
         ftype = 'tsv'
    elif infile.endswith('.parquet'):
       ftype = 'parquet'
    else:
        logging.error('input file must have relevant extension .tsv or .parquet')
        sys.exit(1)
    logging.debug(f'input filetype={ftype}')
    
    if ftype == 'tsv':
        if use_dask:
            df = dd.read_csv(infile, sep='\t')  
        else:
            df = pd.read_csv(infile, sep='\t', index_col=0)

        logging.debug(f'before: dtypes=\n{df.dtypes}')
        
        for col in STR_COLS[fformat]:
            logging.debug(f'converting {col} to string[pyarrow]')
            try:
                df[col] = df[col].astype('string[pyarrow]')
            except KeyError:
                logging.warning(f'no {col} column {fformat}? continue...')
            
        for col in INT_COLS[fformat]:
            logging.debug(f'converting {col} to integer')
            try:
                df[col] = df[col].astype('uint32')    
            except KeyError:
                logging.warning(f'no {col} column {fformat}? continue...')

        for col in CAT_COLS[fformat]:
            logging.debug(f'converting {col} to category')
            try:
                df[col] = df[col].astype('category')    
            except KeyError:
                logging.warning(f'no {col} column {fformat}? continue...')
      
        logging.debug(f'after: dtypes=\n{df.dtypes}')        
    
    elif ftype == 'parquet':
        if use_dask:
            df = dd.read_parquet(infile)  
        else:
            df = pd.read_parquet(infile) 
    return df       
    



def merge_fastq_pairs(config, readfilelist, outdir):
    logging.debug(f'processing {readfilelist}')
    if outdir is None:
        outdir = "."
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logging.debug(f'made outdir={outdir}')
    pairshandled = 0
    for (read1file, read2file) in readfilelist:
        logging.debug(f'read1file = {read1file}')
        logging.debug(f'read2file = {read2file}')
        pairshandled += 1    


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



def process_fastq_pairs_pd_chunked(infilelist, 
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
    logging.info(f'chunksize={chunksize} lines.')
    logging.debug(f'read1[{r1s}:{r1e}] + read2[{r2s}:{r2e}]')
    df = None
    sh = get_default_stats()
    pairnum = 1
    chunknum = 1

    for (read1file, read2file) in infilelist:
        logging.info(f'handling {read1file}, {read2file} ...')
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
            else:
                logging.debug(f'making additional read sequence DF...')
                ndf = pd.DataFrame(columns=['sequence'], dtype="string[pyarrow]")
                logging.info(f'got chunk len={len(chunk1)} slicing...')
                ndf['sequence'] = chunk1[0].str.slice(r1s, r1e) + chunk2[0].str.slice(r2s, r2e) 
                logging.info(f'chunk df type={ndf.dtypes}')                
                df = pd.concat([df, ndf], copy=False, ignore_index=True)
            logging.debug(f'handled chunk number {chunknum}')
            chunknum += 1

        logging.debug(f'handled pair number {pairnum}')
        pairnum += 1
    logging.debug(f'dtypes =\n{df.dtypes}')
    logging.info('Finished processing all input.')
    sh.add_value('/fastq','reads_handled', len(df) )
    sh.add_value('/fastq','pairs_handled', pairnum )
    return df          
        
def get_fh(infile):
    '''
    Opens compressed/uncompressed file as appropriate... 
    '''
    if infile.endswith('.gz'):
        fh = gzip.open(infile, "rt")
        logging.debug('handling gzipped file...')
    else:
        fh = open(infile, 'rt')
    return fh    


def filter_split_pd(df, 
                    max_repeats=None,
                    max_n_bases=None,
                    column='sequence',
                    drop=True,
                    cp=None
                    ):
    '''
    filter sequence by repeats, Ns. 
    split into MAPseq fields
    optionally drop sequence column
    
    '''
    if cp is None:
        cp = get_default_config()

    logging.info(f'Filtering by read quality. Repeats. Ns.')
    df = filter_reads_pd(df, 
                       max_repeats=max_repeats,
                       max_n_bases=max_n_bases, 
                       column=column )
    logging.info(f'Splitting into mapseq fields. ')    
    df = split_mapseq_fields(df, 
                             drop=True)

    return df

 
def split_mapseq_fields(df, column='sequence', drop=False, cp=None):
    '''
    spike_st=24
    spike_end = 32
    libtag_st=30
    libtag_end=32
    umi_st = 32
    umi_end = 44
    ssi_st = 44
    ssi_end = 52
    '''    
    logging.info(f'pulling out MAPseq fields...')
    if cp is None:
        cp = get_default_config()
    vbc_st = int(cp.get('mapseq', 'vbc_st'))
    vbc_end = int(cp.get('mapseq', 'vbc_end'))        
    spike_st=int(cp.get('mapseq', 'spike_st'))
    spike_end = int(cp.get('mapseq', 'spike_end'))
    libtag_st= int(cp.get('mapseq', 'libtag_st'))
    libtag_end= int(cp.get('mapseq', 'libtag_end'))
    umi_st = int(cp.get('mapseq', 'umi_st'))
    umi_end = int(cp.get('mapseq', 'umi_end'))
    ssi_st = int(cp.get('mapseq', 'ssi_st'))
    ssi_end = int(cp.get('mapseq', 'ssi_end'))
    
    df['vbc_read'] = df[column].str.slice(vbc_st,vbc_end).astype('string[pyarrow]')    
    df['spikeseq'] = df[column].str.slice(spike_st,spike_end).astype('string[pyarrow]')
    df['libtag'] = df[column].str.slice(libtag_st,libtag_end).astype('string[pyarrow]')    
    df['umi'] = df[column].str.slice(umi_st,umi_end).astype('string[pyarrow]')
    df['ssi'] = df[column].str.slice(ssi_st,ssi_end).astype('string[pyarrow]')
    sh = get_default_stats()
    
    if drop:
        logging.info(f'dropping {column} column to slim.')
        df.drop( column, axis=1, inplace=True)
    else:
        logging.info(f'drop is false. keeping {column} column.')
    logging.info(f'df done. len={len(df)} returning...')
    # changes in-place, but return as well. 
    return df


def filter_reads_pd(df,
                    max_repeats=None,
                    max_n_bases=None,
                    min_reads=None, 
                    column='sequence',
                    remove=True,
                    cp=None):
    '''

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





            
      


def filter_counts_fasta(cp, infile, 
                        outfile, 
                        datestr=None, 
                        min_count = 2 ):
    '''
    we want to go FASTA to FASTA, but filter by read counts.  
    '''
    sdf = read_fasta_to_df(infile)
    logging.debug(f'got {len(sdf)} sequences.')
    fdf = filter_counts_df(cp, sdf, min_count)
    logging.info(f'got {len(fdf)} filtered sequences. Writing to {outfile} FASTA.')
    write_fasta_from_df(fdf, outfile)
    logging.debug(f'Wrote {outfile}')


def set_counts_df(seqdf, column='', cp=None):
    '''
    assumed sequence column 
    calculates duplicates and sets read_count column
 
    '''
    initlen = len(seqdf)
    logging.debug(f'setting read counts for sequence DF len={len(seqdf)}')
    cdf = make_read_counts_df(cp, seqdf, label=None)
    fmap = cdf.set_index('sequence')['read_count'].to_dict()
    seqdf['read_count'] = seqdf['sequence'].map(fmap)
    # change made to inbound df, but return anyway
    return seqdf

def set_siteinfo(df, sampdf, column='sequence', cp=None):
    '''
    This is a purely utility function to determine site type for
    shoulder plots. 
    
    '''
    if cp is None:
        cp=get_default_config()
    ssi_st = int(cp.get('mapseq','ssi_st') )
    ssi_end = int(cp.get('mapseq','ssi_end') )
    bcfile = os.path.expanduser( cp.get('barcodes','ssifile'))

    df['ssi'] = df[column].str.slice(ssi_st,ssi_end).astype('string[pyarrow]')
    df.drop( column, axis=1, inplace=True)

    # set site
    # map SSIs, set unknown to unmatched.
    logging.debug(f'getting rt labels...')
    labels = get_rtlist(sampdf)
    logging.debug(f'rtlabels={labels}')
    #bcdict = get_barcode_dict(bcfile, labels)
    rtdict = get_rtprimer_dict(bcfile, labels)
   
    #logging.info('filling in labels by SSI sequences...')
    #df['label'] = df['ssi'].map(bcdict)
    #logging.info('labelling unmatched...')
    #df.fillna({'label': 'nomatch'}, inplace=True)

    logging.info('filling in rtprimer number by SSI sequences...')    
    df['rtprimer'] = df['ssi'].map(rtdict)
    logging.info('labelling unmatched...')
    df.fillna({'rtprimer': 'nomatch'}, inplace=True)
    
    sdf = sampdf[['rtprimer','siteinfo']]
    sdf = sdf[sdf['siteinfo'] != '']
    smap = dict(zip(sdf['rtprimer'],sdf['siteinfo']))
    df['site'] = df['rtprimer'].map(smap)
    df.fillna({'site': 'nomatch'}, inplace=True)
    sdf = None    
    
    return df
    




def aggregate_reads_pd(seqdf, pcolumn='sequence'):
    initlen = len(seqdf)
    logging.debug(f'collapsing with read counts for sequence DF len={len(seqdf)}')
    ndf = seqdf.value_counts()
    ndf = ndf.reset_index()
    ndf.rename({'count':'read_count'}, inplace=True, axis=1)
    return ndf

def aggregate_reads_dd(seqdf, column='sequence', outdir=None, min_reads=2, chunksize=50000000, dask_temp=None):
    '''
    ASSUMES INPUT IS DASK DATAFRAME
    
    '''
    if dask_temp is not None:
        logging.info(f'setting dask temp to {os.path.expanduser(dask_temp)} ')
        dask.config.set(temporary_directory=os.path.expanduser(dask_temp))
    else:
        logging.info(f'no dask_temp specified. letting dask use its default')
        
    initlen = len(seqdf)
    logging.debug(f'collapsing with read counts for col={column} len={len(seqdf)}')
    ndf = seqdf[column].value_counts().compute()
    ndf = ndf.reset_index()
    ndf.rename({'count':'read_count'}, inplace=True, axis=1)
  
    if min_reads > 1:
        logging.debug(f'Dropping reads with less than {min_reads} read_count.')
        logging.info(f'Length before read_count threshold={len(ndf)}')
        ndf = ndf[ndf['read_count'] >= min_reads]
        ndf.reset_index(inplace=True, drop=True)
        logging.info(f'Length after read_count threshold={len(ndf)}')    
    else:
        logging.info(f'min_reads = {min_reads} skipping initial read count thresholding.')  
    return ndf

    
def filter_counts_df(cp, countsdf, min_count):
    '''
    Assumes read_count column
    '''      
    countsdf = countsdf[countsdf['read_count'] >= min_count]
    countsdf.reset_index(drop=True, inplace=True)
    logging.info(f'filtering by count: before={initlen} after={len(seqdf)} min_count={min_count} ')
    return cdf    
            



def make_read_counts_dfs(config, filelist, outdir):
    '''
    
    '''
    dflist = []
    for filepath in filelist:
        logging.info(f'calculating counts for file {filepath} ...')    
        dirname = os.path.dirname(filepath)
        
        if outdir is not None:
            dirname = outdir
        
        filename = os.path.basename(filepath)
        (base, ext) = os.path.splitext(filename)   
        logging.debug(f'handling {filepath} base={base}')
        
        # make raw fasta TSV of barcode-splitter output for one barcode. 
        # trim to 44 unique w/ counts. 
        seqdf = make_fasta_df(config, filepath)
        of = os.path.join(dirname , f'{base}.read.seq.tsv')
        seqdf.to_csv(of, sep='\t')
        
        # to calculate threshold we need counts calculated. 
        cdf = make_read_counts_df(config, seqdf, label=base)  
        logging.debug(f'made counts df: {len(cdf)} reads.')
        of = os.path.join(dirname , f'{base}.read.counts.tsv')
        cdf.to_csv(of, sep='\t')
        dflist.append(cdf)
    logging.debug(f'returning list of {len(dflist)} counts DFs...')
    return dflist 
        

def make_clustered_heatmap(df, outprefix, columns=None ):
    '''
    
    Caller should edit columns in order to exclude injection areas from plot. 
    '''
    camp = 'Reds'
    g = sns.clustermap(df, cmap=camp, yticklabels=False, col_cluster=False, standard_scale=0)
    g.fig.subplots_adjust(right=0.7)
    g.ax_cbar.set_position((0.8, .2, .03, .4))
    plt.title(f'{prefix}\nCounts')
    plt.savefig(f'{outprefix}.heatmap.pdf')
    logging.info(f'done making {outprefix}.heatmap.pdf ')
    

def make_read_countplot(config, df, outfile, title=None ): 
    '''
    makes individual read count plot from sequence read_count DF 
    assumes 'sequence' and 'read_count' columns. 
    
    '''   
    plt.set_loglevel (level = 'warning')
    
    logging.debug(f'handling sequence df len={len(df)}')
    outfile = os.path.abspath(outfile)    

    if title is None:
        title='Read count frequence plot.'

    df.sort_values(by='read_count', ascending=False, inplace=True)
    df.reset_index(inplace=True)
    df['index'] = df['index'].astype('int64')

    plt.figure()
    plt.plot(np.log10(df['index']), np.log10(df['read_count']))
    plt.title(title)
    plt.xlabel("log10(index)")
    plt.ylabel("log10(read_count)")
    logging.info(f'Saving count plot to {outfile}')
    plt.savefig(outfile)


def make_read_countplot_sns(cp, df, outfile='count-frequency-plot.pdf' ):
    '''
    makes two figures, one log-normalized, one straight for same counts df. 
    '''
    from matplotlib.backends.backend_pdf import PdfPages as pdfpages
    df.sort_values(by='read_count', ascending=False, inplace=True)
    df.reset_index(inplace=True)
    df['index'] = df.index.astype('int64')

    #page_dims = 
    #  A4 landscape   (11.69,8.27) 
    #  A3 landscape  (16.53,11.69)
    #  A4 portrait  (8.27, 11.69)  
    page_dims = (8.27, 11.69)
    with pdfpages(outfile) as pdfpages:
        axlist = []
        fig,axes = plt.subplots(nrows=2, ncols=1, figsize=page_dims, layout='constrained') 
        fig.suptitle(f'Read counts frequency plots.')
        for a in axes.flat:
            axlist.append(a)
        ax = axlist[0]
        counts_axis_plot_sns(ax, df, scale=None)
        ax = axlist[1]
        counts_axis_plot_sns(ax, df, scale='log10')
        
        pdfpages.savefig(fig)



def make_shoulder_plot_sns(df, title='counts frequency plot', site=None, outfile='count-frequency-plot.pdf' ):
    '''
    makes two figures, one log-normalized, one straight for same counts df. 
    '''
    from matplotlib.backends.backend_pdf import PdfPages as pdfpages
    
    cdf = pd.DataFrame( df[df['site'] == site ]['read_count'].copy(), columns=['read_count'])
    cdf.sort_values(by='read_count', ascending=False, inplace=True)
    cdf.reset_index(inplace=True,  drop=True)
    cdf['index'] = cdf.index.astype('int64')

    #page_dims = 
    #  A4 landscape   (11.69,8.27) 
    #  A3 landscape  (16.53,11.69)
    #  A4 portrait  (8.27, 11.69)  
    page_dims = (8.27, 11.69)
    with pdfpages(outfile) as pdfpages:
        axlist = []
        fig,axes = plt.subplots(nrows=2, ncols=1, figsize=page_dims, layout='constrained') 
        fig.suptitle(f'Read counts freq: site={site}')
        for a in axes.flat:
            axlist.append(a)
        ax = axlist[0]
        counts_axis_plot_sns(ax, cdf, scale=None)
        ax = axlist[1]
        counts_axis_plot_sns(ax, cdf, scale='log10')        
        pdfpages.savefig(fig)

        
def counts_axis_plot_sns(ax, df, scale=None ) :
    '''
    Creates individual axes for single plot within figure. 
    scale = None | log10  | log2
    '''
    s = df['read_count'].sum()
    n = len(df)
    t = df['read_count'].max()
    r = df['index'].max()
    
    title=''
    if scale is None:
        h = calc_freq_threshold(df, fraction=0.9, column = 'read_count')
        sns.lineplot(ax=ax, x=df['index'], y=df['read_count'] )
        title = 'counts frequency plot'
        #lx = 0.2 * t 
        #ly = 0.2 * r
        lx = 1.0
        ly = 1.0
        ax.text(lx, ly, f"n={n}\ntop={t}\nsum={s}\nestimated_threshold={h}", fontsize=11) #add text
        logging.debug(f'made axis without scale.') 
    elif scale == 'log10':
        # avoid divide by 0 runtime warning...
        df['log10index'] = np.log10(df['index'] + 1)
        df['log10counts'] = np.log10(df['read_count'])
        sns.lineplot(ax=ax, x=df['log10index'], y=df['log10counts'] )
        #lx = 0.05 * np.log10(t) 
        #ly = 0.05 * np.log10(r)        
        lx = 0.2
        ly = 0.2
        ax.text(lx, ly, f"n={n}\ntop={t}\nsum={s}\n", fontsize=11) #add text
        title = 'log10(counts) frequency plot'
        logging.debug(f'made axis with log10 scale.')       
    ax.set_title(title, fontsize=10)


def normalize_weight(df, weightdf, columns=None):
    '''   
    Weight values in realdf by spikedf, by column groups
    Assumes matrix index is sequence.
    Assumes matrices have same columns!!  
    If column numbers are mis-matched, will create empty column
    If columns is none, use/weight all columns
    
    !! only use target columns. 
    
    '''
    logging.debug(f'normalizing df=\n{df}\nby weightdf=\n{weightdf}')
    
    # sanity checks, fixes. 
    if len(df.columns) != len(weightdf.columns):
        logging.error(f'mismatched matrix columns df:{len(df.columns)} weightdf: {len(weightdf.columns)} !!')
        
    #which SSI has highest spikein?
    sumlist = []
    for col in weightdf.columns:
        sum = weightdf[col].sum()
        sumlist.append(sum)
    sum_array = np.array(sumlist)
    maxidx = np.argmax(sum_array)
    maxval = sum_array[maxidx]  
    maxcol = weightdf.columns[maxidx]
    logging.debug(f'largest spike sum for {maxcol} sum()={maxval}')
    factor_array =  maxval / sum_array
    logging.debug(f'factor array= {list(factor_array)}')

    max_list = []
    sum_list = []
    for col in df.columns:
        max_list.append(df[col].max())
        sum_list.append(df[col].sum())
    logging.debug(f'real max_list={max_list}')
    logging.debug(f'real sum_list={sum_list}')
    
    normdf = df.copy()
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

    return normdf
       

def normalize_weight_grouped(df, weightdf, columns=None):
    '''
     Weight values in df by weightdf, by column groups
    
    
    df   dataframe to be scaled, weightdf is spikeins
    columns is list of column groupings. weighting will done within the groups.     
   
    Assumes matrix index is sequence.
    Assumes matrices have same columns  
    If column numbers are mis-matched, will create empty column
    If columns is none, use/weight all columns
    Will only return columns somewhere in columns. 

    '''
    logging.debug(f'normalizing df=\n{df}\nby weightdf=\n{weightdf}')
    
    # sanity checks, fixes. 
    if len(df.columns) != len(weightdf.columns):
        logging.error(f'mismatched matrix columns df:{len(df.columns)} weightdf: {len(weightdf.columns)} !!')
    
    if columns is None:
        columns = [ list( df.columns ) ]    
    
    outdf = None
    
    for i, group in enumerate( columns ):
        logging.debug(f'handling group {i} of {len(group)} columns.')
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
    
    logging.info(f'normalizing done.len={len(outdf)}. cols={len(outdf.columns)} ')
    return outdf


def normalize_scale(df, columns = None, logscale='log2', min=0.0, max=1.0):
    '''
    Log scale whole matrix.   log10 or log2 ???
    Set -inf to 0
    Set NaN to 0
    
    
    '''
    #logging.debug(f'making rows sum to one...')
    if logscale == 'log2':
        ldf = np.log2(df)
    elif logscale == 'log10':
        ldf = np.log10(df)
    
    for c in ldf.columns:
        ldf[c] = np.nan_to_num(ldf[c], neginf=0.0)
    return ldf


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


def filter_non_injection(rtdf, ridf, min_injection=1):
    '''
    rtdf and ridf should already be filtered by brain, type, and anything else that might complicate matters.
    remove rows from rtdf that do not have at least <min_injection> value in the row 
    of ridf with the same index (VBC sequence)
    Does an inner join() on the dataframes, keyed on sequence. 
    Keeps values and columns from first argument (rtdf)
    
    '''
    logging.info(f'filtering non-injection. min_injection={min_injection}')
    logging.debug(f'before threshold inj df len={len(ridf)}')
    ridf = ridf[ridf.umi_count >= min_injection]
    ridf.reset_index(inplace=True, drop=True)
    logging.debug(f'before threshold inj df len={len(ridf)}')   
    
    
    mdf = pd.merge(rtdf, ridf, how='inner', left_on='sequence', right_on='sequence')
    incol = mdf.columns
    outcol = []
    selcol =[]
    for c in incol:
        if not c.endswith('_y'):
            selcol.append(c)
            outcol.append(c.replace('_x',''))
    mdf = mdf[selcol]
    mdf.columns = outcol
    logging.debug(f'created merged/joined DF w/ common sequence items.  df=\n{mdf}')
    return mdf


def filter_non_inj_umi(rtdf, ridf, inj_min_umi=1, write_out=False):
    '''
    rtdf and ridf should already be filtered by brain, type, and anything else that might complicate 
    matters.
    remove rows from rtdf that do not have at least <min_injection> value in the row 
    of ridf with the same index (VBC sequence)
    Does an inner join() on the dataframes, keyed on sequence. 
    Keeps values and columns from first argument (rtdf)
    
    '''
    logging.info(f'filtering non-injection. inj_min_umi={inj_min_umi}')
    logging.debug(f'before threshold inj df len={len(ridf)}')
    ridf = ridf[ridf.umi_count >= inj_min_umi]
    ridf.reset_index(inplace=True, drop=True)
    
    logging.debug(f'before threshold inj df len={len(ridf)}')   
    if write_out:
        ridf.to_csv('./ridf.tsv', sep='\t')
        rtdf.to_csv('./rtdf.tsv', sep='\t')

    # get target VBCs that are in injection
    mtdf = merge_and_filter(rtdf, ridf)

    # get injection VBCs that are in at least one target, similarly 
    midf = merge_and_filter(ridf, rtdf)
    
    return ( mtdf, midf)




def merge_and_filter(adf,  bdf, on='vbc_read_col', indicator=True, how='outer'):
    '''
    Return filtered ADF consisting of only entries in ADF that also exist in BDF as joined on 
    the 'on' column.  
    
    
    '''
    #suffixes=('_a','_b')
    # get targets in injection merge using vbc_read_col as common field. 
    #mdf = rtdf.merge(ridf, on='vbc_read_col', indicator=True, how='outer', suffixes=('_t','_i'))
    mdf = adf.merge(bdf, on='vbc_read_col', indicator=True, how='outer', suffixes=('_a','_b'))
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

def process_make_readtable_pd(df,
                          sampdf,
                          bcfile=None, 
                          outdir=None, 
                          cp = None):
    '''
    take split/aligned read file.
    produce fully tagged and filtered data table.
    -- convert SSI column to BCXXX tag. 
    -- add optional region label
    -- classify by site-type
    -- set brain label
    
    aggregate and sum read and UMI counts by type, libtag, and ssi. 
    
    '''
    logging.info(f'inbound df len={len(df)} columns={df.columns}')
    if outdir is None:
        outdir = os.path.abspath('./')
    outdir = os.path.abspath(outdir)    
    os.makedirs(outdir, exist_ok=True)
    
    if cp is None:
        cp = get_default_config()
    
    spikeseq = cp.get('readtable','spikeseq')
    realregex = cp.get('readtable', 'realregex' )
    loneregex = cp.get('readtable', 'loneregex' )
    
    if bcfile is None:
        bcfile = cp.get('barcodes','ssifile')    
    
    logging.debug(f'spikeseq={spikeseq} realregex={realregex} loneregex={loneregex} bcfile={bcfile}')
    
    # map SSIs, set unknown to unmatched.
    logging.debug(f'getting rt labels...')
    labels = get_rtlist(sampdf)
    logging.debug(f'rtlabels={labels}')
    bcdict = get_barcode_dict(bcfile, labels)
    rtdict = get_rtprimer_dict(bcfile, labels)
    logging.debug(f'got {len(bcdict)} barcodes with labels and primer number.')    
    logging.info('filling in labels by SSI sequences...')
    df['label'] = df['ssi'].map(bcdict)
    logging.info('labelling unmatched...')
    df.fillna({'label': 'nomatch'}, inplace=True)

    logging.info('filling in rtprimer number by SSI sequences...')    
    df['rtprimer'] = df['ssi'].map(rtdict)
    logging.info('labelling unmatched...')
    df.fillna({'rtprimer': 'nomatch'}, inplace=True)
    
    # spikeins
    # make temporary column to check spikein. 
    # we have to spikes LAST, because all spikes are real L2s, but not all reals are spikes.
   
    logging.info('identifying reals by libtag...')
    rmap = df['libtag'].str.match('[TC][TC]$')
    df.loc[rmap, 'type'] = 'real'
    
    logging.info('identifying L1s by libtag...')
    lmap = df['libtag'].str.match('[AG][AG]$')
    df.loc[lmap, 'type'] = 'lone'

    logging.info('identifying spikeins by spikeseq matching...')
    smap = df['spikeseq'] == spikeseq
    df.loc[smap, 'type'] = 'spike'
    
    df.fillna({'type': 'nomatch'}, inplace=True)
    df.drop(['spikeseq'], inplace=True, axis=1)

    # set brain
    bdf = sampdf[['rtprimer','brain']]
    bdf = bdf[bdf['brain'] != '']
    bmap = dict(zip(bdf['rtprimer'],bdf['brain']))
    df['brain'] = df['rtprimer'].map(bmap)
    df.fillna({'brain': ''}, inplace=True)
    bdf = None
    
    # set region
    rdf = sampdf[['rtprimer','region']]
    rdf = rdf[rdf['region'] != '']
    rmap = dict(zip(rdf['rtprimer'],rdf['region']))
    df['region'] = df['rtprimer'].map(rmap)
    df.fillna({'region': ''}, inplace=True)
    rdf = None
    
    # set site
    sdf = sampdf[['rtprimer','siteinfo']]
    sdf = sdf[sdf['siteinfo'] != '']
    smap = dict(zip(sdf['rtprimer'],sdf['siteinfo']))
    df['site'] = df['rtprimer'].map(smap)
    df.fillna({'site': ''}, inplace=True)
    sdf = None    
    
    # calc/ collect stats
    sh = get_default_stats()    
    sh.add_value('/readtable','n_full_sequences', str(len(df)) )
    
    #n_spike = (df['type'] == 'spike')['read_count'].sum() 
    #n_lone =  (df['type'] == 'lone')['read_count'].sum()  
    #n_real =  (df['type'] == 'real')['read_count'].sum() 
    #n_nomatch = (df['type'] == 'nomatch')['read_count'].sum()

    # get rid of nans, and calculate template switching value. 
    ndf = df.replace('nomatch',np.nan)
    ndf.dropna(inplace=True, axis=0, ignore_index=True)

    # find and remove ( at least) known template-switch rows from dataframe. 
    tsdf = ndf[ ((ndf['type'] == 'lone') & (ndf['site'].str.startswith('target')) & (ndf['site'] != 'target-control'))]
    of = os.path.join(outdir, 'template_switch.tsv') 
    logging.info(f'Writing template switch DF len={len(df)} Writing to {of}')
    tsdf.to_csv(of, sep='\t')
    n_tswitch = len(tsdf)
    
    ndf.drop(tsdf.index, inplace=True)
    ndf.reset_index(drop=True, inplace=True)
     
    # label tswitch tswitch and filter from final...
    # 

    sh.add_value('/readtable', 'n_tswitch', str(n_tswitch) )
    #sh.add_value('/readtable', 'n_spike', str(n_spike) )     
    #sh.add_value('/readtable', 'n_lone', str(n_lone) )      
    #sh.add_value('/readtable', 'n_real', str(n_real) )     
    #sh.add_value('/readtable', 'n_nomatch', str(n_nomatch) )    

    return df


def make_shoulder_plots(df, outdir=None, cp=None):
    # make shoulder plots. injection, target
    logging.info('making shoulder plots...')
    if outdir is None:
        outdir = os.path.abspath('./')
    logging.getLogger('matplotlib.font_manager').disabled = True
    if len(df[df['site'] == 'injection'] ) > 1:
        make_shoulder_plot_sns(df, site='injection', outfile=f'{outdir}/inj-counts.pdf')
    else:
        logging.info(f'no injection sites, so no plot.')
    if len(df[df['site'] == 'target'] ) > 1:    
        make_shoulder_plot_sns(df, site='target', outfile=f'{outdir}/target-counts.pdf')   
    else:
        logging.info(f'no target sites, so no plot.')
   

def process_make_vbctable_pd(df,
                          outdir=None,
                          inj_min_reads = 3,
                          target_min_reads = 2, 
                          cp = None):
    '''   
    
    primary columns:  label, umi, type  
    derived columns:  brain, region, rtprimer, label, ssi
    -- remove nomatches
    -- threshold read_count by site type
    -- collapse by VBC sequence, calculating UMI count
    
    This should  be prepared for input it make VBC matrices    
    '''
    logging.debug(f'inbound df len={len(df)} columns={df.columns}')
    log_objectinfo(df, 'input-df')
    ndf = df
    
    # remove unneeded columns
    #  libtag:     we now know real/lone
    #  ssi:        we have assigned label
    #  rtprimer:     "  " 
    #  vbc_read:   we have already collapsed
    logging.info(f'dropping redundant readtable columns: libtag ssi rtprimer vbc_read  ')
    ndf.drop(labels=['libtag','ssi','rtprimer'], axis=1, inplace=True)   
    ndf.replace('nomatch',np.nan, inplace=True)
    #ndf = df.replace('nomatch' ,np.nan)
    logging.info(f'DF before removing nomatch: {len(ndf)}')
    log_objectinfo(ndf, 'new-df')
    ndf.dropna(inplace=True, axis=0, ignore_index=True)
    logging.info(f'DF after removing nomatch/NaN: {len(ndf)}')
    log_objectinfo(ndf, 'new-df')
    
    # threshold by read counts, by siteinfo, before starting...
    logging.info(f'thresholding by read count. inj={inj_min_reads} target={target_min_reads} len={len(ndf)}') 
    tdf = ndf[ndf['site'].str.startswith('target')]
    tdf = tdf[tdf['read_count'] >= int(target_min_reads)]
    idf = ndf[ndf['site'].str.startswith('injection')]
    idf = idf[idf['read_count'] >= int(inj_min_reads)]    
    thdf = pd.concat([tdf, idf])
    thdf.reset_index(drop=True, inplace=True)
    logging.info(f'DF after threshold inj={inj_min_reads} tar={target_min_reads}: {len(thdf)}')
    log_objectinfo(thdf, 'threshold-df')    
    udf = thdf.groupby(['vbc_read_col','label','type']).agg( {'umi' : 'nunique',
                                                              'read_count':'sum', 
                                                              'brain':'first',
                                                              'region':'first',
                                                              'site':'first'
                                                              } ).reset_index()
    
    udf.rename( {'umi':'umi_count'}, axis=1, inplace=True)
    logging.info(f'DF after umi/label collapse: {len(udf)}')
    
    sh = get_default_stats()
    sh.add_value('/vbctable','n_vbcs', len(udf) )    
    
    log_objectinfo(udf, 'umi-df')
    return udf


def process_make_matrices_pd(df,
                          outdir=None,
                          exp_id = 'M001',  
                          inj_min_umi = None,
                          target_min_umi = None,
                          target_min_umi_absolute = None, 
                          label_column='label',
                          cp = None):
    '''
    -- per brain, pivot real VBCs (value=umi counts)
    -- create real, real normalized by spikein, and  
    -- use label_column to pivot on, making it the y-axis, x-axis is vbc sequence.  
    -- optionally require VBCs to be in injection to be included
    -- optionally include injections in matrices.
    
    theshold logic. 
    inj_min_umi                VBC UMI must exceed to be kept.
    target_min_umi             if ANY target area exceeds, keep all of that VBC targets. 
    target_min_umi_absolute    hard threshold cutoff
 

    '''
    if cp is None:
        cp = get_default_config()
    if outdir is None:
        outdir = os.path.abspath('./')
    require_injection = cp.getboolean('matrices','require_injection')
    include_injection = cp.getboolean('matrices','include_injection')

    max_negative = 1
    max_water_control = 1

    if inj_min_umi is None:
        inj_min_umi = int(cp.get('matrices','inj_min_umi'))
    if target_min_umi is None:
        target_min_umi = int(cp.get('matrices','target_min_umi'))   
    if target_min_umi_absolute is None:
        target_min_umi_absolute = int(cp.get('matrices','target_min_umi_absolute'))

    use_target_negative=cp.getboolean('matrices','use_target_negative')
    use_target_water_control=cp.getboolean('matrices','use_target_water_control')    
    clustermap_scale = cp.get('matrices','clustermap_scale')
    logging.debug(f'running exp_id={exp_id} inj_min_umi={inj_min_umi} target_min_umi={target_min_umi} use_target_negative={use_target_negative} use_target_water_control={use_target_water_control}')

    df['brain'] = df['brain'].astype('string')
    bidlist = list(df['brain'].dropna().unique())
    bidlist = [ x for x in bidlist if len(x) > 0 ]
    bidlist.sort()
    logging.debug(f'handling brain list: {bidlist}')
    
    sh = get_default_stats()

    norm_dict = {}

    for brain_id in bidlist:
        valid = True
        logging.debug(f'handling brain_id={brain_id}')
        bdf = df[df['brain'] == brain_id]

        # extract and threshold injection areas
        idf = bdf[bdf['site'].str.startswith('injection')]
        ridf = idf[idf['type'] == 'real'] 
        if inj_min_umi > 1:
            before = len(ridf)
            ridf = filter_all_lt(ridf, 'vbc_read_col', 'umi_count', inj_min_umi)            
            #
            # INCORRECT: if require_injection is false, there may not be any injection but thats OK.
            # Let the filter_non_inj_umi() function do its thing...
            #
            #if not len(ridf) > 0:
            #    valid = False
            #    logging.warning(f'No injection VBCs passed min_inj_umi filtering! Skip brain.')
            logging.debug(f'filtering by inj_min_umi={inj_min_umi} before={before} after={len(ridf)}')
        else:
            logging.debug(f'inj_min_umi={inj_min_umi} no filtering.')
                   
        # extract and do absolute thesholding of target areas
        tdf = bdf[bdf['site'].str.startswith('target')]
        rtdf = tdf[tdf['type'] == 'real'] 

        if target_min_umi_absolute > 1:
            before = len(rtdf)
            rtdf = rtdf[rtdf['umi_count'] >= target_min_umi_absolute ]
            rtdf.reset_index(drop=True, inplace=True)
            after = len(rtdf)
            logging.debug(f'filtering by target_min_umi_absolute={target_min_umi_absolute} before={before} after={after}')

        # threshold by min_target or threshold by target-negative
        # if use_target_negative is true, but no target negative site 
        # defined, use min_target and throw warning. 
        if use_target_negative:
            logging.info(f'use_target_negative is {use_target_negative}')
            max_negative = calc_min_umi_threshold(bdf, 'target-negative', cp)
            logging.debug(f'target-negative UMI count = {max_negative}')

        if use_target_water_control:
            logging.info(f'use_target_water_control is {use_target_water_control}')
            max_water_control = calc_min_umi_threshold(bdf, 'target-water-control',cp)
            logging.debug(f'target_water_control UMI count = {max_water_control}')

        target_min_umi = max([target_min_umi, max_negative, max_water_control ])
        logging.debug(f'min_target UMI count after all constraints = {target_min_umi}')   

        # min_target is now either calculated from target-negative, or from config. 
        if target_min_umi > 1:
            before = len(rtdf)
            frtdf = filter_all_lt(rtdf, 'vbc_read_col', 'umi_count', target_min_umi)            
            if not len(frtdf) > 0:
                valid = False
                logging.warning(f'No VBCs passed min_target filtering! Skip brain.')
            logging.debug(f'filtering by target_min_umi={target_min_umi} before={before} after={len(rtdf)}')
        else:
            logging.debug(f'target_min_umi={target_min_umi} no filtering.')
            frtdf = rtdf
  
        # get injection-filtered real target table, and target-filtered real injection table
        # in case either is needed. 
        (ifrtdf, tfridf) = filter_non_inj_umi(frtdf, ridf, inj_min_umi=inj_min_umi)            
        logging.debug(f'{len(ifrtdf)} real target VBCs after injection filtering.')
        logging.debug(f'{len(tfridf)} real injection VBCs after target filtering.')

        if require_injection:
            logging.debug(f'require_injection={require_injection} inj_min_umi={inj_min_umi}')
            if len(ridf) == 0:
                logging.warning('require_injection=True but no real VBCs from any injection site.')
                valid = False
            elif not len(ifrtdf) > 0:
                logging.warning(f'No VBCs passed injection filtering! Skip brain.')
                valid = False
            else:
                # set filtered real to injection-filtered real targets. 
                frtdf = ifrtdf    
        else:
            # frtdf remains unfiltered...
            logging.debug(f'require_injection={require_injection} proceeding...')

        vbcdf = frtdf
        
        if include_injection:
            logging.debug(f'include_injection={include_injection} including mutually present injection VBCs') 
            vbcdf = pd.concat( [frtdf, tfridf], ignore_index=True ) 


        # make matrices if per-brain data is valid... 
        if valid:                               
            vbcdf.to_csv(f'{outdir}/{brain_id}.vbctable.tsv', sep='\t')
            vbcdf.to_parquet(f'{outdir}/{brain_id}.vbctable.parquet')
            
            target_columns = list( vbcdf[vbcdf['site'].str.startswith('target')]['label'].unique())
            injection_columns = list( vbcdf[vbcdf['site'].str.startswith('injection')]['label'].unique())
            
            rbcmdf = vbcdf.pivot(index='vbc_read_col', columns=label_column, values='umi_count')
            scol = natsorted(list(rbcmdf.columns))
            rbcmdf = rbcmdf[scol]
            rbcmdf.fillna(value=0, inplace=True)
            logging.debug(f'brain={brain_id} raw real barcode matrix len={len(rbcmdf)}')            
            
            # filtered reals
            fbcmdf = vbcdf.pivot(index='vbc_read_col', columns=label_column, values='umi_count')
            scol = natsorted(list(fbcmdf.columns))
            fbcmdf = fbcmdf[scol]
            fbcmdf.fillna(value=0, inplace=True)
            logging.debug(f'brain={brain_id} real barcode matrix len={len(fbcmdf)}')

            # spikes
            sdf = tdf[tdf['type'] == 'spike']
            if include_injection: 
                logging.debug(f'include_injection={include_injection} need injection spikes...')
                isdf = idf[idf['type']== 'spike']
                sdf = pd.concat( [sdf, isdf], ignore_index=True  )
                
            sbcmdf = sdf.pivot(index='vbc_read_col', columns=label_column, values='umi_count')
            spcol = natsorted(list(sbcmdf.columns))
            sbcmdf = sbcmdf[spcol]
            sbcmdf.fillna(value=0, inplace=True)    
            logging.debug(f'brain={brain_id} spike barcode matrix len={len(sbcmdf)}')
    
            (fbcmdf, sbcmdf) = sync_columns(fbcmdf, sbcmdf)
            
            if include_injection:
                logging.debug(f'include_injection={include_injection} need grouped normalization.')
                nbcmdf = normalize_weight_grouped(fbcmdf, sbcmdf, columns = [target_columns, injection_columns])
            else:
                logging.debug(f'include_injection={include_injection} ungrouped normalization.')
                nbcmdf = normalize_weight_grouped(fbcmdf, sbcmdf)
            
            # put columns in natural sort order...
            nbcmdf = nbcmdf[ natsorted( list(nbcmdf.columns) ) ]            
            logging.debug(f'nbcmdf.describe()=\n{nbcmdf.describe()}')
            
            norm_dict[brain_id] = nbcmdf
            
            # make scaled normalized matrix
            scbcmdf = normalize_scale(nbcmdf, logscale=clustermap_scale)
            scbcmdf.fillna(value=0, inplace=True)
            logging.debug(f'scbcmdf.describe()=\n{scbcmdf.describe()}')
            
            sh.add_value(f'/matrices/brain_{brain_id}','valid', 'True' )            
            sh.add_value(f'/matrices/brain_{brain_id}','n_vbcs_real', len(rbcmdf) )
            sh.add_value(f'/matrices/brain_{brain_id}','n_vbcs_real_filtered', len(fbcmdf) )
            sh.add_value(f'/matrices/brain_{brain_id}','n_vbcs_spike', len(sbcmdf) )
            
            # raw real
            rbcmdf.to_csv(f'{outdir}/{brain_id}.rbcm.tsv', sep='\t')
            # filtered real
            fbcmdf.to_csv(f'{outdir}/{brain_id}.fbcmdf.tsv', sep='\t')
            # spike-in matrix
            sbcmdf.to_csv(f'{outdir}/{brain_id}.sbcm.tsv', sep='\t')
            # filtered normalized by spike-ins.    
            nbcmdf.to_csv(f'{outdir}/{brain_id}.nbcm.tsv', sep='\t')
            # log-scaled normalized 
            scbcmdf.to_csv(f'{outdir}/{brain_id}.scbcm.tsv', sep='\t')
                       
        else:
            logging.info(f'brain {brain_id} data not valid.')
            sh.add_value(f'/matrices/brain_{brain_id}','valid', 'False' )
            
        logging.info(f'done with brain={brain_id}')    
    logging.info(f'got dict of {len(norm_dict)} normalized barcode matrices. returning.')
    return norm_dict


def process_mapseq_all_XXX(infiles, sampleinfo, bcfile=None, outdir=None, expid=None, cp=None):    
    '''
    DO NOT USE. NEEDS UPDATING FOR LARGE DATA
    
    performs end-to-end default processing. fastq pairs to matrices
    
    '''
    if cp is None:
        cp = get_default_config()
    
    if expid is None:
        expid = 'EXP'    
    
    if bcfile is None:
        bcfile = cp.get('barcodes','ssifile')
    
    if outdir is None:
        outdir = os.path.abspath('./')
        logging.debug(f'outdir not provided, set to {outdir}')
    else:
        outdir = os.path.abspath(outdir)
        logging.debug(f'outdir = {outdir} ')   
   
    logging.debug(f'ensuring outdir: {outdir} ')
    os.makedirs(outdir, exist_ok=True)
    logging.info(f'outdir={outdir}')

   
    # process_fastq_pairs
    if (len(infiles) < 2)  or (len(infiles) % 2 != 0 ):
        parser.print_help()
        print('error: the following arguments are required: 2 or multiple of 2 infiles')
        sys.exit(1)
    infiles = package_pairfiles(infiles) 
    
    outfile = f'{outdir}/{expid}.reads.tsv'
    df = process_fastq_pairs_pd(infiles, 
                                 outdir,  
                                 cp=cp)    

    logging.debug(f'filtering by read quality. repeats. Ns.')
    df = filter_reads_pd(df, 
                           column='sequence' )
    logging.debug(f'calc/set read counts on original reads.')    
    df = set_counts_df(df, column='sequence')
    logging.info(f'dropping sequence column to slim.')
    df.drop('sequence', axis=1, inplace=True)    
    logging.info(f'Got dataframe len={len(df)} Writing to {outfile}')
    logging.debug(f'dataframe dtypes:\n{df.dtypes}\n')
    df.to_csv(outfile, sep='\t')
    dir, base, ext = split_path(outfile)
    outfile = os.path.join( dir , f'{base}.parquet')
    logging.info(f'df len={len(df)} as parquet to {outfile}...')
    df.to_parquet(outfile)
    logging.info(f'process_fastq_pairs done.')
    
    # align_collapse    
    outfile = f'{outdir}/{expid}.collapse.tsv'
    df = align_collapse_pd(df, 
                           outdir=outdir, 
                           cp=cp)
    logging.info(f'Saving len={len(df)} as TSV to {outfile}...')
    df.to_csv(outfile, sep='\t')
    dir, base, ext = split_path(outfile)
    outfile = os.path.join(dir, f'{base}.parquet')
    logging.info(f'df len={len(df)} as parquet to {outfile}...')
    df.to_parquet(outfile)    
    logging.info(f'align_collapse done.')
    
    # make_readtable
    outfile = f'{outdir}/{expid}.readtable.tsv'
    logging.debug(f'loading sample DF...')
    sampdf = load_sample_info(cp, sampleinfo)
    logging.debug(f'\n{sampdf}')
    sampdf.to_csv(f'{outdir}/sampleinfo.tsv', sep='\t')
    
    df = process_make_readtable_pd(df,
                                   sampdf,
                                   bcfile=bcfile, 
                                   outdir=outdir, 
                                   cp=cp)

    logging.info(f'Got dataframe len={len(df)} Writing to {outfile}')
    df.to_csv(outfile, sep='\t')
    dir, base, ext = split_path(outfile)
    outfile = os.path.join(dir, f'{base}.parquet')
    logging.info(f'df len={len(df)} as parquet to {outfile}...')
    df.to_parquet(outfile)
    logging.info(f'make_readtable done.')    

    # make_vbctable
    outfile = f'{outdir}/{expid}.vbctable.tsv'
    df = process_make_vbctable_pd(df,
                               outdir=outdir,
                               cp=cp)

    logging.debug(f'inbound df len={len(df)} columns={df.columns}')
    logging.info(f'Got dataframe len={len(df)} Writing to {outfile}')
    df.to_csv(outfile, sep='\t')
    
    dir, base, ext = split_path(outfile)
    outfile = os.path.join( dir , f'{base}.parquet')
    logging.info(f'df len={len(df)} as parquet to {outfile}...')
    df.to_parquet(outfile)    
    
    # make_matrices
    process_make_matrices_pd(df,
                               exp_id = expid,  
                               outdir=outdir, 
                               cp=cp)
    
    


def build_seqmapdict(udf, components, column='vbc_read'):
    '''
    Create mappings from all unique sequences to component sequence
    '''
    seqmapdict = {}
    comphandled = 0
    comphandled_interval = 100000
    comp_len = len(components)
    
    for comp in components:
        ser = list(udf[column].iloc[comp])
        t = ser[0]
        for s in ser: 
            seqmapdict[s]= t    
        if comphandled % comphandled_interval == 0:
            logging.debug(f'[{comphandled}/{comp_len}] t = {t}')
        comphandled += 1
    return seqmapdict


def build_seqmapdict_pd(udf, components, column='vbc_read', pcolumn='count'):
    '''
    Create mappings from all unique sequences to component sequence
    Can we do this faster? Choose most abundant variant?
    dict should be oldsequence -> newsequence
    
    '''
    seqmapdict = {}
    comphandled = 0
    pcolumn='count'
    logging.debug(f'udf len={len(udf)} components len={len(components)} column={column} pcolumn={pcolumn} ')
    comphandled_interval = 1000
    comp_len = len(components)    
    for i, indexlist in enumerate( components):
        cdf = udf[[column, pcolumn]].iloc[indexlist]
        cdf.reset_index(inplace=True, drop=True)
        #logging.debug(f'component [{i}/{comp_len}]: len={len(cdf)}')
        maxid = int(cdf[pcolumn].idxmax())
        t = cdf[column].iloc[maxid]
        for compseq in list(cdf[column]): 
            #logging.debug(f'compseq={compseq} -> {t}')
            seqmapdict[compseq]= t    
        if comphandled % comphandled_interval == 0:
            logging.debug(f'[{i}/{comp_len}]: len={len(cdf)} seq = {t} ')
        comphandled += 1
    return seqmapdict


def collapse_by_components_pd(fulldf, uniqdf, components, column, pcolumn, outdir=None):
    #
    # *** Assumes components are multi-element components only ***
    #
    # components consist of indices within the uniqdf. 
    # assumes component elements are integers, as they are used as dataframe indices. 
    # create map of indices in original full DF that correspond to all members of a (multi-)component
    # from alignment
    # hash key is index of first appearance of sequence in full DF (i.e. its index in uniqdf.) 
    # 
    # returns copy of input fulldf with sequence column collapsed to most-common sequence in component. 
    #
    #  SETUP
    #  infile is all reads...
    #      fdf = read_fasta_to_df(infile, seq_length=32)
    #      udf = pd.DataFrame(fdf['sequence'].unique(), columns=['sequence'])
    #      components = remove_singletons(components)
    
    logging.debug(f'multi-element components len={len(components)}')
    logging.debug(f'fulldf length={len(fulldf)} uniqdf length={len(uniqdf)} {len(components)} components.')
    logging.info(f'building seqmapdict {len(uniqdf)} unique seqs, {len(components)} components, for {len(fulldf)} raw sequences. ')
    smd = build_seqmapdict_pd(uniqdf, components,  column, pcolumn)
      
    # Make new full df:
    logging.info('seqmapdict built. Applying.')
    if outdir is not None:
        outfile = f'{outdir}/seqmapdict.json'
        logging.debug(f'writing seqmapdict len={len(smd)} tp {outfile}')
        with open(outfile, 'w') as fp:
            json.dump(smd, fp, indent=4)
    else:
        logging.debug(f'no outdir given.')
    logging.info(f'applying seqmapdict...')
    # make deep copy of original sequence column
    fulldf.loc[:, f'{column}_col'] = fulldf.loc[:, column]    
    # map old to new
    fulldf[f'{column}_col'] = fulldf[column].map(smd, na_action='ignore')
    # fill in NaN with original values. 
    # /Users/hover/git/mapseq-processing/mapseq/core.py:1930: FutureWarning: A value is trying to 
    # be set on a copy of a DataFrame or Series through chained assignment using an inplace method.
    # The behavior will change in pandas 3.0. This inplace method will never work because the 
    # intermediate object on which we are setting values always behaves as a copy.
    # For example, when doing 'df[col].method(value, inplace=True)', try using 'df.method({col: value}, inplace=True)' 
    #  or df[col] = df[col].method(value) instead, to perform the operation inplace on the original object.
    # fulldf['vbc_read_col'].fillna(fulldf['vbc_read'], inplace=True)
    
    # initial method:
    #  fulldf['vbc_read_col'].fillna(fulldf['vbc_read'], inplace=True)
    
    fulldf.fillna( { 'vbc_read_col' : fulldf['vbc_read'] }, inplace=True)    
    logging.info(f'New collapsed df = \n{fulldf}')
    log_objectinfo(fulldf, 'fulldf')
    return fulldf



def align_collapse_pd(df,
                      column='vbc_read',
                      pcolumn='read_count', 
                      max_mismatch=None, 
                      outdir=None, 
                      datestr=None, 
                      cp=None):
    '''
    Assumes dataframe with sequence and read_count columns
    Use read_count to choose parent sequence.
    Speed up use of pandas, map() function rather than apply()
    
    relationship between read_count on full_read sequences (52nt), and 
        value_counts() on unique VBCs (30nts)
    
    '''
    # housekeeping...
    if cp is None:
        cp = get_default_config()
    
    aligner = cp.get('collapse','tool')    
    if max_mismatch is None:
        max_mismatch = int(cp.get('collapse', 'max_mismatch'))
    else:
        max_mismatch = int(max_mismatch)
    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    if outdir is None:
        outdir = './'
    outdir = os.path.abspath(outdir)    
    os.makedirs(outdir, exist_ok=True)
    logging.debug(f'collapse: aligner={aligner} max_mismatch={max_mismatch} outdir={outdir}')    
    
    sh = get_default_stats()      
    sh.add_value('/collapse','n_full_sequences', len(df) )

    # get reduced dataframe of unique head sequences
    logging.info('Getting unique DF...')    
    #udf = pd.DataFrame(df['sequence'].unique(), columns=['sequence'])
    udf = df[column].value_counts().reset_index() 
    sh.add_value('/collapse','n_unique_sequences', len(udf) )    

    of = os.path.join( outdir , f'{column}.unique.tsv')
    logging.info(f'Writing unique DF to {of}')
    udf.to_csv(of, sep='\t') 
    
    of = os.path.join( outdir , f'{column}.unique.fasta')
    logging.info(f'Writing uniques as FASTA to {of}')
    seqfasta = write_fasta_from_df(udf, outfile=of, sequence=[column])    

    of = os.path.join(outdir, f'{column}.fulldf.tsv')
    logging.info(f'Writing slimmed full DF to {of}')    
    df.to_csv(of, sep='\t', columns=[column, pcolumn])

    # run allXall bowtiex
    of = os.path.join( outdir , f'unique_sequences.bt2.sam')
    logging.info(f'Running {aligner} on {seqfasta} file to {of}')
    btfile = run_bowtie(cp, seqfasta, of, tool=aligner)
    logging.info(f'Bowtie done. Produced output {btfile}. Creating btdf dataframe...')
    btdf = make_bowtie_df(btfile, max_mismatch=max_mismatch, ignore_self=True)
    of = os.path.join( outdir , f'unique_sequences.btdf.tsv')
    btdf.to_csv(of, sep='\t') 
    sh.add_value('/collapse','n_bowtie_entries', len(btdf) )

    # perform collapse...      
    logging.info('Calculating Hamming components...')
    edgelist = edges_from_btdf(btdf)
    btdf = None  # help memory usage
    sh.add_value('/collapse','n_edges', len(edgelist) )
    
    of = os.path.join( outdir , f'edgelist.txt')
    writelist(of, edgelist)
    logging.debug(f'edgelist len={len(edgelist)}')
    
    components = get_components(edgelist)
    logging.debug(f'all components len={len(components)}')
    sh.add_value('/collapse','n_components', len(components) )
    of = os.path.join( outdir , f'components.txt')
    writelist(of, components)
    edgelist = None  # help memory usage    

    # assess components...
    # components is list of lists.
    data = [ len(c) for c in components]
    data.sort(reverse=True)
    ccount = pd.Series(data)
    of = os.path.join( outdir , f'component_count.tsv')
    ccount.to_csv(of, sep='\t')            

    mcomponents = remove_singletons(components)
    logging.debug(f'multi-element components len={len(mcomponents)}')
    sh.add_value('/collapse','n_multi_components', len(mcomponents) )
    of = os.path.join( outdir , f'multi_components.json')
    logging.debug(f'writing components len={len(components)} t {of}')
    with open(of, 'w') as fp:
        json.dump(mcomponents, fp, indent=4)

    logging.info(f'Collapsing {len(components)} components...')
    newdf = collapse_by_components_pd(df, udf, mcomponents, column=column, pcolumn=pcolumn, outdir=outdir)
    # newdf has sequence and newsequence columns, rename to orig_seq and sequence
    newcol = f'{column}_col'        
    newdf.rename( {'newsequence': newcol}, inplace=True, axis=1)    
    logging.info(f'Got collapsed DF. len={len(newdf)}')

    rdf = newdf[[ column, newcol, 'read_count' ]]
    of = os.path.join( outdir , f'read_collapsed.tsv')
    logging.info(f'Writing reduced mapping TSV to {of}')
    rdf.to_csv(of, sep='\t')
    
    newdf.drop(column, inplace=True, axis=1)
    #joindf = pd.DataFrame( newdf['new_seq'] + newdf['tail'], columns=['sequence'])
    #of = os.path.join( outdir , f'collapsed.fasta')
    #logging.info(f'Writing collapsed fasta to {of}')
    #write_fasta_from_df(joindf, of)        
    return newdf        

       

def collapse_only_pd(fdf, udf, mcomponents, column, pcolumn, outdir, cp ):
    '''
    Handle partially completed align-collapse. 
    Inputs:
       fulldf TSV
       uniqudf TSV
       multi-components list. 
    
    Output:
        collapsed DF. 
        df = collapse_pd(fulldf,
                     uniquedf,
                     mcomponents,   
                     column=args.column,
                     pcolumn=args.parent_column,
                     outdir=outdir, 
                     cp=cp)
    '''
  
    # assess components...
    sh = get_default_stats()
    logging.debug(f'multi-element components len={len(mcomponents)}')
    sh.add_value('/collapse','n_multi_components', len(mcomponents) )

    logging.info(f'Collapsing {len(mcomponents)} mcomponents...')
    newdf = collapse_by_components_pd(fdf, udf, mcomponents, column=column, pcolumn=pcolumn, outdir=outdir)

    # newdf has sequence and newsequence columns, rename to orig_seq and sequence
    newcol = f'{column}_col'        
    newdf.rename( {'newsequence': newcol}, inplace=True, axis=1)    
    logging.info(f'Got collapsed DF. len={len(newdf)}')

    rdf = newdf[[ column, newcol, 'read_count' ]]
    of = os.path.join( outdir , f'read_collapsed.tsv')
    logging.info(f'Writing reduced mapping TSV to {of}')
    rdf.to_csv(of, sep='\t')
    logging.debug(f'dropping {column}')
    newdf.drop(column, inplace=True, axis=1)   
    return newdf        
