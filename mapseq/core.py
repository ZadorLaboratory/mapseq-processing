import copy
import glob
import gzip
import itertools
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
from natsort import natsorted, index_natsorted, order_by_index

import scipy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

import dask
import dask.dataframe as dd
from dask.dataframe import from_pandas
from dask.distributed import Client

from mapseq.utils import *
from mapseq.bowtie import *
from mapseq.barcode import *
from mapseq.stats import *
from mapseq.plotting import *

#
# STANDARD FORMATS AND DATATYPES
#

FFORMATS = ['reads','aggregated','filtered','readtable','collapsed','vbctable']

STR_COLS = {
    'reads'      : ['sequence'],
    'aggregated' : ['sequence'],
    'filtered'   : ['vbc_read', 'spikeseq', 'libtag', 'umi',  'ssi'],
    'readtable'  : ['vbc_read', 'umi' ],
    #'collapsed'  : ['vbc_read','spikeseq', 'libtag', 'umi',  'ssi'],
    'collapsed'  : ['vbc_read', 'umi'],
    'vbctable'   : ['vbc_read'],      
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
    'reads'      : ['source'],
    'aggregated' : ['source'],
    'filtered'   : ['source'],
    'collapsed'   : ['source'],
    'readtable'  : ['label','site','type','brain','region','source','ourtube','rtprimer'],
    'vbctable'   : ['label','site','type','brain','region','ourtube'],        
    }

FMT_DTYPES = {      'read_count'    : 'int64',
                    'umi_count'     : 'int64',
                    'sequence'      : 'string[pyarrow]',
                    'vbc_read'      : 'string[pyarrow]',
                    'vbc_read_col'  : 'string[pyarrow]',        
                    'spikeseq'      : 'string[pyarrow]',    
                    'libtag'        : 'string[pyarrow]',    
                    'umi'           : 'string[pyarrow]',    
                    'ssi'           : 'string[pyarrow]',
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
#
# TARGET
# target-negative            user-defined    treated but expected to be low

# CONTROLS
# target-negative-control    user-provided   untreated sample
# target-wt-control          core added      untreated biological sample
# target-water-control       core added      empty sample

# injection-water-control    core added      empty sample

CONTROL_SITES=['target-negative-control', 
               'target-wt-control',
               'target-water-control',
               'injection-water-control']

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
            df = fix_mapseq_df_types(df, fformat=fformat, use_arrow=use_arrow) 
    
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



def fix_mapseq_df_types(df, fformat='reads', use_arrow=True):
    '''
    confirm that columns are proper types. 
    use string[pyarrow] for strings, else string
    '''
    logging.info(f'old dataframe dtypes=\n{df.dtypes}')
    
    if fformat in FFORMATS:
        if use_arrow:
            tt = 'string[pyarrow]'
        else:
            tt = 'string[python]'
            
        for scol in STR_COLS[fformat]:
            dt = df[scol].dtype 
            if dt != tt:
                logging.debug(f'converting col={scol} from {dt} to {tt} ...')
                df[scol] = df[scol].astype(tt) 
            
        for icol in INT_COLS[fformat]:
            dt = df[icol].dtype
            if dt != 'int64':
                logging.debug(f"converting col={icol} from {dt} to 'int64' ...")
                df[icol] = df[icol].astype('int64')                
            
        for ccol in CAT_COLS[fformat]:
            dt = df[ccol].dtype
            if dt != 'category':
                logging.debug(f"converting col={ccol} from {dt} to 'category' ...")
                df[ccol] = df[ccol].astype('category')        
    else:
        logging.warning('unrecognized mapseq format. return original')
    logging.info(f'new dataframe dtypes=\n{df.dtypes}')
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
    
    sample_columns = ['usertube', 'ourtube', 'samplename', 'siteinfo', 'si_ratio', 'min_reads', 'rtprimer', 'brain', 'region', 'matrixcolumn'] 

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
        logging.info(f'loaded DF from Excel {file_name}')
        
    elif file_name.endswith('.tsv'):
        sdf = pd.read_csv(file_name, sep='\t', index_col=0, keep_default_na=False, dtype =str, comment="#")
        sdf = sdf.astype('str', copy=False)    
    else:
        logging.error(f'file {file_name} neither .xlsx or .tsv')
        sdf = None
        
    logging.debug(f'created reduced sample info df:\n{sdf}')
    return sdf


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


  
#
#  BOWTIE/TARJAN FUNCTIONS
#
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


####################################################################
#    PIPELINE PHASES
####################################################################

#
#    FASTQ  HANDLING 
#
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
    vcs = df.value_counts( column )
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
    Fuction to be handed to Dask client for scalable calculation. 
    
    '''
    return x['sequence'].value_counts()


def aggregate_reads_dd_clientXXX(df, 
                              column='sequence', 
                              outdir=None, 
                              min_reads=1, 
                              chunksize=50000000, 
                              dask_temp='./temp'):
    '''
    
    WARNING dask 
    
    ASSUMES INPUT IS DASK DATAFRAME
    retain other columns and keep first value

    Suggestions for partitions. 
    https://stackoverflow.com/questions/44657631/strategy-for-partitioning-dask-dataframes-efficiently

    Limiting resource usage (CPU, memory)
    https://stackoverflow.com/questions/59866151/limit-dask-cpu-and-memory-usage-single-node
    
    '''
    if dask_temp is not None:
        logging.info(f'setting dask temp to {os.path.expanduser(dask_temp)} ')
        dask.config.set(temporary_directory=os.path.abspath( os.path.expanduser(dask_temp)))
    else:
        logging.info(f'no dask_temp specified. letting dask use its default')


    client = Client( memory_limit='16GB', 
                     processes=True,
                     n_workers=8, 
                     threads_per_worker=2)
       
    before_len = len(df)
    logging.debug(f'collapsing with read counts for col={column} len={before_len}')
    sent = client.submit(sequence_value_counts, df)
    ndf = sent.result().compute()
    client.close()
    
    ndf.reset_index(inplace=True, drop=True)
    ndf.rename({'count':'read_count'}, inplace=True, axis=1)
    logging.info(f'computed counts. new DF len={len(ndf)}')
    
    # get back all other columns from original set. e.g. 'source'
    logging.info(f'merging to recover other columns from original DF')
    ndf = dd.from_pandas(ndf)
    #ndf = pd.merge(ndf, seqdf.drop_duplicates(subset=column,keep='first'),on=column, how='left')  
    result = ndf.merge(seqdf.drop_duplicates(subset=column,keep='first'), on=column, how='left')  
    ndf = result.compute()
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


#
#    SPLIT MAPSEQ COLUMNS, FILTER             
#

def split_mapseq_fields(df, column='sequence', drop=False, cp=None):
    '''
    Used by filter_split
    
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
        drop_mismatch = cp.getboolean('readfilter','drop_mismatch')

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
    @arg cp             ConfigParser object with [collapse] section.
    
    '''
    if gcolumn is not None:
        logging.info(f'Grouped align_collapse. Group column = {gcolumn}')
        df = align_collapse_pd_grouped( df, 
                                        column = column,
                                        pcolumn = pcolumn,
                                        gcolumn = gcolumn,
                                        max_mismatch=max_recursion,
                                        max_recursion=max_recursion, 
                                        outdir=outdir, 
                                        datestr=datestr,
                                        force=force, 
                                        cp=cp )
    else:
        logging.info(f'Global align_collapse.')
        df = align_collapse_pd( df = df,
                                column = column,
                                pcolumn = pcolumn,
                                max_mismatch=max_recursion,
                                max_recursion=max_recursion, 
                                outdir=outdir, 
                                datestr=datestr,
                                force=force, 
                                cp=cp ) 
    logging.info(f'Done. returning DF len={len(df)}')
    return df
                            

def align_collapse_pd(df,
                      column='vbc_read',
                      pcolumn='read_count', 
                      max_mismatch=None,
                      max_recursion=None, 
                      outdir=None, 
                      datestr=None,
                      force=False, 
                      cp=None):
    '''
    Assumes dataframe with sequence and read_count columns
    Use read_count to choose parent sequence.
    max_recursion appears to be important, rather than memory reqs. 40000 or 50000 
        novaseq (1.5B reads) may be needed. 
    By default, picking up wherever the algorithm left off, depending on
        file availability. Force overrides and recalculates from beginning. 
      
    TODO:
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

    if max_recursion is not None:
        rlimit = int(max_recursion)
        logging.info(f'set new recursionlimit={rlimit}')
        sys.setrecursionlimit(rlimit)
    else:
        rlimit = int(cp.get('collapse','max_recursion'))
        logging.info(f'set new recursionlimit={rlimit}')
        sys.setrecursionlimit(rlimit)        

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

    #of = os.path.join( outdir , f'{column}.unique.tsv')
    #logging.info(f'Writing unique DF to {of}')
    #udf.to_csv(of, sep='\t') 
    
    of = os.path.join( outdir , f'{column}.unique.fasta')
    logging.info(f'Writing uniques as FASTA to {of}')
    seqfasta = write_fasta_from_df(udf, outfile=of, sequence=[column])    

    #of = os.path.join(outdir, f'{column}.fulldf.tsv')
    #logging.info(f'Writing slimmed full DF to {of}')    
    #df.to_csv(of, sep='\t', columns=[column, pcolumn])

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
    newdf.rename( { newcol : column }, inplace=True, axis=1) 
    #joindf = pd.DataFrame( newdf['new_seq'] + newdf['tail'], columns=['sequence'])
    #of = os.path.join( outdir , f'collapsed.fasta')
    #logging.info(f'Writing collapsed fasta to {of}')
    #write_fasta_from_df(joindf, of)        
    return newdf        


def align_collapse_pd_grouped(df,
                              column='vbc_read',
                              pcolumn='read_count',
                              gcolumn='brain', 
                              max_mismatch=None,
                              max_recursion=None, 
                              outdir=None, 
                              datestr=None,
                              force=False, 
                              cp=None):
    '''
    Groups alignment and collapse by gcolumn value. [brain]
    Uses sub-directories for standard intermediate output/scratch. 
    Assumes dataframe with sequence (vbc_read) and read_count columns
    Use read_count to choose parent sequence.
        
    max_recursion appears to be important, rather than memory reqs. 40000 or 50000 
        novaseq (1.5B reads) may be needed. 
    By default, picking up wherever the algorithm left off, depending on
        file availability. Force overrides and recalculates from beginning. 
      
    TODO:
    Speed up use of pandas, map() function rather than apply()
    relationship between read_count on full_read sequences (52nt), and 
        value_counts() on unique VBCs (30nts)
    
    '''
    # housekeeping...
    if cp is None:
        cp = get_default_config()
        
    if max_mismatch is None:
        max_mismatch = int(cp.get('collapse', 'max_mismatch'))
    else:
        max_mismatch = int(max_mismatch)

    if max_recursion is not None:
        rlimit = int(max_recursion)
        logging.info(f'set new recursionlimit={rlimit}')
        sys.setrecursionlimit(rlimit)
    else:
        rlimit = int(cp.get('collapse','max_recursion'))
        logging.info(f'(from config) set new recursionlimit={rlimit}')
        sys.setrecursionlimit(rlimit)        

    aligner = cp.get('collapse','tool')

    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    
    if outdir is None:
        outdir = './'
    outdir = os.path.abspath(outdir)    
    os.makedirs(outdir, exist_ok=True)

    logging.debug(f'collapse: aligner={aligner} max_mismatch={max_mismatch} outdir={outdir}')    
    sh = get_default_stats()      
    sh.add_value('/collapse','n_full_sequences', len(df) )

    # Get list of ids to group collapse by...
    df[gcolumn] = df[gcolumn].astype('string')
    gidlist = list( df[gcolumn].dropna().unique() )
    gidlist = [ x for x in gidlist if len(x) > 0 ]
    gidlist.sort()
    logging.debug(f'handling group list: {gidlist}')

    logging.info(f'pulling out no group for merging at end...')
    nogroup_df = df[ df[gcolumn] == '' ]
    sh.add_value('/collapse','n_no_group', len(nogroup_df) )

    gdflist = []

    for gid in gidlist:
        logging.info(f"collapsing '{column}' by {gcolumn} = '{gid}' ")
        gdir = os.path.join( outdir, f'{gcolumn}.{gid}' )
        os.makedirs(gdir, exist_ok=True)
        
        gdf = df[df[gcolumn] == gid]
        gdf.reset_index(inplace=True, drop=True)
        initial_len = len(gdf)  
        logging.info(f'[{gcolumn}:{gid}] initial len={len(gdf)} subdir={gdir}')        
        
        # get reduced dataframe of unique head sequences
        logging.info('Getting unique DF...')    
        udf = gdf[column].value_counts().reset_index() 
        sh.add_value(f'/collapse/{gcolumn}_{gid}','n_unique_sequences', len(udf) )    
    
        of = os.path.join( gdir , f'{column}.unique.tsv')
        logging.info(f'Writing unique DF to {of}')
        udf.to_csv(of, sep='\t') 
        
        of = os.path.join( gdir , f'{column}.unique.fasta')
        logging.info(f'Writing uniques as FASTA to {of}')
        seqfasta = write_fasta_from_df(udf, outfile=of, sequence=[column])    
    
        #of = os.path.join(outdir, f'{column}.fulldf.tsv')
        #logging.info(f'Writing slimmed full DF to {of}')    
        #df.to_csv(of, sep='\t', columns=[column, pcolumn])
    
        # run allXall bowtiex
        of = os.path.join( gdir , f'unique_sequences.bt2.sam')
        logging.info(f'Running {aligner} on {seqfasta} file to {of}')
        btfile = run_bowtie(cp, seqfasta, of, tool=aligner)
        logging.info(f'Bowtie done. Produced output {btfile}. Creating btdf dataframe...')
        btdf = make_bowtie_df(btfile, max_mismatch=max_mismatch, ignore_self=True)
        of = os.path.join( gdir , f'unique_sequences.btdf.tsv')
        btdf.to_csv(of, sep='\t') 
        sh.add_value(f'/collapse/{gcolumn}_{gid}','n_bowtie_entries', len(btdf) )
    
        # perform collapse...      
        logging.info('Calculating Hamming components...')
        edgelist = edges_from_btdf(btdf)
        btdf = None  # help memory usage
        sh.add_value(f'/collapse/{gcolumn}_{gid}','n_edges', len(edgelist) )
        
        of = os.path.join( gdir , f'edgelist.txt')
        writelist(of, edgelist)
        logging.debug(f'edgelist len={len(edgelist)}')
        sh.add_value(f'/collapse/{gcolumn}_{gid}','n_edges', len(edgelist) )
        
        components = get_components(edgelist)
        logging.debug(f'all components len={len(components)}')
        sh.add_value(f'/collapse/{gcolumn}_{gid}','n_components', len(components) )
        of = os.path.join( gdir , f'components.txt')
        writelist(of, components)
        edgelist = None  # help memory usage    
    
        # assess components...
        # components is list of lists.
        data = [ len(c) for c in components]
        data.sort(reverse=True)
        ccount = pd.Series(data)
        of = os.path.join( gdir , f'component_count.tsv')
        ccount.to_csv(of, sep='\t')            
    
        mcomponents = remove_singletons(components)
        logging.debug(f'multi-element components len={len(mcomponents)}')
        sh.add_value(f'/collapse/{gcolumn}_{gid}','n_multi_components', len(mcomponents) )
        of = os.path.join( gdir , f'multi_components.json')
        logging.debug(f'writing components len={len(components)} t {of}')
        with open(of, 'w') as fp:
            json.dump(mcomponents, fp, indent=4)
    
        logging.info(f'Collapsing {len(components)} components...')
        newdf = collapse_by_components_pd(gdf, udf, mcomponents, column=column, pcolumn=pcolumn, outdir=gdir)
        # newdf has sequence and newsequence columns, rename to orig_seq and sequence
        newcol = f'{column}_col'        
        newdf.rename( {'newsequence': newcol}, inplace=True, axis=1)    
        logging.info(f'Got collapsed DF. len={len(newdf)}')
    
        rdf = newdf[[ column, newcol, 'read_count' ]]
        of = os.path.join( gdir , f'read_collapsed.tsv')
        logging.info(f'Writing reduced mapping TSV to {of}')
        rdf.to_csv(of, sep='\t')
        
        newdf.drop(column, inplace=True, axis=1)
        newdf.rename( { newcol : column }, inplace=True, axis=1) 
        gdflist.append(newdf)
        #joindf = pd.DataFrame( newdf['new_seq'] + newdf['tail'], columns=['sequence'])
        #of = os.path.join( outdir , f'collapsed.fasta')
        #logging.info(f'Writing collapsed fasta to {of}')
        #write_fasta_from_df(joindf, of)        
    
    # merge all brains into one dataframe...
    logging.debug(f'sizes={[ len(x) for x in gdflist ]} adding nogroup len={len(nogroup_df)}')    
    gdflist.append(nogroup_df)
    outdf = pd.concat(gdflist, ignore_index = True)
    outdf.reset_index(inplace=True, drop=True)
    logging.info(f'All groups. Final DF len={len(outdf)}')
    outdf = fix_mapseq_df_types(outdf, fformat='readtable')
    return outdf
       

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
    
    fulldf.fillna( { 'vbc_read_col' : fulldf['vbc_read'] }, inplace=True)    
    logging.info(f'New collapsed df = \n{fulldf}')
    log_objectinfo(fulldf, 'fulldf')
    return fulldf

       

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
    n_initial = len(df)
    if outdir is None:
        outdir = os.path.abspath('./')
    outdir = os.path.abspath(outdir)    
    os.makedirs(outdir, exist_ok=True)
    
    if cp is None:
        cp = get_default_config()
    
    spikeseq = cp.get('readtable','spikeseq')
    realregex = cp.get('readtable', 'realregex' )
    loneregex = cp.get('readtable', 'loneregex' )
    use_libtag = cp.getboolean('readtable','use_libtag')
    
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
    #logging.info('filling in label by SSI sequence...')
    #df['label'] = df['ssi'].map(bcdict)

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

    badsitedf = df[ df.isna().any(axis=1) ]
    n_badsite = len(badsitedf)
    df.drop(badsitedf.index, inplace=True)
    df.reset_index(inplace=True, drop=True)    

    of = os.path.join( outdir, 'bad_site.tsv')
    badsitedf.reset_index(inplace=True, drop=True)
    start = dt.datetime.now()
    badsitedf.to_csv(of, sep='\t')
    end = dt.datetime.now()
    delta_seconds = (dt.datetime.now() - start).seconds
    log_transferinfo(of, delta_seconds)
    
    logging.info(f'Wrote bad_site DF len={len(badsitedf)} to {of}')
    badsitedf = None   

    if use_libtag:
        # L1/L2 libtag
        logging.info('identifying L2 (reals) by libtag...')
        rmap = df['libtag'].str.match(realregex)
        df.loc[rmap, 'type'] = 'real'
        
        logging.info('identifying L1s by libtag...')
        lmap = df['libtag'].str.match(loneregex)
        df.loc[lmap, 'type'] = 'lone'
        
        logging.info(f'Identifying spikeins by spikeseq={spikeseq}')
        smap = df['spikeseq'] == spikeseq
        df.loc[smap, 'type'] = 'spike'
        
    else:
        logging.info(f'Identifying spikeins by spikeseq={spikeseq}')
        smap = df['spikeseq'] == spikeseq
        df.loc[smap, 'type'] = 'spike'
        
        logging.info('ignoring L1/L2 libtag. All non-spikes are real.')        
        nsmap = df['type'] != 'spike'
        df.loc[nsmap, 'type'] = 'real'

    # identify bad type rows.
    # must not be spikein, and libtag must not match L1 or L2 
    # i.e. neither all purines or all pyrimidenes
    logging.debug('Identifying bad type rows.')  
    badtypedf = df[ df.isna().any(axis=1) ]
    n_badtype = len(badtypedf)
    df.drop(badtypedf.index, inplace=True)
    df.reset_index(inplace=True, drop=True)    

    of = os.path.join( outdir, 'bad_type.tsv')
    badtypedf.reset_index(inplace=True, drop=True)
    badtypedf.to_csv(of, sep='\t')
    logging.info(f'Wrote bad_type DF len={len(badtypedf)} to {of}')
    badtypedf = None   
    
    logging.debug('Dropping redundant sequence fields (spikeseq, libtag).')
    df.drop(['spikeseq','libtag'], inplace=True, axis=1)

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
    
    # calc/ collect stats
    sh = get_default_stats()    
    sh.add_value('/readtable','n_full_sequences', str(len(df)) )
    
    # find and remove (at least) known template-switch rows from dataframe.
    # template switch type is L1 (from libtab) but is a valid target (from SSI) 
    tsdf = df[ ((df['type'] == 'lone') & ( df['site'].str.startswith('target'))) ]
    of = os.path.join(outdir, 'template_switch.tsv') 
    logging.info(f'Writing template switch DF len={len(tsdf)} Writing to {of}')
    tsdf.to_csv(of, sep='\t')
    n_tswitch = len(tsdf)
    
    # remove known template switch from readtable
    df.drop(tsdf.index, inplace=True)
    df.reset_index(drop=True, inplace=True)
    n_final = len(df)
    
    sh.add_value('/readtable', 'n_initial', str(n_initial) ) 
    sh.add_value('/readtable', 'n_badssi', str(n_badssi) )
    sh.add_value('/readtable', 'n_badtype', str(n_badtype) )
    sh.add_value('/readtable', 'n_badsite', str(n_badsite) )
    sh.add_value('/readtable', 'n_tswitch', str(n_tswitch) )
    sh.add_value('/readtable', 'n_final', str(n_final) )     
    
    df = fix_mapseq_df_types(df, fformat='readtable')
    return df

#
#            VBCTABLE
#

def process_make_vbctable_pd(df,
                          outdir=None,
                          inj_min_reads = 2,
                          target_min_reads = 2, 
                          cp = None):
    '''   
    
    primary columns:  label, umi, type  
    derived columns:  brain, region, rtprimer, label, ssi
    -- remove nomatches
    -- threshold read_count by site type
    -- collapse by VBC sequence, calculating UMI count
    
    Drops source field since aggregation would make it ambiguous.
    
    This should  be prepared for input it make VBC matrices    
    '''
    if cp is None:
        cp = get_default_config()
    project_id = cp.get('project','project_id')
    
    logging.debug(f'inbound df len={len(df)} columns={df.columns}')
    log_objectinfo(df, 'input-df')
    ndf = df
    
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
    udf = thdf.groupby(['vbc_read','label','type'], observed=True ).agg( {'umi' : 'nunique',
                                                                              'read_count':'sum', 
                                                                              'brain':'first',
                                                                              'region':'first',
                                                                              'site':'first',
                                                                              'ourtube':'first',    
                                                                              } ).reset_index()
    udf.rename( {'umi':'umi_count'}, axis=1, inplace=True)
    logging.info(f'DF after umi/label collapse: {len(udf)}')

    # remove internal controls that end user needn't see
    # remove L1s
    # pull out L1 VBCs expected to have L1 libtag. 
    lones = udf[ (udf['site'] == 'target-lone') & (udf['type'] == 'lone') ]
    udf = udf[ udf['site'] != 'target-lone' ]

    
    # pull out remaining VBCs with matching SSIs that we did not expect to have L1 libtags. 
    anomalies = udf[ udf['type'] == 'lone' ]
    udf = udf[ udf['type'] != 'lone' ]
    anomalies.reset_index(inplace=True, drop=True)  
    anomalies.to_csv(f'{outdir}/{project_id}.anomalies.tsv', sep='\t')

    # output L1s and anomalies. 
    lones.reset_index(inplace=True, drop=True)
    lones.to_csv(f'{outdir}/{project_id}.lones.tsv', sep='\t')

    # output controls by SSI/site, save to TSV
    controls = udf[ udf['site'].isin( CONTROL_SITES ) ]
    controls = controls[ controls['type'] == 'real' ]
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

    require_injection = cp.getboolean('vbcfilter','require_injection')
    include_injection = cp.getboolean('vbcfilter','include_injection')
    include_controls = cp.getboolean('vbcfilter','include_controls')
    use_target_negative=cp.getboolean('vbcfilter','use_target_negative')
    use_target_water_control=cp.getboolean('vbcfilter','use_target_water_control')

    df['brain'] = df['brain'].astype('string')
    bidlist = list(df['brain'].dropna().unique())
    bidlist = [ x for x in bidlist if len(x) > 0 ]
    bidlist.sort()
    logging.debug(f'handling brain list: {bidlist}')
    
    sh = get_default_stats()

    norm_dict = {}

    bdflist = []

    for brain_id in bidlist:
        valid = True
        ndf = None 
        logging.debug(f'handling brain_id={brain_id}')
        bdf = df[df['brain'] == brain_id]
        bdf.reset_index(inplace=True, drop=True)
        initial_len = len(bdf)  
        logging.info(f'[{brain_id}] initial len={len(bdf)}')
        max_negative = 1
        max_water_control = 1
        
        logging.debug(f'[{brain_id}]: inj_min_umi={inj_min_umi} target_min_umi={target_min_umi}')
        logging.debug(f'[{brain_id}]:use_target_negative={use_target_negative} use_target_water_control={use_target_water_control}')
        sh = get_default_stats() 
    
        # select all controls by SSI/site, save to TSV
        controls = bdf[ bdf['site'].isin( CONTROL_SITES ) ]
        controls.reset_index(inplace=True, drop=True)
        controls.to_csv(f'{outdir}/{brain_id}.removed_controls.tsv', sep='\t')
        
        # optionally keep/remove for inclusion in each brain matrix. 
        if not include_controls:
            bdf = bdf[ bdf['site'].isin( CONTROL_SITES ) == False]
    
        # remove spikes and save them 
        # Spikes to NOT get thresholded by UMI, or restricted by injection/target presence 
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
    df = pd.concat(bdflist, ignore_index = True)
    df.reset_index(inplace=True, drop=True)
    logging.info(f'All brains. Final merged filtered DF len={len(df)}')
    df = fix_mapseq_df_types(df, fformat='vbctable')
    return df


def filter_targets_min_umi_any(targets, min_umi ):
    '''
    retain all VBC rows for which any target label umi_count exceeds target_min_umi
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

#
#    MATRICES
#


def process_make_matrices(df,
                          outdir=None,
                          exp_id = None,  
                          label_column='label',
                          cp = None):
    '''
    
    Simplified matrix creation. 
    NOTE: Assumes ALL input is valid, thresholded, and to be included.
    
    -- per brain, pivot real VBCs (value=umi counts)
    -- create real, real normalized by spike-in, and scaled normalized.   
    -- use label_column to pivot on, making it the y-axis, x-axis is vbc sequence.  
    
    '''
    if cp is None:
        cp = get_default_config()
    if outdir is None:
        outdir = os.path.abspath('./')
    if exp_id is None:
        exp_id = cp.get('project','project_id')
   
    logging.debug(f'running exp_id={exp_id} ')

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
            nbcmdf = nbcmdf[ natsorted( list(nbcmdf.columns) ) ]            
            logging.debug(f'nbcmdf.describe()=\n{nbcmdf.describe()}')
            norm_dict[brain_id] = nbcmdf
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
    return norm_dict


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


def normalize_weight_grouped(df, weightdf, columns=None):
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
    creates excel sheet of unique VBC counts for 
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
    
    vdf = df.groupby(by=['label','type'],observed=False).agg( {'vbc_read':'nunique'} )
    vdf.reset_index(inplace=True, drop=False)
    vdf.sort_values(by='label', inplace=True, key=lambda x: np.argsort(index_natsorted( vdf['label'])))
    vdf.reset_index(inplace=True, drop=True)

    sdf = df[df['type'] == 'spike']
    sdf = sdf.groupby(by=['label','type'],observed=True).agg( {'vbc_read':'nunique'} )
    sdf.reset_index(inplace=True, drop=False)
    sdf.sort_values(by='label', inplace=True, key=lambda x: np.argsort(index_natsorted( sdf['label'])))
    sdf.reset_index(inplace=True, drop=True)    

    rdf = df[df['type'] == 'real']
    rdf = rdf.groupby(by=['label','type'],observed=True).agg( {'vbc_read':'nunique'} )
    rdf.reset_index(inplace=True, drop=False)
    rdf.sort_values(by='label', inplace=True, key=lambda x: np.argsort(index_natsorted( rdf['label'])))
    rdf.reset_index(inplace=True, drop=True)

    rcdf = df.groupby(by=['label','type'],observed=True).agg( {'read_count':'sum'} )
    rcdf.reset_index(inplace=True,drop=False)
    rcdf.sort_values(by='label', inplace=True, key=lambda x: np.argsort(index_natsorted( rcdf['label'])))
    rcdf.reset_index(inplace=True, drop=True)    

    with pd.ExcelWriter(outfile) as writer:
        vdf.to_excel(writer, sheet_name='Unique VBC')
        sdf.to_excel(writer, sheet_name='Spike Unique VBC')
        rdf.to_excel(writer, sheet_name='Real Unique VBC')
        rcdf.to_excel(writer, sheet_name='Read Count Sum')

    logging.info(f'Wrote XLSX report: {outfile} ')
    


def make_vbctable_qctables(df, 
                           outdir=None, 
                           cp=None, 
                           cols=['site','type'], 
                           vals = ['vbc_read','label','umi_count','read_count'], 
                           sort_by='umi_count' ):
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
    #sh.add_value('/vbctable','n_vbcs', len(udf) )    
        
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
    
    xlout = os.path.join(outdir, f'{project_id}.controls.xlsx')
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
                sh.add_value('/vbctable',f'n_{fname}_vbcs', len(ndf) )
                # If it is a real control, add to control report. 
                for s in CONTROL_SITES:
                    if ( ( s in fname) and ('real' in fname)):
                        logging.debug(f"writing control '{fname}' len={len(ndf)} to {xlout}")
                        ndf.to_excel(writer, sheet_name=fname)
                    else:
                        pass
                        #logging.debug(f"'{fname}' tuple not {s} and 'real'.")                    
            else:
                logging.info(f'no entries for {fname}')
        logging.debug(f'done with all tuples in {combinations}')
        
        
def make_vbctable_parameter_report_xlsx(df, 
                                   outdir=None, 
                                   cp=None, 
                                   params=None
                                   ):
    '''
        Create single XLSX with table of number of VBCs per SSI at the various thrsholds
        for injection, target UMI counts. 
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
        params = eval(  cp.get( 'vbctable','test_params') ) 
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
    testcp.set('vbcfilter','include_injection', 'True')
    
    for i, (inj_min_umi, target_min_umi) in enumerate(params):
        logging.debug(f'inj_min_umi = {inj_min_umi} target_min_umi = {target_min_umi} ')
        colname = f'inj:{inj_min_umi}, tar:{target_min_umi}'
        fdf = process_filter_vbctable(df, 
                                      inj_min_umi=inj_min_umi, 
                                      target_min_umi = target_min_umi, 
                                      target_min_umi_absolute=1, 
                                      outdir = outdir, 
                                      cp=testcp)
        fdf = fdf[fdf['type'] == 'real']
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
    logging.info(f'Wrote XLSX report: {outfile} ')

    return outdf
        
        
        
        

