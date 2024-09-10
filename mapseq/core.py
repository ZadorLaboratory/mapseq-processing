import glob
import gzip
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
    

def process_ssifasta_files(config, sampleinfo, infilelist, numthreads=1, outdir=None):
    '''
    Process each infile in separate process.
    Produce {outdir}/<BASE>.all.tsv 
    These are then combined to product outfile.   
    
    '''
    logging.debug(f'called: sampleinfo={sampleinfo} outdir={outdir}')
    if outdir is None:
        afile = infilelist[0]
        filepath = os.path.abspath(afile)    
        dirname = os.path.dirname(filepath)
        outdir = dirname
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
        logging.debug(f'made outdir={outdir}')

    ncpus, threads = calc_thread_count(numthreads)   
    logging.info(f'using {threads} of {ncpus} CPUs in parallel...')
    proglog = '-v'
    if logging.getLogger().level == logging.INFO:
        proglog = '-v'
    elif logging.getLogger().level == logging.DEBUG:
        proglog = '-d'
    
   
    cfilename =  f'{outdir}/process_ssifasta.config.txt'
    configfile = write_config(config, cfilename, timestamp=True)
    
    prog = os.path.expanduser('~/git/mapseq-processing/scripts/process_ssifasta_single.py')

    jstack = JobStack()
    outfilelist = []
    for infile in infilelist:
        base = get_mainbase(infile)
        logfile = f'{outdir}/{base}.log'
        outfile = f'{outdir}/{base}.all.tsv'
        outfilelist.append(outfile)
        cmd = [ prog, 
               proglog,
               '-c', configfile , 
               '-L', logfile ,
               '-s', sampleinfo ,  
               '-O' , outdir ,
               '-o' , outfile,
               infile
               ]
        jstack.addjob(cmd)
    jset = JobSet(max_processes = threads, jobstack = jstack)
    jset.runjobs()

    logging.info(f'finished jobs. merging output files:{outfilelist} ')    
    outdf = merge_tsvs(outfilelist)
    logging.info(f'made final DF len={len(outdf)}')
    return outdf


def process_ssifasta(config, infile, outdir=None, site=None, datestr=None):
    '''
    by default, outdir will be same dir as infile
    assumes infile fasta has already been trimmed to remove SSI
    
    site = ['target-control','injection-control','target','target-negative',target-lone']   
    Will use relevant threshold. If None, will use default threshold
    
    '''
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    
    if outdir is None:
        outdir = "."
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logging.debug(f'made outdir={outdir}')

    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    sh = StatsHandler(config, outdir=outdir, datestr=datestr)
    
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    logging.debug(f'handling {filepath} base={base}')
    
    # make raw fasta TSV of barcode-splitter output for one barcode. 
    # trim to 44 nt since we know last 8 are SSI  
    logging.debug('calc counts...')
    seqdf = make_fasta_df(config, infile)
    of = os.path.join(dirname , f'{base}.read.seq.tsv')
    seqdf.to_csv(of, sep='\t')
    
    # to calculate threshold we need counts calculated. 
    cdf = make_read_counts_df(config, seqdf, label=base)  
    logging.debug(f'initial counts df {len(cdf)} all reads.')
    
    # these are ***READ*** counts
    of = os.path.join(dirname , f'{base}.read.counts.tsv')
    cdf.to_csv(of, sep='\t') 
        
    threshold = get_read_count_threshold(config, cdf, site)
    logging.debug(f'got threshold={threshold} for site {site}')
    tdf = threshold_read_counts(config, cdf, threshold=threshold)
    logging.info(f'at threshold={threshold} {len(tdf)} unique molecules.')
    
    # thresholded raw counts. duplicates are all one UMI, so set counts to 1. 
    # each row (should be) a distinct UMI, so trim to 32. 
    #tdf['counts'] = 1          
    tdf['sequence'] = tdf['sequence'].str[:32]
    tdf['umi_count'] = 1
    
    # this contains duplicate VBCs with *different* UMIs
    of = os.path.join(dirname , f'{base}.umi.seq.tsv')
    tdf.to_csv(of, sep='\t') 
    
    # now we have actual viral barcode df with *unique molecule counts.*
    vbcdf = make_umi_counts_df(config, tdf)
    of = os.path.join(dirname , f'{base}.umi.counts.tsv')
    vbcdf.to_csv(of, sep='\t')     
    
    # split out spike, real, lone, otherwise same as 32.counts.tsv    
    #spikedf, realdf, lonedf = split_spike_real_lone_barcodes(config, vbcdf)
    spikedf, realdf, lonedf, unmatched = split_spike_real_lone_barcodes(config, vbcdf)
    
    # write out this step...
    realdf.to_csv(os.path.join(outdir , f'{base}.real.counts.tsv'), sep='\t')
    lonedf.to_csv(os.path.join(outdir , f'{base}.lone.counts.tsv'), sep='\t')
    spikedf.to_csv(os.path.join(outdir , f'{base}.spike.counts.tsv'), sep='\t')
    unmatched.to_csv(os.path.join(outdir , f'{base}.unmatched.counts.tsv'), sep='\t')

    # remove homopolymers in real sequences.
    #
    #  base repeats already assumed removed in newer processing. 
    #
    # max_homopolymer_run=int(config.get('ssifasta', 'max_homopolymer_run')) 
    # realdf = remove_base_repeats(realdf, col='sequence', n=max_homopolymer_run)
     
    # align and collapse all.         

    # add labels for merging...
    realdf['type'] = 'real'
    spikedf['type'] = 'spike'
    lonedf['type'] = 'lone'
    unmatched['type'] = 'unmatched'

    outdf = merge_dfs([ realdf, spikedf, lonedf, unmatched ])
    outdf['label'] = base
    outdf.sort_values(by = ['type', 'umi_count'], ascending = [True, False], inplace=True)
    outdf.reset_index(drop=True, inplace=True)
    return outdf



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


def calc_freq_threshold(df, fraction, column):
    '''
    sorts column of input column
    calculates index of point at which <fraction> of data points are less 
    returns column value at that point. 
    '''
    ser = df[column].copy()
    ser.sort_values(ascending = False, inplace=True)
    return 122


def threshold_read_counts(config, df, threshold=1):
    '''
    
    '''
    logging.debug(f'threshold counts threshold={threshold}')
    threshold= int(threshold)   
    df = df[df['read_count'] >= threshold].copy()
    return df


def filter_low_complexity(config, seqdf):
    return seqdf


def split_spike_real_lone_barcodes(config, df):
    '''
    df has  sequence  counts
    should be length 32 ( 30 + YY ) Y= C or T

    # T or C = YY
    # last 2nt of spike-ins AND reals
    
    # A or G = RR
    # last 2nt of L1 controls 
    
    spikeinregex= CGTCAGTC$
    realregex = [TC][TC]$
    loneregex = [AG][AG]$
          
    '''
    #  df[df["col"].str.contains("this string")==False]
    sire = config.get('ssifasta', 'spikeinregex')
    realre = config.get('ssifasta','realregex')
    lonere = config.get('ssifasta', 'loneregex')
    
    logging.debug(f'before filtering: {len(df)}')   
    logging.debug(f"spike-in regex = '{sire}' ")
    simap = df['sequence'].str.contains(sire, regex=True) == True
        
    spikedf = df[simap].copy()
    spikedf.reset_index(inplace=True, drop=True)
    
    remaindf = df[~simap]
    logging.debug(f'spikeins={len(spikedf)} remaindf={len(remaindf)}')
    
    # split real from L1s, and track any that fail both. 
    logging.debug(f"realre = '{realre}' lonere = '{lonere}' ")
    realmap = remaindf['sequence'].str.contains(realre, regex=True) == True
    realdf = remaindf[realmap].copy()
    realdf.reset_index(inplace=True, drop=True)
    
    # remove reals from input. 
    remaindf = remaindf[~realmap]
    logging.debug(f'realdf={len(realdf)} remaindf={len(remaindf)}')
        
    lonemap = remaindf['sequence'].str.contains(lonere, regex=True) == True 
    lonedf = remaindf[lonemap].copy()
    lonedf.reset_index(inplace=True, drop=True)
    logging.debug(f'lonedf={len(lonedf)} remaindf={len(remaindf)}')
    
    unmatcheddf = remaindf[~lonemap].copy()
    unmatcheddf.reset_index(inplace=True, drop=True)
        
    logging.info(f'initial={len(df)} spikeins={len(spikedf)} real={len(realdf)} lone={len(lonedf)} unmatched={len(unmatcheddf)}')    
    return (spikedf, realdf, lonedf, unmatcheddf)


def load_sample_info(config, file_name):
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
        sheet_name = 'Sample information'
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
                logging.warning(traceback.format_exc(None))
                
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
                logging.warn(f'no column {scol}, required. Creating...')
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




def process_fastq_pairs(cp, infilelist, outfile , 
                        force=False,  
                        datestr=None, 
                        max_repeats=None, 
                        max_n_bases=None):
    '''
    parses paired-end fastq files. merges from r1 and r2 at selected positions. 
    outputs to single fasta file. 

    config expectation. 
    [fastq]
    r1start = 0
    r1end = 32
    r2start = 0
    r2end = 20
    # reporting intervals for verbose output, defines how often to log. 
    seqhandled_interval = 1000000
    
    QC drop all sequences with any N in them. 
    Drop all sequences with homopolymer runs > max_homopolymers

    Create read-count frequency TSV and shoulder plot. <outfile>.counts.tsv <outfile>.plot.pdf
    
    '''
    outfile = os.path.abspath(outfile)    
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
        logging.debug(f'made outdir={outdir}')
    
    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")

    sh = StatsHandler(cp, outdir=outdir, datestr=datestr)  

    r1s = int(cp.get('fastq','r1start'))
    r1e = int(cp.get('fastq','r1end'))
    r2s = int(cp.get('fastq','r2start'))
    r2e = int(cp.get('fastq','r2end'))    

    if max_repeats is None:
        max_repeats = int(cp.get('fastq','max_repeats')) 

    if max_n_bases is None:
        max_n_bases = int(cp.get('fastq','max_n_bases'))
        
    seqhandled_interval = int(cp.get('fastq','seqhandled_interval')) 
    
    seqlist = []
    
    if ( not os.path.exists(outfile) ) or force:
        of = open(outfile, 'w')
        pairshandled = 0
        num_handled_total = 0
        num_has_n = 0
        num_has_repeats = 0
        seq_id = 0

        # handle all pairs of readfiles from readfilelist
        for (read1file, read2file) in infilelist:
            pairshandled += 1
            num_handled = 0
            logging.debug(f'handling file pair {pairshandled}')
            if read1file.endswith('.gz'):
                r1f = gzip.open(read1file, "rt")
            else:
                r1f = open(read1file)
            
            if read2file.endswith('.gz'):
                r2f = gzip.open(read2file, "rt")         
            else:
                r2f = open(read2file)
        
            while True:
                try:
                    meta1 = r1f.readline()
                    if len(meta1) == 0:
                        raise StopIteration
                    seq1 = r1f.readline().strip()
                    sep1 = r1f.readline()
                    qual1 = r1f.readline().strip()

                    meta2 = r2f.readline()
                    if len(meta2) == 0:
                        break
                    seq2 = r2f.readline().strip()
                    sep2 = r2f.readline()
                    qual2 = r2f.readline().strip()                    

                    sub1 = seq1[r1s:r1e]
                    sub2 = seq2[r2s:r2e]
                    fullread = sub1 + sub2
                    
                    
                    #has_n = 'N' in fullread
                    has_n = has_n_bases(fullread, n=max_n_bases)
                    if has_n:
                        num_has_n += 1
                        
                    has_repeats = has_base_repeats(fullread, n=max_repeats)
                    if has_repeats:
                        num_has_repeats +=1
                    
                    if has_repeats or has_n:
                        pass
                        #logging.debug(f'seq {num_handled_total} dropped for QC.')
                    else:
                        seqlist.append(fullread)
                        of.write(f'>{seq_id}\n{fullread}\n')
                        seq_id += 1

                    num_handled += 1
                    num_handled_total += 1

                    # report progress...                    
                    if num_handled % seqhandled_interval == 0: 
                        logging.info(f'handled {num_handled} reads from pair {pairshandled}.')
                        logging.info(f'QC: num_has_n={num_has_n} has_repeats={num_has_repeats}')
                
                except StopIteration as e:
                    logging.debug('iteration stopped')    
                    break

            logging.debug(f'finished with pair {pairshandled} saving pair info.')
            # collect pair-specific stats.
            sh.add_value(f'/fastq/pair{pairshandled}','reads_total', num_handled)
        
        of.close()
        # collect global statistics..
        sh.add_value('/fastq','reads_handled', num_handled_total )
        sh.add_value('/fastq','max_repeats', max_repeats  )
        sh.add_value('/fastq','num_has_repeats', num_has_repeats )
        sh.add_value('/fastq','num_has_n', num_has_n )
        sh.add_value('/fastq','num_kept', seq_id )
        logging.info(f'handled {num_handled_total} sequences. {pairshandled} pairs.')
        logging.info(f'QC dropped num_has_n={num_has_n} has_repeats={num_has_repeats} total={num_handled_total} left={seq_id}')
        
        logging.info(f'Creating sequence TSV.')
        sdf = pd.DataFrame(seqlist, columns=['sequence'])
        cdf = make_read_counts_df(cp, sdf)
        write_tsv(cdf, f'{outfile}.counts.tsv')
        make_read_countplot_sns(cp, cdf, f'{outfile}.countsplot.pdf' )
    
    else:
        logging.warn('All FASTA output exists and force=False. Not recalculating.')
        
        

def process_fasta(config, infile, bclist, outdir, force=False, datestr=None):
    '''
    Take paired FASTA file as input rather than raw FASTQ
    Use combined barcode sorter structure to minimize character comparisons. 

    [splitfasta]
    ssifile = ~/git/mapseq-processing/etc/barcode_v2.txt
    seqhandled_interval = 1000000
    matched_interval = 1000000
    unmatched_interval = 1000000

    '''
    if outdir is None:
        outdir = "."
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logging.debug(f'made outdir={outdir}')
    
    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    
    output_exists = check_output(bclist)
    cfilename = f'{outdir}/process_fasta.config.txt'
    bc_length=len(bclist[0].barcode)
    unmatched = os.path.abspath(f'{outdir}/unmatched.fasta')
    
    seqhandled_interval = int(config.get('splitfasta','seqhandled_interval')) 
    matched_interval = int(config.get('splitfasta','matched_interval'))
    unmatched_interval = int(config.get('splitfasta','unmatched_interval'))
    
    logging.info(f'performing split. outdir={outdir} output_exists={output_exists} force={force}')
    
    
    if ( not output_exists ) or force:
        write_config(config, cfilename, timestamp=True, datestring=datestr)
        sh = StatsHandler(config, outdir=outdir, datestr=datestr)        
        # create structure to search for all barcode simultaneously
        labeldict = {}
        for bco in bclist:
            labeldict[bco.barcode] = bco.label
        matchdict, seqdict, unmatched_file = build_bcmatcher(bclist) 

        pairshandled = 0
        num_handled_total = 0
        num_matched_total = 0
        num_unmatched_total = 0

        # handle all sequences in input 
        with open(infile) as f:
            num_handled = 0
            num_matched = 0
            num_unmatched = 0
            
            logging.debug(f'handling file {infile}')
            while True:
                try:
                    line = f.readline()
                    if line.startswith('>'):
                        pass
                    else:
                        if len(line) == 0:
                            raise StopIteration
                        
                        fullread = line.strip()
                        # handle sequence
                        matched = False
                        seq = fullread[-bc_length:]
                        # id, seq, matchdict, fullseq, unmatched
                        matched = do_match(num_handled, seq, matchdict, fullread, unmatched_file)
                        num_handled += 1
                        num_handled_total += 1
    
                        # report progress...                    
                        if not matched:
                            num_unmatched += 1
                            num_unmatched_total += 1                        
                            if num_unmatched % unmatched_interval == 0:
                                logging.debug(f'{num_unmatched} unmatched so far.')
                        else:
                            num_matched += 1
                            num_matched_total += 1
                            if num_matched % matched_interval == 0:
                                logging.debug(f'match {num_matched}: found SSI in {fullread}!')                        
    
                        if num_handled % seqhandled_interval == 0: 
                            logging.info(f'handled {num_handled} matched={num_matched} unmatched={num_unmatched}')
                    
                except StopIteration as e:
                    logging.debug('iteration stopped')    
                    break
            
        logging.debug(f'finished with {infile}')
        sh.add_value('/fasta','reads_handled', num_handled_total )
        sh.add_value('/fasta','reads_unmatched', num_unmatched_total )
        f.close()
        
        matchrate = 0.0
        if num_matched_total > 0: 
            unmatchrate = num_unmatched_total / num_matched_total              
            matchrate = 1.0 - unmatchrate
        logging.info(f'handled {num_handled_total} sequences. {num_matched_total} matched. {num_unmatched_total} unmatched matchrate={matchrate}')
    else:
        logging.warn('All FASTA output exists and force=False. Not recalculating.')
    
   



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


def make_read_countplot_sns(cp, df, outfile='plot.pdf' ):
    '''
    makes two figures, one log-normalized, one straight for same counts df. 
    '''
    from matplotlib.backends.backend_pdf import PdfPages as pdfpages
    
    df.sort_values(by='read_count', ascending=False, inplace=True)
    df.reset_index(inplace=True)
    df['index'] = df['index'].astype('int64')

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
        
        
def counts_axis_plot_sns(ax, df, scale=None ) :
    '''
    Creates individual axes for single plot within figure. 
    scale = None | log10  | log2
    
    
    '''
    s = df['read_count'].sum()
    n = len(df)
    t = df['read_count'].max()
    
    title=''
    if scale is None:
        sns.lineplot(ax=ax, x=df['index'], y=df['read_count'] )
        title = 'counts frequency plot' 
    
    elif scale == 'log10':
        # avoid divide by 0 runtime warning...
        df['log10index'] = np.log10(df['index'] + 1)
        df['log10counts'] = np.log10(df['read_count'])
        sns.lineplot(ax=ax, x=df['log10index'], y=df['log10counts'] )
        title = 'log10(counts) frequency plot' 
         
    ax.set_title(title, fontsize=10)
    ax.text(0.15, 0.2, f"n={n}\ntop={t}\nsum={s}\n", fontsize=9) #add text
   
   

def normalize_weight(df, weightdf, columns=None):
    '''
    Weight values in realdf by spikedf
    Assumes matrix index is sequence.
    Assumes matrices have same columns!!  
    If column numbers are mis-matched, will create empty column
    If columns is none, use/weight all columns, otherwise ignore unlisted columns

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


def process_merged(config, infile, outdir=None, expid=None, recursion=200000, label_column='region' ):
    '''
     takes in combined 'all' TSVs. columns=(sequence, counts, type, label, brain, site) 
     outputs brain-specific SSI x target matrix dataframes, with counts normalized to spikeins by target.  
     writes all output to outdir (or current dir). 
     
    '''
    #sh = get_default_stats()
    logging.debug(f'infile={infile}')
    alldf = load_df(infile)
    logging.debug(f'alldf=\n{alldf}')
    if outdir is None:
        outdir = './'

    require_injection = config.getboolean('merged','require_injection')
    min_injection = int(config.get('merged','min_injection'))
    min_target = int(config.get('merged','min_target'))   
    use_target_negative=config.getboolean('merged','use_target_negative')
    use_target_water_control=config.getboolean('merged','use_target_water_control')
    
    clustermap_scale = config.get('analysis','clustermap_scale')
      
    if expid is None:
        expid = 'M000'
    
    logging.debug(f'running exp={expid} min_injection={min_injection} min_target={min_target} use_target_negative={use_target_negative} use_target_water_control={use_target_water_control}')

    alldf['brain'] = alldf['brain'].astype('string')
    bidlist = list(alldf['brain'].dropna().unique())
    bidlist = [ x for x in bidlist if len(x) > 0 ]
    bidlist.sort()
    logging.debug(f'handling brain list: {bidlist}')
    
    for brain_id in bidlist:
        valid = True
        logging.debug(f'handling brain_id={brain_id}')
        bdf = alldf[alldf['brain'] == brain_id]
                    
        # handle target areas...
        tdf = bdf[bdf['site'].str.startswith('target')]
        rtdf = tdf[tdf['type'] == 'real'] 

        # threshold by min_target or threshold by target-negative
        # if use_target_negative is true, but no target negative site 
        # defined, use min_target and throw warning. 
        if use_target_negative:
            logging.info(f'use_target_negative is {use_target_negative}')
            max_negative = calc_min_threshold(config, bdf, 'target-negative')
            logging.debug(f'target-negative UMI count = {max_negative}')

        if use_target_water_control:
            logging.info(f'use_target_water_control is {use_target_water_control}')
            max_water_control = calc_min_threshold(config, bdf, 'target-water-control')
            logging.debug(f'target_water_control UMI count = {max_water_control}')
        
        min_target = max([ min_target, max_negative, max_water_control ])
        logging.debug(f'min_target UMI count after all constraints = {min_target}')   

        # min_target is now either calculated from target-negative, or from config. 
        if min_target > 1:
            before = len(rtdf)
            frtdf = filter_all_lt(rtdf, 'sequence', 'umi_count', min_target)            
            if not len(frtdf) > 0:
                valid = False
                logging.warning(f'No VBCs passed min_target filtering! Skip brain.')
            logging.debug(f'filtering by min_target={min_target} before={before} after={len(rtdf)}')
        else:
            logging.debug(f'min_target={min_target} no filtering.')
            frtdf = rtdf

        if require_injection:
            # extract and filter injection areas.
            logging.debug(f'require_injection={require_injection} min_injection={min_injection}') 
            idf = bdf[bdf['site'].str.startswith('injection')]
            ridf = idf[idf['type'] == 'real']  
            if len(ridf) == 0:
                logging.warning('require_injection=True but no real VBCs from any injection site.')
            logging.debug(f'{len(frtdf)} real target VBCs before filtering.')      
            frtdf = filter_non_injection(frtdf, ridf, min_injection=min_injection)
            logging.debug(f'{len(frtdf)} real target VBCs after injection filtering.')
            if not len(frtdf) > 0:
                logging.warning(f'No VBCs passed injection filtering! Skip brain.')
                valid = False
        else:
            logging.debug(f'require_injection={require_injection} proceeding...')
               
        # make matrices if brain data is valid... 
        if valid:       
            
            # raw reals
            rbcmdf = rtdf.pivot(index='sequence', columns=label_column, values='umi_count')
            scol = natsorted(list(rbcmdf.columns))
            rbcmdf = rbcmdf[scol]
            rbcmdf.fillna(value=0, inplace=True)
            logging.debug(f'brain={brain_id} raw real barcode matrix len={len(rbcmdf)}')            
            
            # filtered reals
            fbcmdf = frtdf.pivot(index='sequence', columns=label_column, values='umi_count')
            scol = natsorted(list(fbcmdf.columns))
            fbcmdf = fbcmdf[scol]
            fbcmdf.fillna(value=0, inplace=True)
            logging.debug(f'brain={brain_id} real barcode matrix len={len(fbcmdf)}')

            # spikes
            sdf = tdf[tdf['type'] == 'spike']
            sbcmdf = sdf.pivot(index='sequence', columns=label_column, values='umi_count')
            spcol = natsorted(list(sbcmdf.columns))
            sbcmdf = sbcmdf[spcol]
            sbcmdf.fillna(value=0, inplace=True)    
            logging.debug(f'brain={brain_id} spike barcode matrix len={len(sbcmdf)}')
    
            (fbcmdf, sbcmdf) = sync_columns(fbcmdf, sbcmdf)
            
            
            nbcmdf = normalize_weight(fbcmdf, sbcmdf)
            logging.debug(f'nbcmdf.describe()=\n{nbcmdf.describe()}')
            
            scbcmdf = normalize_scale(nbcmdf, logscale=clustermap_scale)
            scbcmdf.fillna(value=0, inplace=True)
            logging.debug(f'scbcmdf.describe()=\n{scbcmdf.describe()}')
                        
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
            
        logging.info(f'done with brain={brain_id}')
        
        


def align_collapse_fasta(config, infile, seq_length=None, max_mismatch=None, outdir=None, datestr=None):
    '''
    Algorithm:
    
    take input FASTA to FULL dataframe. ( already removed Ns and homopolymer runs ) 
    split first seq_length part (VBC) and remainder to 2 columns. 
    
    get UNIQUE DF dataframe on (VBC) sequence      
    save VBCs to fasta
    align all to all bowtie
    make bowtie DF
    drop mismatch > max_mismatch
        for each component, create map from all variants back to a single virtual parent sequence
        (this is not necessarily the true biological parent sequence, but it doesn't matter. they just need to all
        be the same).
    in FULL DF, set all variant sequences to virtual parent sequence
    re-assemble VBC + remainder sequence
    output to adjusted FULL FASTA file ready for input to SSI splitting.  
   
    can use iloc and range of indices to set value for all elements of components:
    
    testdf = pd.DataFrame([[0, 2, 3], [0, 4, 1], [10, 20, 30]],columns=['A', 'B', 'C'])
    testdf
            A   B   C
        0   0   2   3
        1   0   4   1
        2  10  20  30
    testdf.iloc[[0,2],[0]] = 15
            A   B   C
        0  15   2   3
        1   0   4   1
        2  15  20  30

https://stackoverflow.com/questions/46204521/pandas-get-unique-values-from-column-along-with-lists-of-row-indices-where-the
    
    '''
    aligner = config.get('collapse','tool')
    
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    head = base.split('.')[0]
    
    if seq_length is None:
        seq_length = int(config.get('collapse', 'seq_length'))
    
    if max_mismatch is None:
        max_mismatch = int(config.get('collapse', 'max_mismatch'))
    
    if outdir is None:
        outdir = dirname    
    os.makedirs(outdir, exist_ok=True)

    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")

    sh = StatsHandler(config, outdir=outdir, datestr=datestr)
    #  sh.add_value('/collapse','n_components', len(components) )
    
    # need to strip to seq_length
    logging.info(f'Reading {infile} to df...')
    fdf = read_fasta_to_df(infile, seq_length)
    logging.debug(f'fdf=\n{fdf}')
    of = os.path.join( outdir , f'{base}.fulldf.tsv')
    logging.info(f'Writing full DF to {of}')
    fdf.to_csv(of, sep='\t')
    sh.add_value('/collapse','n_full_sequences', len(fdf) )

    # get reduced dataframe of unique sequences
    logging.info('Getting unique DF...')    
    udf = pd.DataFrame(fdf['sequence'].unique(), columns=['sequence'])
    
    sh.add_value('/collapse','n_unique_sequences', len(udf) )
    
    logging.debug(f'udf = \n{udf}')
    of = os.path.join( outdir , f'{base}.udf.tsv')
    logging.info(f'Writing unique DF to {of}')
    udf.to_csv(of, sep='\t') 
    
    of = os.path.join( outdir , f'{base}.udf.fasta')
    logging.info(f'Writing uniques as FASTA to {of}')
    seqfasta = write_fasta_from_df(udf, outfile=of)
    
    # run allXall bowtiex
    of = os.path.join( outdir , f'{base}.bt2.sam')
    logging.info(f'Running {aligner} on {seqfasta} file to {of}')
    # switch to generic bowtie later... JRH
    #afile = run_bowtie(config, seqfasta, seqfasta, of, tool=aligner)
    btfile = run_bowtie(config, seqfasta, of, tool=aligner)
    logging.info(f'Bowtie done. Produced output {btfile}. Creating btdf dataframe...')
    btdf = make_bowtie_df(btfile, max_mismatch=max_mismatch, ignore_self=True)
    of = os.path.join( outdir , f'{base}.btdf.tsv')
    btdf.to_csv(of, sep='\t') 
    sh.add_value('/collapse','n_bowtie_entries', len(btdf) )
    
    # perform collapse...      
    logging.info('Calculating Hamming components...')
    edgelist = edges_from_btdf(btdf)
    btdf = None  # help memory usage
    sh.add_value('/collapse','n_edges', len(edgelist) )
    
    logging.debug(f'edgelist len={len(edgelist)}')
    components = get_components(edgelist)
    logging.debug(f'all components len={len(components)}')
    sh.add_value('/collapse','n_components', len(components) )
    edgelist = None  # help memory usage

    # components is list of lists.
    data = [ len(c) for c in components]
    data.sort(reverse=True)
    ccount = pd.Series(data)
    of = os.path.join( outdir , f'{base}.comp_count.tsv')
    ccount.to_csv(of, sep='\t') 
 
    
    
    components = remove_singletons(components)
    logging.debug(f'multi-element components len={len(components)}')
    sh.add_value('/collapse','n_multi_components', len(components) )
    of = os.path.join( outdir , f'{base}.multi.components.txt')
    writelist(of, components)
    

    logging.info(f'Collapsing {len(components)} components...')
    newdf = collapse_by_components(fdf, udf, components)
    del fdf 
    del udf
    del components
        
    of = os.path.join( outdir , f'{base}.collapsed.tsv')
    logging.info(f'Got collapsed DF. Writing to {of}')
    newdf.to_csv(of, sep='\t')     
    logging.info('Done. Calculating fasta.')
    of = os.path.join( outdir , f'{base}.collapsed.fasta')
    cdf = pd.DataFrame( newdf['sequence'] + newdf['tail'], columns=['sequence'])
    logging.info(f'Writing fasta to {of}')
    write_fasta_from_df(cdf, of)
    logging.info(f'Wrote re-joined sequences to {of}')



def apply_setcompseq(row, seqmapdict ):
    ''' 
      Set Component Sequence
      Applies mapping from old to new sequence from dict. 
      
      seqmapdict = { oldseq : newseq }
      
      Usage example:
            df['sequence'] = df.apply(apply_setcompseq, axis=1, seqmapdict=seqmapdict)
      
    '''
    try:
        #logging.debug(f"getting mapping for {row['sequence']}")
        a = seqmapdict[row['sequence']]
    except KeyError:
        # A KeyError means we don't want to remap, just use original value
        a = row['sequence']
    return a


def build_seqmapdict(udf, components):
    '''
    Create mappings from all unique sequence to component sequence
    '''
    seqmapdict = {}
    comphandled = 0
    comphandled_interval = 100000
    comp_len = len(components)
    
    for comp in components:
        ser = list(udf['sequence'].iloc[comp])
        t = ser[0]
        for s in ser: 
            seqmapdict[s]= t    
        if comphandled % comphandled_interval == 0:
            logging.debug(f'[{comphandled}/{comp_len}] t = {t}')
        comphandled += 1
    return seqmapdict


def collapse_by_components(fulldf, uniqdf, components):
    #
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
    seqmapdict = build_seqmapdict(uniqdf, components)
      
    # Make new full df:
    logging.info('seqmapdict built. Applying.')
    fulldf['sequence'] =fulldf.apply(apply_setcompseq, axis=1, seqmapdict=seqmapdict)
    logging.info(f'New collapsed df = \n{fulldf}')
    log_objectinfo(fulldf, 'fulldf')
    return fulldf

