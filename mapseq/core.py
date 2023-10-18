import glob
import gzip
import logging
import math
import os
import subprocess
import sys
import traceback

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

from kneed import KneeLocator

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from mapseq.utils import dataframe_to_seqlist, write_fasta_from_df, remove_base_repeats
from mapseq.utils import run_command_shell, NonZeroReturnException, merge_dfs 
from mapseq.utils import merge_tsvs, setup_logging
from mapseq.utils import JobRunner, JobStack, JobSet

from mapseq.bowtie import run_bowtie, make_bowtie_df

from mapseq.barcode import *

def fix_columns_int(df, columns):
    '''
    forces column in dataframe to be an integer. NaNs become '0'
    Only floating points can be NaN. No good solution for integers...
    '''
    for col in columns:
        try:
            logging.debug(f'trying to fix col {col}')
            fixed = np.array(df[col], np.int16)
            logging.debug(f'fixed=\n{fixed}')
            df[col] = fixed
                
        except ValueError:
            logging.debug(f'invalid literal in {col}')
    return df

def fix_columns_str(df, columns):
    '''
    forces column in dataframe to be string NaNs become ''
    '''
    for col in columns:
        try:
            logging.debug(f'trying to fix col {col}')
            df[col].replace(0,'',inplace=True)
            df[col].replace(np.nan,'', inplace=True)
                 
        except Exception as ex:
            logging.error(f'error while handling {col} ')
            logging.warning(traceback.format_exc(None))
    return df


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
            y = int(x)
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
    rtprimer_num = int(rtprimer_num)
    logging.debug(f'base={base} head={head} guessing rtprimer={rtprimer_num} sampdf=\n{sampdf}')
    
    df = sampdf[sampdf['rtprimer'] == rtprimer_num]
    df.reset_index(inplace=True, drop=True)
    site = None
    if len(df)> 0:
        try:
            site = df['siteinfo'][0]
        except:
            logging.warning(f'unable to get siteinfo for {infile}')
            site = 'target'  # default to target. 
        
        try:        
            brain = df['brain'][0]
        except:
            logging.warning(f'unable to get brain info for {infile}')
            brain = '1'

        try:        
            region = df['region'][0]
        except:
            logging.warning(f'unable to get region info for {infile}')
            region = str(rtprimer_num) # default to SSI label      
            
        logging.debug(f'got site={site} brain={brain} region={region} for rtprimer={rtprimer_num}') 

    logging.debug(f'got site={site} for rtprimer guessed from {infile}')
    return (rtprimer_num, site, brain, region )
    
    
    
def process_ssifasta(config, infile, outdir=None, site=None):
    '''
    by default, outdir will be same dir as infile
    assumes infile fasta has already been trimmed to remove SSI
    
    site = ['target-control','injection-control','target','target-negative',target-lone']   
    Will use relevant threshold. If None, will use default threshold
    
    '''
    aligner = config.get('ssifasta','tool')
    
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    if outdir is not None:
        dirname = outdir
    
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    logging.debug(f'handling {filepath} base={base}')
    
    # make raw fasta TSV of barcode-splitter output for one barcode. 
    # trim to 44 unique w/ counts. 
    logging.debug('calc counts...')
    seqdf = make_fasta_df(config, infile)
    of = os.path.join(dirname , f'{base}.44.seq.tsv')
    seqdf.to_csv(of, sep='\t')
    
    # to calculate threshold we need counts calculated. 
    cdf = make_counts_df(config, seqdf, label=base)  
    logging.debug(f'initial counts df {len(cdf)} all reads.')
    of = os.path.join(dirname , f'{base}.44.counts.tsv')
    cdf.to_csv(of, sep='\t') 
        
    threshold = get_threshold(config, cdf, site)
    logging.debug(f'got threshold={threshold} for site {site}')
    tdf = threshold_counts(config, cdf, threshold=threshold)
    logging.info(f'at threshold={threshold} {len(tdf)} unique molecules.')
    
    # now that we've only kept bc + umi, the counts is distinct molecules. trim to viral barcode only  
    tdf['counts'] = 1          
    tdf['sequence'] = tdf['sequence'].str[:32]
    
    of = os.path.join(dirname , f'{base}.32.seq.tsv')
    tdf.to_csv(of, sep='\t') 
    
    # now have actual viral barcode df with *unique molecule counts.*
    vbcdf = make_counts_df(config, tdf)
    of = os.path.join(dirname , f'{base}.32.counts.tsv')
    vbcdf.to_csv(of, sep='\t')     
    
    # split out spike, real, lone, otherwise same as 32.counts.tsv    
    spikedf, realdf, lonedf = split_spike_real_lone_barcodes(config, vbcdf)
    
    # write out this step...
    realdf.to_csv(os.path.join(dirname , f'{base}.real.seq.tsv'), sep='\t')
    lonedf.to_csv(os.path.join(dirname , f'{base}.lone.seq.tsv'), sep='\t')
    spikedf.to_csv(os.path.join(dirname , f'{base}.spike.seq.tsv'), sep='\t')

    # remove homopolymers in real sequences.
    max_homopolymer_run=int(config.get('ssifasta', 'max_homopolymer_run')) 
    realdf = remove_base_repeats(realdf, col='sequence', n=max_homopolymer_run)
   
    # make counts df
    realcdf = make_counts_df(config, realdf)
    spikecdf = make_counts_df(config, spikedf)
    lonecdf = make_counts_df(config, lonedf)    

    # these are counts of copies of the same UMI, so same original molecule
    realcdf.to_csv(os.path.join(dirname , f'{base}.real.counts.tsv'), sep='\t')
    lonecdf.to_csv(os.path.join(dirname , f'{base}.lone.counts.tsv'), sep='\t')
    spikecdf.to_csv(os.path.join(dirname , f'{base}.spike.counts.tsv'), sep='\t')    
    
    # align and collapse all.         
    acrealdf = align_and_collapse(config, realcdf, dirname, base, 'real')
    acspikedf = align_and_collapse(config, spikecdf, dirname, base, 'spike')
    aclonedf = align_and_collapse(config, lonecdf, dirname, base, 'lone')
    
    acrealdf.to_csv(os.path.join(dirname , f'{base}.real.tsv'), sep='\t')
    acspikedf.to_csv(os.path.join(dirname , f'{base}.spike.tsv'), sep='\t')
    aclonedf.to_csv(os.path.join(dirname , f'{base}.lone.tsv'), sep='\t')

    # add labels for merging...
    acrealdf['type'] = 'real'
    acspikedf['type'] = 'spike'
    aclonedf['type'] = 'lone'
    outdf = merge_dfs([ acrealdf, acspikedf, aclonedf ])
    outdf['label'] = base
    outdf.sort_values(by = ['type', 'counts'], ascending = [True, False], inplace=True)
    outdf.reset_index(drop=True, inplace=True)
    return outdf


def align_and_collapse(config, countsdf, outdir, base, label):
    '''
    countsdf  'sequence' and 'counts' columns
    outdir    working dir or temp dir. 
    base      leading file name, e.g. barcode label, e.g. 'SSI4'
    label     type of sequence, e.g. real, spike, L1 (lone)
    
    '''
    newdf = None
    logging.debug(f'handling {base} {label}s...')
    aligner = config.get('ssifasta','tool')
    logging.info(f'{label} {len(countsdf)} sequences, representing {countsdf.counts.sum()} reads.')      
    of = os.path.join( outdir , f'{base}.{label}.seq.fasta')
    logging.debug(f'make fasta for {aligner} = {of}') 
    seqfasta = write_fasta_from_df(config, countsdf, outfile=of)
    of = os.path.join(outdir , f'{base}.{label}.{aligner}')
    logging.debug(f'running {aligner}...')
    try:
        afile = run_bowtie(config, seqfasta, of, tool=aligner )  
        logging.debug(f'handle {aligner} align file: {afile}')
        btdf = make_bowtie_df(afile)
        of = os.path.join(outdir , f'{base}.{label}.btdf.tsv')
        btdf.to_csv(of, sep='\t') 
        edgelist = edges_from_btdf(btdf)
        components = get_components(edgelist)
        logging.debug(f'countdf columns are {countsdf.columns}')
        newdf = collapse_counts_df(countsdf, components)
        logging.debug(f'orig len={len(countsdf)}, {len(components)} components, collapsed len={len(newdf)}')

    except NonZeroReturnException:
        logging.warning(f'NonZeroReturn Exception. Probably no {label}s found. ')
        newdf = pd.DataFrame(columns = ['sequence','counts'])
    return newdf


def max_hamming(sequence, sequencelist):
    '''
    calculates maximum mismatch between sequence and all sequences in sequencelist. 
    assumes all sequences are same length
    no indels, just substitutions. 
    '''
    #logging.debug(f'seq={sequence}')
    #logging.debug(f'sequencelist={sequencelist}')
    max_dist = 0
    for s in sequencelist:
        dist = 0
        for i in range(0,len(s)):
            if sequence[i] != s[i]:
                dist += 1
        if dist > max_dist:
            max_dist = dist
    return max_dist


def unique_df(seqdf):
    '''
    filters for only unique sequences, sets counts to 1
    '''
    pass



def collapse_counts_df(countsdf, components):
    '''
    takes components consisting of indices
    determines sequence with largest count columns: 'sequence', 'counts'
    collapses all other member components to the sequence of the largest.
    adds their counts to that of that sequence.

    retain columns and values for highest counts row. 
     
    '''
    # list of lists to collect values..
    logging.debug(f'collapsing countsdf len={len(countsdf)} w/ {len(components)} components.')
    lol = []
    colnames = list(countsdf.columns)
    for component in components:
        logging.debug(f'component={component}')
        # make new df of only component sequence rows
        cdf = countsdf.iloc[component].reset_index(drop=True)
        logging.debug(f'cdf=\n{cdf}')
        # which sequence has highest count?
        maxid = cdf.counts.idxmax()
        # extract sequence and count as python list
        row = list(cdf.iloc[maxid])
        # set counts as sum of all collapse sequences. 
        row[1] = cdf.counts.sum()
        lol.append(row)
    
        if logging.getLogger().level <= logging.DEBUG:
            slist = list(cdf['sequence'])
            logging.debug(f'slist={slist}')
            if len(slist) > 1:
                s = row[0]
                maxdiff = max_hamming(s, slist)
                logging.debug(f'max_hamming = {maxdiff} n_seqs={len(slist)}')
            else:
                 logging.debug(f'skip distance calc, one sequence in component.')
        
    newdf = pd.DataFrame(data=lol, columns=colnames)
    logging.debug(f'original len={len(countsdf)} collapsed len={len(newdf)}')
    newdf.sort_values('counts',ascending=False, inplace=True)
    newdf.reset_index(drop=True, inplace=True)
    return newdf


def edges_from_btdf(btdf):
    readlist = btdf.name_read.values.tolist()
    alignlist = btdf.name_align.values.tolist()  
    edgelist = [ list(t) for t in zip(readlist, alignlist)]
    return edgelist

def get_components(edgelist):
    complist = []
    logging.debug(f'getting connected components from edgelist len={len(edgelist)}')
    if len(edgelist) < 100:
        logging.debug(f'{edgelist}')
    for g in tarjan(from_edges(edgelist)):
        complist.append(g)
    logging.debug(f'{len(complist)} components.')
    if len(edgelist) < 100:
        logging.debug(f'{complist}')
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

def make_counts_df(config, seqdf, label=None):
    '''
    input dataframe with 'sequence' column
    make counts column for identical sequences.  
    optionally assign a label to set in new column
    
    '''
    logging.debug(f'seqdf=\n{seqdf}')
    ser = seqdf['sequence'].value_counts()
    df = pd.DataFrame(columns=['sequence','counts'])
    df['sequence'] = ser.index
    df['counts'] = ser.values
    logging.debug(f'counts df = \n{df}')
    if label is not None:
        df['label'] = label
    return df


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


def trim_fasta(config, infile, outdir=None, length=44):
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    if outdir is not None:
        dirname = os.path.abspath(outdir)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)
    head = filename.split('.')[0]    
    logging.debug(f'handling {filepath}')
    
    ofpath = f'{dirname}/{head}.{length}.fasta'
    logging.debug(f'opening {ofpath}...')
    outfile = open(ofpath, 'w')    
    trimmed = []
    sfa = SeqIO.parse(filepath, "fasta")
    for sr in sfa:
        tseq = sr.seq[:length]
        tsr = SeqRecord( tseq, id=sr.id, name=sr.name, description=sr.description)
        trimmed.append(tsr)
    SeqIO.write(trimmed, outfile, 'fasta')
    logging.debug(f'wrote {len(trimmed)} to {ofpath}')
    return ofpath


def cumulative_fract_idx_naive(ser, fract):
    '''
    value at index of row that marks cumulative fraction of total. 
    assumes series sorted in descending order. 
    starts with largest value. 
    '''
    sum_total = ser.sum()
    fraction_int = int(fract * sum_total)
    cum_total = 0
    val = 0
    idx = 0
    for idx in range(0, len(ser)):
        val = ser[idx]
        cum_total = cum_total + val
        fraction =  cum_total / sum_total
        if cum_total > fraction_int:
            break
        else:
            logging.debug(f'at idx={idx} fraction is {fraction}')
    logging.debug(f'val={val} idx={idx} cum_total={cum_total} ')
    return val    
        

def cumulative_fract_idx(ser, fract):
    '''
    value at index of row that marks cumulative fraction of total. 
    assumes series sorted in descending order. 
    starts with largest value. 
    '''
    sum_total = ser.sum()
    fraction_int = int(fract * sum_total)
    cumsum = ser.cumsum()
    ltser = cumsum[cumsum < fraction_int]
    if len(ltser) < 1:
        idx = 0
        val = cumsum
    else:
        idx =  ltser.index[-1]
        val = ser[idx]
    logging.debug(f'val={val} idx={idx} ')
    return val    

def calc_kneed_idx(x, y , inflect, poly=2, sense=4):
    '''
    assumes convex, then concave, decreasing curve.
    inflect = 'knee'|'elbow'
    
    '''
    if inflect == 'knee':
        kl = KneeLocator(x=x, y=y, S=sense, curve='convex',direction='decreasing',interp_method='polynomial',polynomial_degree=poly)
        val = kl.knee
        logging.debug(f'got value {val} for knee from kneed...')
    elif inflect == 'elbow':
        # not validated!
        kl = KneeLocator(x=x, y=y, S=sense, curve='convex',direction='decreasing',interp_method='polynomial',polynomial_degree=poly)        
        val = kl.elbow
        logging.debug(f'got value {val} for knee from kneed...')
    return val



def calc_final_thresholds(config, threshdf):
    '''
    take threshold df for all sites, and derive final thresholds df for
    
    threshdf columns used:   site  count_threshold   
    
    target_threshold = 100
    target-control_threshold = 1000
    target-negative_threshold = 100
    target-lone_threshold = 100
    injection_threshold = 2
    injection-control_threshold=2
    
    
        'site'  'threshold'
    
    
    '''
    tdf = pd.concat( [threshdf[threshdf['site'] == 'target-negative'], 
                      threshdf[threshdf['site'] == 'target']] )
    idf = threshdf[threshdf['site'] == 'injection']
    
    target_thresh = int(tdf['count_threshold'].min())
    inj_thresh = int(idf['count_threshold'].min())

    finaldf = pd.DataFrame( data=[ ['target', target_thresh ],['injection', inj_thresh] ], 
                            columns= ['site','threshold'] )
                  
    return finaldf
    




def calc_thresholds_all(config, sampdf, filelist, fraction=None ):
    '''
    reads in all counts.df (assumes counts column).
     
    calculates thresholds for 'target' and 'injection'
    
    returns 2 dfs. one general info, one with final thresholds
    '''
    if fraction is not None:
        config.set('ssifasta','count_threshold_fraction', fraction)
    
    outlist = []
    
    for filename in filelist:
       logging.debug(f'handling {filename}') 
       (rtprimer, site, brain, region) = guess_site(filename, sampdf)
       cdf = pd.read_csv(filename ,sep='\t', index_col=0)
       (count_threshold, label, clength, counts_max, counts_min)   = calculate_threshold(config, cdf, site )
       outlist.append( [rtprimer, site, count_threshold, label, clength, counts_max, counts_min    ])
    threshdf = pd.DataFrame(data=outlist, columns=['rtprimer', 'site', 'count_threshold', 'label', 'counts_length', 'counts_max', 'counts_min'  ])
    finaldf = calc_final_thresholds(config, threshdf)   
    
    return (finaldf, threshdf)
     
    

def calculate_threshold(config, cdf, site=None):
    '''
    takes counts dataframe (with 'counts' column) 
    if 'label', use that. 
    and calculates 'shoulder' threshold
    site = ['control','injection','target']   
        Will use relevant threshold. If None, will use default threshold
    
    target_threshold=100
    target_ctrl_threshold=1000
    inj_threshold=2
    inj_ctrl_threshold=2
    
    '''
    count_pct = float(config.get('ssifasta','count_threshold_fraction'))
    min_threshold = int(config.get('ssifasta','count_threshold_min'))
    label = 'BCXXX'
    
    try:
        label = cdf['label'].unique()[0]
    except:
        logging.warn(f'no SSI label in DF')
        
    # assess distribution.
    counts = cdf['counts']
    clength = len(counts)
    counts_max = counts.max()
    counts_min = counts.min()
    counts_mean = counts.mean()
    logging.info(f'handling {label} length={clength} max={counts_max} min={counts_min} ')
    
    val =  cumulative_fract_idx(counts, count_pct)
    if val < min_threshold:
        logging.warning(f'calc threshold < min...')
    else:
        logging.debug(f'calculated count threshold={val} for SSI={label}')
    count_threshold=max(val, min_threshold)
    
    
    #if site is None:
    #    count_threshold = int(config.get('ssifasta', 'default_threshold'))
    #else:
    #    count_threshold = int(config.get('ssifasta', f'{site}_threshold'))
    #logging.debug(f'count threshold for {site} = {count_threshold}')
    return (count_threshold, label, clength, counts_max, counts_min)



def calculate_threshold_kneed(config, cdf, site=None, inflect=None ):
    '''
    takes counts dataframe (with 'counts' column) 
    if 'label', use that. 
    and calculates 'knee' or 'elbow' threshold
    site = ['control','injection','target']   
        Will use relevant threshold. If None, will use default threshold
    
    target_threshold=100
    target_ctrl_threshold=1000
    inj_threshold=2
    inj_ctrl_threshold=2
    
    '''
    if inflect is None:
        inflect = config.get('ssifasta','threshold_heuristic')
    min_threshold = int(config.get('ssifasta','count_threshold_min'))
    label = 'BCXXX'
    
    try:
        label = cdf['label'].unique()[0]
    except:
        logging.warn(f'no SSI label in DF')
        
    # assess distribution.
    counts = cdf['counts']
    clength = len(counts)
    counts_max = counts.max()
    counts_min = counts.min()
    counts_mean = counts.mean()
    logging.info(f'handling {label} length={clength} max={counts_max} min={counts_min} ')
    
    val = calc_kneed_idx(cdf.index, cdf.counts, inflect='knee'  )
    if val < min_threshold:
        logging.warning(f'kneed calc threshold < min...')
    else:
        logging.debug(f'calculated count threshold={val} for SSI={label}')
    count_threshold=max(val, min_threshold)
        
    #if site is None:
    #    count_threshold = int(config.get('ssifasta', 'default_threshold'))
    #else:
    #    count_threshold = int(config.get('ssifasta', f'{site}_threshold'))
    #logging.debug(f'count threshold for {site} = {count_threshold}')
    return (count_threshold, label, clength, counts_max, counts_min)


def get_threshold(config, cdf, site=None):
    '''
    site = ['control','injection','target']   
        Will use relevant threshold. If None, will use default threshold
    
    target_threshold=100
    target_ctrl_threshold=1000
    inj_threshold=2
    inj_ctrl_threshold=2
    
    '''
    count_threshold=2
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


def threshold_counts(config, df, threshold=None):
    '''
    
    '''
    logging.debug(f'threshold counts threshold={threshold}')
    threshold= int(threshold)   
    df = df[df['counts'] >= threshold].copy()
    return df


def filter_low_complexity(config, seqdf):
    return seqdf


def split_spike_real_lone_barcodes(config, df):
    '''
    df has  sequence  counts
    should be length 32 ( 30 + YY ) Y= C or T
          
    '''
    #  df[df["col"].str.contains("this string")==False]
    sire = config.get('ssifasta', 'spikeinregex')
    realre = config.get('ssifasta','realregex')
    lonere = config.get('ssifasta', 'loneregex')
    
    logging.debug(f'before filtering: {len(df)}')
    
    # get spikeins
    siseq = 'CGTCAGTC'
    logging.debug(f"spike-in regex = '{sire}' ")
    #simap = df['sequence'].str.contains(sire, regex=True) == True
    simap = df['sequence'].str.endswith(siseq) == True
    spikedf = df[simap]
    spikedf.reset_index(inplace=True, drop=True)
    remaindf = df[~simap]
    logging.debug(f'spikeins={len(spikedf)} remaindf={len(remaindf)}')
    
    # split real/L1
    logging.debug(f"realre = '{realre}' lonere = '{lonere}' ")
    realmap = remaindf['sequence'].str.contains(realre, regex=True) == True
    lonemap = remaindf['sequence'].str.contains(lonere, regex=True) == True 
    
    realdf = remaindf[realmap]
    realdf.reset_index(inplace=True, drop=True)
    lonedf = remaindf[lonemap]
    lonedf.reset_index(inplace=True, drop=True)
    logging.info(f'initial={len(df)} spikeins={len(spikedf)} real={len(realdf)} lone={len(lonedf)}')    
    return (spikedf, realdf, lonedf)


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
        }
    
    sample_columns = ['usertube', 'ourtube', 'samplename', 'siteinfo', 'rtprimer', 'brain', 'region'] 
    int_sample_col = ['usertube', 'ourtube', 'rtprimer']     # brain is often not a number. 
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

            #sdf[sample_columns[i]] = pd.Series(np.nan, np.arange(len(edf))) 
        for scol in sample_columns:
            try:
                ser = sdf[scol]
            except KeyError as ke:
                logging.warn(f'no column {scol}, required. Creating...')
                if scol == 'samplename':
                    sdf[scol] = sdf['ourtube']
                
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


   

def process_fastq_pairs(config, sampdf, readfilelist, bclist, outdir, force=False, countsplots=True):

    # if all the output files for bclist exist, don't recalc unless force=True. 
    if outdir is None:
        outdir = "."
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logging.debug(f'made outdir={outdir}')
    output_exists = check_output(bclist)
    logging.debug(f'output_exists={output_exists} force={force}')
    
    if ( not output_exists ) or force:
        outfile = os.path.abspath(f'{outdir}/unmatched.fasta')
        pairedfile = os.path.abspath(f'{outdir}/paired.txt')
        umf = open(outfile, 'w')
        pf = open(pairedfile, 'w')
        r1s = int(config.get('fastq','r1start'))
        r1e = int(config.get('fastq','r1end'))
        r2s = int(config.get('fastq','r2start'))
        r2e = int(config.get('fastq','r2end'))
        
        seqhandled_interval = int(config.get('fastq','seqhandled_interval')) 
        matched_interval = int(config.get('fastq','matched_interval'))
        unmatched_interval = int(config.get('fastq','unmatched_interval'))

        seqshandled = 0
        pairshandled = 0
        unmatched = 0
        didmatch = 0
    
        #
        # handle pairs of readfiles from readfilelist
        #
        for (read1file, read2file) in readfilelist:
            pairshandled += 1
            logging.debug(f'handling file pair {pairshandled}')
            if read1file.endswith('.gz'):
                 read1file = gzip.open(read1file, "rt")
            if read2file.endswith('.gz'):
                 read2file = gzip.open(read2file, "rt")         
                
            recs1 = SeqIO.parse(read1file, "fastq")
            recs2 = SeqIO.parse(read2file, "fastq")
        
            while True:
                try:
                    r1 = next(recs1)
                    r2 = next(recs2)
                    sub1 = r1.seq[r1s:r1e]
                    sub2 = r2.seq[r2s:r2e]
                    fullread = sub1 + sub2
                    pf.write(f'{fullread}\n')
                    
                    matched = False
                    for bch in bclist:
                        r = bch.do_match(seqshandled, fullread)
                        if r:
                            didmatch += 1
                            if didmatch % matched_interval == 0:
                                logging.debug(f'match {didmatch}: found SSI {bch.label} in {fullread}!')
                            matched = True
                            break
                    if not matched:
                        unmatched += 1
                        if unmatched % unmatched_interval == 0:
                            logging.debug(f'{unmatched} unmatched so far.')
                        id = str(seqshandled)
                        sr = SeqRecord( fullread, id=id, name=id, description=id)
                        SeqIO.write([sr], umf, 'fasta')
                    
                    seqshandled += 1
                    if seqshandled % seqhandled_interval == 0: 
                        logging.debug(f'handled {seqshandled} reads from pair {pairshandled}. matched={didmatch} unmatched={unmatched}')
                
                except StopIteration as e:
                    logging.debug(f'iteration stopped')
                    break
                
        
        umf.close()
        pf.close()
        for bch in bclist:
            bch.finalize()    
        # close possible gzip filehandles??
        #max_mismatch = bclist[0].max_mismatch
        logging.info(f'handled {seqshandled} sequences. {pairshandled} pairs. {didmatch} matched. {unmatched} unmatched')
    else:
        logging.warn('all output exists and force=False. Not recalculating.')
    
    filelist = []
    for bch in bclist:
        filelist.append(bch.filename)
    logging.info(f'Making counts df for {filelist} in {outdir}')
    make_counts_dfs(config, filelist, outdir)

    # by default make countsplots 
    if countsplots:
        logging.info('Making combined countsplots PDF...')
        countsfilelist = []
        for bch in bclist:
            dirname = os.path.dirname(bch.filename)
            filename = os.path.basename(bch.filename)
            (base, ext) = os.path.splitext(filename)   
            of = os.path.join(dirname , f'{base}.44.seq.tsv')
            countsfilelist.append(of)
        make_countsplot_combined_sns(config, sampdf, countsfilelist, outfile=None, expid=None )



def calc_thread_count(nthreads):
    ncpus = os.cpu_count()
    threads = 1 # safe default
    if nthreads > 1:
        # use nthreads CPUS up to ncpus.
        threads = min(nthreads, ncpus)
        logging.debug(f'nthreads positive. use {threads}') 
    elif nthreads < 0:
        # use all but -N CPUs
        threads = ncpus - abs(nthreads)
        threads = max(threads, 1)
        logging.debug(f'nthreads negative. Use all but {abs(nthreads)}')
    else:
        # ntrheads is 0
        logging.debug(f'nthreads = 0, use all CPUS: {ncpus}')
        threads = ncpus
    return (ncpus, threads)


def process_fastq_pairs_parallel(config, readfilelist, bclist, outdir, nthreads, force=False):
    '''
    
    nthreads:    use this number of CPUs. 0 means all. -1 means all but 1. 3 means 3. 
    
    '''
    ncpus, threads = calc_thread_count(nthreads)   
    logging.info(f'using {threads} of {ncpus} CPUs in parallel...')

    prog = os.path.expanduser('~/git/mapseq-processing/scripts/process_fastq_thread.py')
    readfilestr = ""
    for (a,b) in readfilelist:
        readfilestr += f" {a} {b} "
    
    logging.debug(f'readfilestr = {readfilestr}')
    
    # from cshlwork.utils import JobRunner, JobStack, JobSet
    
    jstack = JobStack()
    
    for bco in bclist:
        cmd = [ prog, 
               '-d',
               '-B', bco.label , 
               '-O' , outdir ,
               readfilestr 
               ]
        jstack.addjob(cmd)
    jset = JobSet(max_processes = threads, jobstack = jstack)
    jset.runjobs()
    
    
    filelist = []
    for bco in bclist:
        filelist.append(bco.filename)
    logging.info(f'Making counts df for {filelist} in {outdir}')
    make_counts_dfs(config, filelist, outdir)


def process_fastq_pairs_single(config, readfilelist, bclist, outdir, force=False):

    # if all the output files for bclist exist, don't recalc unless force=True. 
    if outdir is None:
        outdir = "."
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logging.debug(f'made outdir={outdir}')
    output_exists = check_output(bclist)
    logging.debug(f'output_exists={output_exists} force={force}')
    
    # list should have only one...
    bcho = bclist[0]
    
    if ( not output_exists ) or force:
        #outfile = os.path.abspath(f'{outdir}/unmatched.fasta')
        #pairedfile = os.path.abspath(f'{outdir}/paired.txt')
        #umf = open(outfile, 'w')
        #pf = open(pairedfile, 'w')
        r1s = int(config.get('fastq','r1start'))
        r1e = int(config.get('fastq','r1end'))
        r2s = int(config.get('fastq','r2start'))
        r2e = int(config.get('fastq','r2end'))
        
        seqhandled_interval = int(config.get('fastq','seqhandled_interval')) 
        matched_interval = int(config.get('fastq','matched_interval'))
        unmatched_interval = int(config.get('fastq','unmatched_interval'))

        seqshandled = 0
        pairshandled = 0
        unmatched = 0
        didmatch = 0
    
        #
        # handle pairs of readfiles from readfilelist
        #
        for (read1file, read2file) in readfilelist:
            pairshandled += 1
            logging.debug(f'handling file pair {pairshandled}')
            if read1file.endswith('.gz'):
                 read1file = gzip.open(read1file, "rt")
            if read2file.endswith('.gz'):
                 read2file = gzip.open(read2file, "rt")         
                
            recs1 = SeqIO.parse(read1file, "fastq")
            recs2 = SeqIO.parse(read2file, "fastq")
        
            while True:
                try:
                    r1 = next(recs1)
                    r2 = next(recs2)
                    sub1 = r1.seq[r1s:r1e]
                    sub2 = r2.seq[r2s:r2e]
                    fullread = sub1 + sub2
                    #pf.write(f'{fullread}\n')
                    
                    matched = False
                    r = bcho.do_match(seqshandled, fullread)
                    if r:
                        didmatch += 1
                        if didmatch % matched_interval == 0:
                            logging.debug(f'match {didmatch}: found SSI {bcho.label} in {fullread}!')
                        matched = True
                    else:
                        unmatched += 1
                        # when processing single, unmatched number not useful. 
                        #if unmatched % unmatched_interval == 0:
                        #    logging.debug(f'{unmatched} unmatched so far.')
                        #id = str(seqshandled)
                        #sr = SeqRecord( fullread, id=id, name=id, description=id)
                        #SeqIO.write([sr], umf, 'fasta')
                    
                    seqshandled += 1
                    if seqshandled % seqhandled_interval == 0: 
                        logging.debug(f'handled {seqshandled} reads from pair {pairshandled}. matched={didmatch} unmatched={unmatched}')
                
                except StopIteration as e:
                    logging.debug(f'iteration stopped?')
                    logging.warning(traceback.format_exc(None))
                    break
                
        #umf.close()
        #pf.close()
        #for bch in bclist:
        bcho.finalize()    
        # close possible gzip filehandles??
        #max_mismatch = bclist[0].max_mismatch
        logging.info(f'handled {seqshandled} sequences. {pairshandled} pairs. {didmatch} matched. {unmatched} unmatched')
    else:
        logging.warn('all output exists and force=False. Not recalculating.')
     

def make_counts_dfs(config, filelist, outdir):
    '''
    
    '''
    dflist = []
    for filepath in filelist:
        logging.debug(f'calculating counts for  file {filepath} ...')    
        dirname = os.path.dirname(filepath)
        
        if outdir is not None:
            dirname = outdir
        
        filename = os.path.basename(filepath)
        (base, ext) = os.path.splitext(filename)   
        logging.debug(f'handling {filepath} base={base}')
        
        # make raw fasta TSV of barcode-splitter output for one barcode. 
        # trim to 44 unique w/ counts. 
        seqdf = make_fasta_df(config, filepath)
        of = os.path.join(dirname , f'{base}.44.seq.tsv')
        seqdf.to_csv(of, sep='\t')
        
        # to calculate threshold we need counts calculated. 
        cdf = make_counts_df(config, seqdf, label=base)  
        logging.debug(f'initial counts df {len(cdf)} all reads.')
        of = os.path.join(dirname , f'{base}.44.counts.tsv')
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
    

def make_countsplots(config, filelist ): 
    '''
    makes individual read counts plots from 44.counts.tsv files. 
    
    '''   
    for bcfile in filelist:
        logging.debug(f'handling {bcfile}')
        filepath = os.path.abspath(bcfile)    
        dirname = os.path.dirname(filepath)   
        filename = os.path.basename(filepath)
        (base, ext) = os.path.splitext(filename)
        base = base.split('.')[0] 
        
        bcdata = pd.read_csv(bcfile, sep='\t')
        plt.figure()
        plt.plot(np.log10(bcdata['Unnamed: 0']), np.log10(bcdata['counts']))
        plt.title(base)
        plt.xlabel("log10(BC index)")
        plt.ylabel("log10(BC counts)")
        plt.savefig(bcfile.replace('tsv', 'pdf'))


def counts_axis_plot_sns(ax, bcdata, labels):
    '''
    Creates individual axes for single plot within figure. 
    
    '''
    bcdata['log10index'] = np.log10(bcdata.index)
    bcdata['log10counts'] = np.log10(bcdata['counts'])
    sns.lineplot(ax=ax, x=bcdata['log10index'], y=bcdata['log10counts'] )
    s = bcdata.counts.sum()
    n = len(bcdata)
    t = bcdata.counts.max()
    

    title = f"{bcdata['label'][0]}"      
    ax.set_title(title, fontsize=10)
    ax.text(0.15, 0.2, f"site={labels['site']}\nn={n}\ntop={t}\nsum={s}\nthreshold={labels['threshold']}", fontsize=9) #add text
    
    #sns.move_legend(ax, "lower left")
    #ax.set_xlabel("log10(BC index)", fontsize=5)
    #ax.set_ylabel("log10(BC counts)",fontsize=5)


def make_countsplot_combined_sns(config, sampdf, filelist, outfile=None, expid=None ):    
    '''
     makes combined figure with all plots. 
     assumes column 'label' for title. 
     
    '''
    min_ssi_count = int(config.get('analysis','min_ssi_count')) 
    
    from matplotlib.backends.backend_pdf import PdfPages as pdfpages
    
    if outfile is None:
        outfile = 'countsplots.pdf'
        if expid is not None:
            outfile = f'{expid}_{outfile}'
    
    # do nine per figure...
    page_dims = (11.7, 8.27)
    with pdfpages(outfile) as pdfpages:
        #fig_n = math.ceil( math.sqrt(len(filelist)) )
        #fig, axes = plt.subplots(nrows=fig_n, ncols=fig_n, figsize=a4_dims,  layout='constrained')
        plots_per_page = 9
        num_figs = float(len(filelist)) / float(plots_per_page)
        if num_figs % 9 == 0:
            num_figs = int(num_figs)
        else:
            num_figs = int(num_figs) + 1
        logging.debug(f'with {plots_per_page} plots/page, need {num_figs} for {len(filelist)} file plots.')
        
        figlist = []
        axlist = []
        for i in range(0,num_figs):
            fig,axes = plt.subplots(nrows=3, ncols=3, figsize=page_dims,  layout='constrained') 
            if expid is not None:
                fig.suptitle(f'{expid} read counts frequency plots.')
            else:
                fig.suptitle(f'Read counts frequency plots')
            figlist.append(fig)
            # numpy.flatirator doesn't handle indexing
            for a in axes.flat:
                axlist.append(a)
        logging.debug(f'created {len(figlist)} figures to go on {num_figs} pages. ')
                  
        #fig.set_xlabel("log10(BC index)")
        #fig.set_ylabel("log10(BC counts)")
        filelist = natsorted(filelist)
        logging.debug(f'handling {len(filelist)} files...')            
        for i, bcfile in enumerate(filelist):
            logging.debug(f'handling {bcfile}')
            bcdata = pd.read_csv(bcfile, sep='\t')
            if len(bcdata) > min_ssi_count:
                (rtprimer_num, site, brain, region ) = guess_site(bcfile, sampdf )           
                count_threshold, label, clength, counts_max, counts_min = calculate_threshold(config, bcdata)
                labels = {'rtprimer':rtprimer_num,
                          'site':site,
                          'brain':brain,
                          'region': region,
                          'threshold' : count_threshold
                          }
                
                ax = axlist[i]
                counts_axis_plot_sns(ax, bcdata, labels=labels)
            else:
                ax = axlist[i]
                # make empty axis?
    
        for f in figlist:
            pdfpages.savefig(f)
    logging.info(f'saved plot PDF to {outfile}')
    



def normalize_weight(df, weightdf, columns=None):
    '''
    Weight values in realdf by spikedf
    Assumes matrix index is sequence.
    Assumes matrices have same columns!!  
    If column numbers are mis-matched, will create empty column
    If columns is none, use/weight all columns, otherwise ignore unlisted columns
    
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
    logging.debug(f'before threshold inj df len={len(ridf)}')
    ridf = ridf[ridf.counts >= min_injection]
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

def filter_min_target(df, min_target=1):
    '''
    
    '''
    logging.debug(f'before threshold inj df len={len(ridf)}')
    ridf = ridf[ridf.counts >= min_injection]
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



def process_merged(config, filelist, outdir=None, expid=None, recursion=200000, combined_pdf=True, label_column='region' ):
    '''
     takes in combined 'all' TSVs. columns=(sequence, counts, type, label, brain, site) 
     outputs brain-specific SSI x target matrix DF, with counts normalized to spikeins by target.  
     writes all output to outdir (or current dir). 
     
    '''
    
    from matplotlib.backends.backend_pdf import PdfPages as pdfpages
    sys.setrecursionlimit(recursion)    
    
    logging.debug(f'{filelist}')
    
    alldf = merge_tsvs(filelist)
    logging.debug(f'alldf len={len(alldf)}')
    
    cmap = config.get('plots','heatmap_cmap')
    require_injection = config.getboolean('analysis','require_injection')
    min_injection = int(config.get('analysis','min_injection'))
    min_target = int(config.get('analysis','min_target'))   
    clustermap_scale = config.get('plots','clustermap_scale') # log10 | log2
    
      
    if expid is None:
        expid = 'MAPseq'
            
    outfile = f'{expid}.all.heatmap.pdf'
    if require_injection:
        outfile = f'{expid}.all.{min_injection}.{min_target}.{clustermap_scale}.pdf'
    else:
        outfile = f'{expid}.all.noinj.{min_target}.{clustermap_scale}.pdf'
    
    logging.debug(f'running exp={expid} min_injection={min_injection} min_target={min_target} cmap={cmap} clustermap_scale={clustermap_scale} ')
    
    page_dims = (11.7, 8.27)
    with pdfpages(outfile) as pdfpages:
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

            #threshold by min_target ...
            if min_target > 1:
                before = len(rtdf)
                rtdf = rtdf[rtdf['counts'] >= min_target]
                rtdf.reset_index(inplace=True, drop=True)
                logging.debug(f'filtering by min_target={min_target} before={before} after={len(rtdf)}')
            else:
                logging.debug(f'min_target={min_target} no filtering.')
            
            if require_injection:
                # extract and filter injection areas.
                logging.debug(f'require_injection={require_injection} min_injection={min_injection}') 
                idf = bdf[bdf['site'].str.startswith('injection')]
                ridf = idf[idf['type'] == 'real']  
                if len(ridf) == 0:
                    logging.warning('require_injection=True but no real VBCs from any injection site.')
                logging.debug(f'{len(rtdf)} real target VBCs before filtering.')      
                frtdf = filter_non_injection(rtdf, ridf, min_injection=min_injection)
                logging.debug(f'{len(rtdf)} real target VBCs after injection filtering.')
                if not len(rtdf) > 0:
                    logging.warning(f'No VBCs passed injection filtering! Skip brain.')
                    valid = False
            else:
                logging.debug(f'require_injection={require_injection} proceeding...')
                frtdf = rtdf
            
            
            # make 
            if valid:       
                rbcmdf = frtdf.pivot(index='sequence',columns=label_column, values='counts')
                scol = natsorted(list(rbcmdf.columns))
                rbcmdf = rbcmdf[scol]
                rbcmdf.fillna(value=0, inplace=True)
                logging.debug(f'brain={brain_id} real barcode matrix len={len(rbcmdf)}')
                # spikes
                sdf = tdf[tdf['type'] == 'spike']
                sbcmdf = sdf.pivot(index='sequence', columns=label_column, values='counts')
                spcol = natsorted(list(sbcmdf.columns))
                sbcmdf = sbcmdf[spcol]
                sbcmdf.fillna(value=0, inplace=True)    
                logging.debug(f'brain={brain_id} spike barcode matrix len={len(sbcmdf)}')
        
                (rbcmdf, sbcmdf) = sync_columns(rbcmdf, sbcmdf)
                
                nbcmdf = normalize_weight(rbcmdf, sbcmdf)
                logging.debug(f'nbcmdf.describe()=\n{nbcmdf.describe()}')
                scbcmdf = normalize_scale(nbcmdf, logscale=clustermap_scale)
                scbcmdf.fillna(value=0, inplace=True)
                logging.debug(f'scbcmdf.describe()=\n{scbcmdf.describe()}')
                
                rbcmdf.to_csv(f'{outdir}/{brain_id}.rbcm.tsv', sep='\t')
                sbcmdf.to_csv(f'{outdir}/{brain_id}.sbcm.tsv', sep='\t')    
                nbcmdf.to_csv(f'{outdir}/{brain_id}.nbcm.tsv', sep='\t')
                scbcmdf.to_csv(f'{outdir}/{brain_id}.scbcm.tsv', sep='\t')
                
                # check to ensure no columns are missing barcodes.
                droplist = []
                for c in scbcmdf.columns:
                    if not scbcmdf[c].sum() > 0:
                        logging.warn(f'columns {c} for brain {brain_id} has no barcodes, dropping...')
                        droplist.append(c)
                logging.debug(f'dropping columns {droplist}')
                scbcmdf.drop(droplist,inplace=True, axis=1 )    
                
                
                try:
                    kws = dict(cbar_kws=dict(orientation='horizontal'))  
                    g = sns.clustermap(scbcmdf, cmap=cmap, yticklabels=False, col_cluster=False, standard_scale=1, **kws)
                    #g.ax_cbar.set_title('scaled log10(cts)')
                    x0, _y0, _w, _h = g.cbar_pos
                    #g.ax_cbar.set_position((0.8, .2, .03, .4))
                    g.ax_cbar.set_position([x0, 0.9, g.ax_row_dendrogram.get_position().width, 0.05])
                    g.fig.suptitle(f'{expid} {brain_id}')
                    g.ax_heatmap.set_title(f'Scaled {clustermap_scale}(counts)')
                    plt.savefig(f'{outdir}/{brain_id}.{clustermap_scale}.clustermap.pdf')
                    if combined_pdf:
                        logging.info(f'saving plot to {outfile} ...')
                        pdfpages.savefig(g.fig)
                except Exception as ee:
                    logging.warning(f'Unable to clustermap plot for {brain_id}. Message: {ee}')
                    
            logging.info(f'done with brain={brain_id}')



def process_qc(config, exp_dir ):
    '''
    consume files in standard directories and generate qc stats and analysis report..
    
    '''
    pass


def process_mapseq_dir(exp_id, loglevel, force):
    '''
    
    process_fastq.py -d/-v -force -b barcode_v2.<dir>.txt -s <dir>_sampleinfo.jrh.xlsx -O fastq.out fastq/*.fastq*  
    process_ssifasta.py -d/-v -s <dir>_sampleinfo.jrh.xlsx -o <dir>.all.tsv -O ssi.out 
    process_merged.py -d/-v  -e <dir> --combined -s <dir>_sampleinfo.jrh.xlsx -O merged.out 
           
    '''   
    d = os.path.abspath(exp_id)
    if not os.path.exists(d):
        sys.exit(f'Experiment directory {d} does not exist.')

    logging.info(f'processing experiment dir: {d}')
    
    config = get_default_config()
    expconfig = f'{d}/mapseq.conf'

    if os.path.exists(expconfig):
         config.read(expconfig)
         logging.debug(f'read {expconfig}')
    
         
    try:
       samplefile = f'{d}/{exp_id}_sampleinfo.jrh.xlsx'
       sampdf = load_sample_info(config, samplefile)
       rtlist = get_rtlist(sampdf)
       
       outdir = f'{d}/fastq.out'
       bcfile = f'{d}/barcode_v2.{exp_id}.txt'       
       bclist = load_barcodes(config, bcfile, labels=rtlist, outdir=outdir)
       readfilelist = package_pairfiles( glob.glob(f'{d}/fastq/*.fastq*'))

       logging.info(f'running process_fastq_pairs. readfilelist={readfilelist} outdir={outdir}')
       process_fastq_pairs(config, readfilelist, bclist, outdir, force=False)
       
         
       #process_ssifasta(config, infile, outdir=None, site=None)
       #process_merged(config, filelist, outdir=None, expid=None, recursion=100000, combined_pdf=True)
       #process_qc(config, exp_dir)

    except Exception as ex:
        logging.error(f'error while handling {d} ')
        logging.warning(traceback.format_exc(None))
        



