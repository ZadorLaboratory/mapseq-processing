import gzip
import logging
import os
import sys
import traceback

from configparser import ConfigParser
from collections import defaultdict

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

import pandas as pd
import numpy as np
from natsort import natsorted

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from cshlwork.utils import dataframe_to_seqlist, write_fasta_from_df
from cshlwork.utils import run_command_shell, NonZeroReturnException, merge_dfs 
from cshlwork.utils import merge_tsvs, setup_logging, fix_columns_int
from alignment.bowtie import run_bowtie, make_bowtie_df

from mapseq.barcode import check_output


def get_default_config():
    dc = os.path.expanduser('~/git/mapseq-processing/etc/mapseq.conf')
    cp = ConfigParser()
    cp.read(dc)
    return cp

def process_ssifasta(config, infile, outdir=None):
    '''
    by default, outdir will be same dir as infile
    assumes infile fast has already been trimmed to remove SSI
    '''
    aligner = config.get('ssifasta','tool')
    
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    if outdir is not None:
        dirname = outdir
    
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    logging.debug(f'handling {filepath}')
    
    # make raw fasta TSV of barcode-splitter output for one barcode. 
    # trim to 44 unique w/ counts. 
    logging.info('calc counts...')
    seqdf = make_fasta_df(config, infile)
    of = os.path.join(dirname , f'{base}.seq.raw.tsv')
    seqdf.to_csv(of, sep='\t')
    
    # to calculate threshold we need counts calculated. 
    cdf = make_counts_df(config, seqdf, label=base)  
    logging.info(f'initial counts df {len(cdf)} all reads.')
    of = os.path.join(dirname , f'{base}.44.counts.tsv')
    cdf.to_csv(of, sep='\t') 
        
    threshold = calculate_threshold(config, cdf)
    tdf = threshold_counts(config, cdf)
    logging.info(f'at threshold={threshold} {len(tdf)} unique molecules.')
    
    # now that we've only kept bc + umi, the counts is distinct molecules.  
    tdf['counts'] = 1          
    # trim to just viral barcodes.  ?
    tdf['sequence'] = tdf['sequence'].str[:32]
    
    of = os.path.join(dirname , f'{base}.32.raw.tsv')
    tdf.to_csv(of, sep='\t') 
    
    # now have actual viral barcode df with unique molecule counts.
    bcdf = make_counts_df(config, tdf)
    of = os.path.join(dirname , f'{base}.32.counts.tsv')
    bcdf.to_csv(of, sep='\t')     
    
    # split out spike, real, lone    
    spikedf, realdf, lonedf = split_spike_real_lone_barcodes(config, bcdf)

    realdf.to_csv(os.path.join(dirname , f'{base}.real.raw.tsv'), sep='\t')
    lonedf.to_csv(os.path.join(dirname , f'{base}.lone.raw.tsv'), sep='\t')
    spikedf.to_csv(os.path.join(dirname , f'{base}.spike.raw.tsv'), sep='\t')
    
    realcdf = make_counts_df(config, realdf)
    spikecdf = make_counts_df(config, spikedf)
    lonecdf = make_counts_df(config, lonedf)    

    realcdf.to_csv(os.path.join(dirname , f'{base}.real.counts.tsv'), sep='\t')
    lonecdf.to_csv(os.path.join(dirname , f'{base}.lone.counts.tsv'), sep='\t')
    spikecdf.to_csv(os.path.join(dirname , f'{base}.spike.counts.tsv'), sep='\t')    
            
    acrealdf = align_and_collapse(config, realdf, dirname, base, 'real')
    acspikedf = align_and_collapse(config, spikedf, dirname, base, 'spike')
    aclonedf = align_and_collapse(config, lonedf, dirname, base, 'lone')
    
    acrealdf.to_csv(os.path.join(dirname , f'{base}.real.tsv'), sep='\t')
    acspikedf.to_csv(os.path.join(dirname , f'{base}.spike.tsv'), sep='\t')
    aclonedf.to_csv(os.path.join(dirname , f'{base}.lone.tsv'), sep='\t')

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
    logging.info(f'handling {base} {label}s...')
    aligner = config.get('ssifasta','tool')
    pcthreshold = config.get('ssifasta','post_threshold')
    logging.info(f'{label} {len(countsdf)} sequences, representing {countsdf.counts.sum()} reads.')      
    of = os.path.join( outdir , f'{base}.{label}.seq.fasta')
    logging.debug(f'make fasta for {aligner} = {of}') 
    seqfasta = write_fasta_from_df(config, countsdf, outfile=of)
    of = os.path.join(outdir , f'{base}.{label}.{aligner}')
    logging.info(f'running {aligner}...')
    try:
        afile = run_bowtie(config, seqfasta, of, tool=aligner )  
        logging.info(f'handle {aligner} align file: {afile}')
        btdf = make_bowtie_df(afile)
        of = os.path.join(outdir , f'{base}.{label}.btdf.tsv')
        btdf.to_csv(of, sep='\t') 
        edgelist = edges_from_btdf(btdf)
        components = get_components(edgelist)
        logging.info(f'countdf columns are {countsdf.columns}')
        newdf = collapse_counts_df(countsdf, components)
        newdf = threshold_counts(config, newdf, threshold=pcthreshold)
        logging.info(f'orig len={len(countsdf)}, {len(components)} components, collapsed len={len(newdf)}')
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
    



def collapse_counts_df(countsdf, components):
    '''
    takes components consisting of indices
    determines sequence with largest count columns: 'sequence', 'counts'
    collapses all other member components to the sequence of the largest.
    adds their counts to that of that sequence.

    retain columns and values for highest counts row. 
     
    '''
    # list of lists to collect values..
    logging.info(f'collapsing countsdf len={len(countsdf)} w/ {len(components)} components.')
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
    logging.info(f'original len={len(countsdf)} collapsed len={len(newdf)}')
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
    logging.info(f'getting connected components from edgelist len={len(edgelist)}')
    if len(edgelist) < 100:
        logging.debug(f'{edgelist}')
    for g in trajan(from_edges(edgelist)):
        complist.append(g)
    logging.info(f'{len(complist)} components.')
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
    
    
def trajan(V):
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
    logging.info(seqdf)
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
    logging.info(f"kept {len(slist)} sequences out of {handled}")    
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
    logging.info(f'wrote {len(trimmed)} to {ofpath}')
    return ofpath


def calculate_threshold(config, df, minimum=2):
    '''
    takes counts dataframe (with 'counts' column) 
    and calculates 'shoulder' threshold
    
    '''
    return minimum

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
    if threshold is None:
        count_threshold = int(config.get('ssifasta', 'countthreshold'))
    else:
        count_threshold = int(threshold)
    logging.info(f'thresh = {count_threshold}')    
    df = df[df.counts >= count_threshold].copy()
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


def match_strings(a, b, max_mismatch=0):
    '''
    attempt at efficient comparison of two same-length strings. 
    
    '''
    if len(a) != len(b):
        logging.error(f'two strings must be same length')
    if len(a) <= max_mismatch:
        logging.error(f'length less than mismatch, all will match.')
    mm = 0
    is_match = True
    for i,v in enumerate(a):
        if b[i] == v:
            pass
        else:
            mm += 1
            if mm > max_mismatch:
                is_match = False
                break
    return is_match
      



def load_sample_info(config, file_name):
    print('running new')    
    #['Tube # by user', 'Our Tube #', 'Sample names provided by user',
    #   'Site information', 'RT primers for MAPseq', 'Brain ', 'Column#']
    sheet_columns = ['Tube # by user', 'Our Tube #', 'Sample names provided by user', 'Site information', 'RT primers for MAPseq', 'Brain' ]
    sample_columns = ['usertube', 'ourtube','samplename','siteinfo','rtprimer','brain'] 
    int_sample_col = ['usertube', 'ourtube','rtprimer','brain']
    sheet_name = 'Sample information'
    edf = pd.read_excel(file_name, sheet_name=sheet_name, header=1)        
    sdf = pd.DataFrame()
    
    for i,sc in enumerate(sheet_columns):
        try:
            cser = edf[sc]
            logging.debug(f'column for {sc}:\n{cser}')
            sdf[sample_columns[i]] = cser
        except:
            sdf[sample_columns[i]] = pd.Series(np.nan, np.arange(len(edf))) 
    
    sdf = fix_columns_int(sdf, columns=int_sample_col)
        
    logging.debug(f'created reduced sample info df:\n{sdf}')
    return sdf



   

def process_fastq_pairs(config, readfilelist, bclist, outdir, force=False):

    # if all the output files for bclist exist, don't recalc unless force=True. 
    if outdir is None:
        outdir = "."
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logging.info(f'made outdir={outdir}')
    output_exists = check_output(bclist)
    logging.info(f'output_exists={output_exists} force={force}')
    
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
            logging.info(f'handling file pair {pairshandled}')
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
                        logging.info(f'handled {seqshandled} reads from pair {pairshandled}. matched={didmatch} unmatched={unmatched}')
                
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


def normalize_by_spikeins(realdf, spikedf):
    '''
    Weight values in realdf by spikedf
    
    '''
    logging.debug(f'normalizing arg1 by arg2')
    normdf = realdf
    return normdf

def normalize_target_sum(normdf, columns = None):
    '''
    Normalize so values across rows sum to one across columns 
    given, or all columns if none. 
    
    '''
    logging.debug(f'making rows sum to one...')
    return normdf




def process_merge_areas(config, filelist, outdir=None ):
    '''
     merges SSI-specific real, spike, lone DFs with counts. 
     outputs SSI x target matrix DF, with counts normalized to spikeins by target.  
     
    '''
    logging.info(f'{filelist}')
    alldf = merge_tsvs(filelist)
    
    logging.info(f'alldf len={len(alldf)}')
      
    rdf = alldf[alldf['type'] == 'real']      
    bcm = rdf.pivot(index='sequence', columns='label', values='counts')
    bcm.reset_index(inplace=True)
    bcm.drop(labels=['sequence'], axis=1, inplace=True)
    scol = natsorted(list(bcm.columns))
    bcm = bcm[scol]
    bcm.fillna(value=0, inplace=True)
    logging.info(f'real barcode matrix len={len(bcm)}')
    
    
    sdf = alldf[alldf.type == 'spike']
    sbcm = sdf.pivot(index='sequence', columns='label', values='counts')
    sbcm.reset_index(inplace=True)
    sbcm.drop(labels=['sequence'], axis=1, inplace=True)
    spcol = natsorted(list(sbcm.columns))
    sbcm = sbcm[spcol]
    sbcm.fillna(value=0, inplace=True)    
    logging.info(f'spike barcode matrix len={len(sbcm)}')
        
    return (bcm, sbcm)

