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

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from cshlwork.utils import dataframe_to_seqlist, write_fasta_from_df, run_command_shell, NonZeroReturnException, setup_logging
from alignment.bowtie import run_bowtie, make_bowtie_df



def get_default_config():
    dc = os.path.expanduser('~/git/mapseq-processing/etc/mapseq.conf')
    cp = ConfigParser()
    cp.read(dc)
    return cp


def process_bcfasta(config, infile, outdir=None):
    '''
    by default, outdir will be same dir as infile
    assumes infile fast has already been trimmed to remove SSI
    '''
    aligner = config.get('bcfasta','tool')
    
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    if outdir is not None:
        dirname = outdir
    
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    logging.debug(f'handling {filepath}')
    
    logging.info('calc counts...')
    df = make_counts_df(config, infile)
    logging.info('save initial counts df...')
    of = os.path.join(dirname , f'{base}.counts.tsv')
    df.to_csv(of, sep='\t') 
        
    threshold = calculate_threshold(config, df)
    logging.info(f'threshold={threshold}')
    df = threshold_counts(config, df)
    df['sequence'] = df['sequence'].str[:32]
    spikedf, realdf, lonedf = split_spike_real_lone_barcodes(config, df)

    realdf.to_csv(os.path.join(dirname , f'{base}.real.raw.tsv'), sep='\t')
    lonedf.to_csv(os.path.join(dirname , f'{base}.lone.raw.tsv'), sep='\t')
    spikedf.to_csv(os.path.join(dirname , f'{base}.spike.raw.tsv'), sep='\t')

    
    logging.info('handling reals...')
    logging.info(f'real {len(realdf)} sequences, representing {realdf.counts.sum()} reads.')      
    of = os.path.join(dirname , f'{base}.real.seq.fasta')
    logging.debug(f'make fasta for {aligner} = {of}') 
    seqfasta = write_fasta_from_df(config, realdf, outfile=of)
    of = os.path.join(dirname , f'{base}.real.{aligner}')
    logging.info(f'running {aligner}...')
    afile = run_bowtie(config, seqfasta, of, tool=aligner )  
    logging.info(f'handle bowtie align file: {afile}')
    btdf = make_bowtie_df(afile)
    of = os.path.join(dirname , f'{base}.real.btdf.tsv')
    btdf.to_csv(of, sep='\t') 
    edgelist = edges_from_btdf(btdf)
    components = get_components(edgelist)
    crealdf = collapse_counts_df(realdf, components)
    logging.info(f'oldlen={len(realdf)}, collapsed={len(crealdf)}')
    
    
    # process spikeins
    logging.info(f'spike {len(spikedf)} sequences, representing {spikedf.counts.sum()} reads.')      
    of = os.path.join(dirname , f'{base}.spike.seq.fasta')
    logging.debug(f'make fasta for {aligner} = {of}') 
    seqfasta = write_fasta_from_df(config, spikedf, outfile=of)
    of = os.path.join(dirname , f'{base}.spike.{aligner}')
    logging.info(f'running {aligner}...')
    afile = run_bowtie(config, seqfasta, of, tool=aligner )  
    logging.info(f'handle bowtie align file: {afile}')
    btdf = make_bowtie_df(afile)
    of = os.path.join(dirname , f'{base}.spike.btdf.tsv')
    btdf.to_csv(of, sep='\t') 
    edgelist = edges_from_btdf(btdf)
    components = get_components(edgelist)
    cspikedf = collapse_counts_df(spikedf, components)
    logging.info(f'oldlen={len(spikedf)}, collapsed={len(cspikedf)}')    
    

    # process L1 controls
    logging.info(f'lone {len(lonedf)} sequences, representing {lonedf.counts.sum()} reads.')      
    of = os.path.join(dirname , f'{base}.lone.seq.fasta')
    logging.debug(f'make fasta for {aligner} = {of}') 
    seqfasta = write_fasta_from_df(config, lonedf, outfile=of)
    of = os.path.join(dirname , f'{base}.lone.{aligner}')
    logging.info(f'running {aligner}...')
    clonedf = None
    try:
        afile = run_bowtie(config, seqfasta, of, tool=aligner )  
        logging.info(f'handle bowtie align file: {afile}')
        btdf = make_bowtie_df(afile)
        of = os.path.join(dirname , f'{base}.lone.btdf.tsv')
        btdf.to_csv(of, sep='\t')     
        edgelist = edges_from_btdf(btdf)
        components = get_components(edgelist)
        clonedf = collapse_counts_df(lonedf, components)
        logging.info(f'oldlen={len(lonedf)}, collapsed={len(clonedf)}')
    
    except NonZeroReturnException:
        logging.warning(f'NonZeroReturn Exception. Probably no L1s found. ')
      
    
    #logging.info('handling reals...')
    #logging.info(f'lone {len(lonedf)} sequences, representing {lonedf.counts.sum()} reads.')      
    #of = os.path.join(dirname , f'{base}.lone.seq.fasta')
    #logging.debug(f'make fasta for {aligner} = {of}') 
    #seqfasta = write_fasta_from_df(config, lonedf, outfile=of)
    #of = os.path.join(dirname , f'{base}.lone.{aligner}')
    #logging.info(f'running {aligner}...')
    #try:
    #    afile = run_bowtie(config, seqfasta, of, tool=aligner )  
    #    logging.info(f'handle bowtie align file: {afile}')
    #    btdf = make_bowtie_df(afile)
    #    of = os.path.join(dirname , f'{base}.lone.btdf.tsv')
    #    btdf.to_csv(of, sep='\t')
    #except NonZeroReturnException:
    #    logging.warning(f'NonZeroReturn Exception. Probably no L1s found. ')

    crealdf.to_csv(os.path.join(dirname , f'{base}.real.tsv'), sep='\t')
    cspikedf.to_csv(os.path.join(dirname , f'{base}.spike.tsv'), sep='\t')
    if clonedf is not None:
        clonedf.to_csv(os.path.join(dirname , f'{base}.lone.tsv'), sep='\t')

    return (crealdf, cspikedf)


def collapse_counts_df(countsdf, components):
    '''
    takes components consisting of indices
    determines sequence with largest count columns: 'sequence', 'counts'
    collapses all other member components to largest. 
    '''
    # list of lists to collect values..
    logging.info(f'collapsing countsdf len={len(countsdf)} w/ {len(components)} components.')
    lol = []
    for component in components:
        #logging.debug(f'component={component}')
        cdf = countsdf.iloc[component].reset_index(drop=True)
        maxid = cdf.counts.idxmax()
        row = list(cdf.iloc[maxid])
        row[1] = cdf.counts.sum()
        lol.append(row)
    
    newdf = pd.DataFrame(data=lol, columns=['sequence','counts'])
    logging.info(f'collapsed df len={len(newdf)}')
    return newdf


def edges_from_btdf(btdf):
    readlist = btdf.name_read.values.tolist()
    alignlist = btdf.name_align.values.tolist()  
    edgelist = [ list(t) for t in zip(readlist, alignlist)]
    return edgelist

def get_components(edgelist):
    complist = []
    logging.info(f'getting connected components from edgelist len={len(edgelist)}')
    for g in trajan(from_edges(edgelist)):
        complist.append(g)
    logging.info(f'{len(complist)} components.')
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


def make_adjacency_df(bowtiedf):
    labels = np.unique(btdf[['name_read','name_align']])
    sdf = btdf.filter( ['name_read','name_align'], axis=1 )
    sdf['val'] = 1
    mdf = sdf.pivot(index = 'name_read', 
                    columns='name_align', 
                    values='val').reindex(columns=labels, index=labels, fill_value=0)
    mdf.fillna(0, inplace=True)
    return mdf


def make_counts_df(config, infile, bclabel=None):
    '''
    input fasta 
    ignore 'N' sequences. 
    '''   
    slist = []
    rcs = SeqIO.parse(infile, "fasta")
    handled = 0
    for sr in rcs:
        s = sr.seq
        if 'N' in sr:
            pass
        else:
            slist.append(str(s))
        handled += 1    
    logging.info(f"kept {len(slist)} non-'N' sequences out of {handled}")    
    df = pd.DataFrame(slist, columns=['sequence'] )
    ser = df.sequence.value_counts()
    df = pd.DataFrame()
    df['sequence'] = ser.index
    df['counts'] = ser.values
    if bclabel is not None:
        df['bc_label'] = bclabel
    logging.debug(f'counts df = \n{df}')
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


def threshold_counts(config, df, threshold=None):
    '''
    
    '''
    if threshold is None:
        count_threshold = int(config.get('bcfasta', 'countthreshold'))
    else:
        count_threshold = int(threshold)
    logging.info(f'thresh = {count_threshold}')    
    df = df[df.counts > count_threshold]
    return df


def split_spike_real_lone_barcodes(config, df):
    '''
    df has  sequence  counts
    should be length 32 ( 30 + YY ) Y= C or T
          
    '''
    #  df[df["col"].str.contains("this string")==False]
    sire = config.get('bcfasta', 'spikeinregex')
    realre = config.get('bcfasta','realregex')
    lonere = config.get('bcfasta', 'loneregex')
    
    logging.debug(f'before filtering: {len(df)}')
    
    # get spikeins
    logging.debug(f"spike-in regex = '{sire}' ")
    simap = df['sequence'].str.contains(sire, regex=True) == True
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
      

class BarCodeHandler(object):
    '''
    Basically implements fastx_barcode_handler.pl. 
    Matches against given barcode sequence, and
    writes out target fasta (minus barcode) to barcode-specific file.
    
    check end of line of sequence, length of barcode only.   
    
        
    '''
    def __init__(self, label, barcode, outdir, eol=True, max_mismatch=0):
        self.barcode = barcode
        self.label = label
        self.filename = os.path.abspath(f'{outdir}/{label}.fasta')
        self.eol = True
        self.max_mismatch = max_mismatch
        self.of = None
        self.dataframe = None
        if outdir is None:
            outdir = "."

    def do_match(self, id, seq ):
        '''
        test sequence exactly against this barcode.
        just test end of sequence the same length as barcode. 
        only write part of sequence up to barcode/SSI.
        EOL testing only now, add BOL...
         
        '''
        if self.of is None:
            logging.debug(f'open file for {self.label}')
            self.of = open(self.filename, 'w')
            self.of.write(f';sequences for  barcode {self.label}\n')
        
        r = False
        if self.max_mismatch == 0:
            if self.barcode == seq[-len(self.barcode):]:
            #if match_strings(self.barcode , seq[-len(self.barcode):], max_mismatch=self.max_mismatch):        
                id = str(id)
                sr = SeqRecord( seq[:-len(self.barcode)], id=id, name=id, description=id)
                try:
                    SeqIO.write([sr], self.of, 'fasta')
                except Exception as ex:
                    logging.warning(f'problem with {self.label}')
                    logging.warning(traceback.format_exc(None))            
                r = True
        else:
            #if self.barcode == seq[-len(self.barcode):]:
            if match_strings(self.barcode , seq[-len(self.barcode):], max_mismatch=self.max_mismatch):        
                id = str(id)
                sr = SeqRecord( seq[:-len(self.barcode)], id=id, name=id, description=id)
                try:
                    SeqIO.write([sr], self.of, 'fasta')
                except Exception as ex:
                    logging.warning(f'problem with {self.label}')
                    logging.warning(traceback.format_exc(None))            
                r = True
        return r
        
    def finalize(self):
        logging.debug(f'closing file for {self.label}')
        self.of.close()



    @classmethod
    def make_all_bch_counts(cls, config, bcolist, outdir=None):
        '''
        Makes counts dataframes for all barcode handler objects. 
        Saves 
        '''
        alldf = pd.DataFrame({'bc_label': pd.Series(dtype='str'),
                   'sequence': pd.Series(dtype='str'),
                   'counts': pd.Series(dtype='int')})
        for bch in bcolist:
            bcfile = bch.ofname
            filepath = os.path.abspath(bcfile)    
            dirname = os.path.dirname(filepath)
            if outdir is not None:
                dirname = os.path.abspath(outdir)
            filename = os.path.basename(filepath)
            (base, ext) = os.path.splitext(filename)
            cdf = make_counts_df(config, bcfile)
            cdf['logcounts'] = np.log(cdf.counts)
            cdf['bc_label'] = bch.label       
            bch.dataframe = cdf
            of = os.path.join(dirname , f'{base}.counts.tsv')
            cdf.to_csv(of, sep='\t')
            alldf = pd.concat([alldf, cdf], ignore_index = True, copy=False)
        return alldf

    @classmethod
    def merge_counts(cls, config, dflist, outfile=None):
        '''
        Makes counts dataframe  objects. 
        Saves 
        '''
        alldf = pd.DataFrame({'bc_label': pd.Series(dtype='str'),
                   'sequence': pd.Series(dtype='str'),
                   'counts': pd.Series(dtype='int')})
        for df in dflist:
            alldf = pd.concat([alldf, df], ignore_index = True, copy=False)
        return alldf


def load_barcodes(config, bcfile, labels = None, outdir=None, eol=True, max_mismatch=0):
    '''
     labellist is a python list of ids. barcode labels will be checked against
     BC<id> from labellist. 
    '''
    codelist = [ f'BC{x}' for x in labels]
    bclist = []
    with open(bcfile) as fh:
        while True:
            ln = fh.readline()
            if len(ln) < 2:
                break
            (label, bcseq) = ln.split()
            if labels is None or label in codelist:
                bch = BarCodeHandler(label, bcseq, outdir, eol, max_mismatch)
                bclist.append(bch)
    logging.debug(f'made list of {len(bclist)} barcode handlers.')
    return bclist


def load_sample_info(config, file_name):    
    #['Tube # by user', 'Our Tube #', 'Sample names provided by user',
    #   'Site information', 'RT primers for MAPseq', 'Brain ', 'Column#']
    sample_columns = ['usertube', 'ourtube','samplename','siteinfo','rtprimer','brain','col_num'] 
    
    
    dfs = pd.read_excel(file_name, sheet_name=None)        
    sidf = dfs['Sample information']
    # set columns names to first row
    sidf.columns = sidf.iloc[0]
    # drop old row of column names. 
    sidf = sidf.iloc[pd.RangeIndex(len(sidf)).drop(0)]
    # drop all columsn beyond 7
    sidf = sidf.iloc[:,:7]
    sidf.columns = sample_columns
    logging.debug(f'created reduced sample info df: {sidf}')

    return sidf


def check_output(bclist):
    '''
    check to see if valid output exists for all input barcode handlers, in 
    order to facility short-circuiting. 
    '''
    output_exists = True
    missing = []
    for bch in bclist:
        if os.path.exists(bch.filename) and ( os.path.getsize(bch.filename) >= 1 ):
            logging.debug(f'Non-empty BC{bch.label}.fasta exists.')
        else:
            logging.info(f"{bch.filename} doesn't exist. output_exists=False")
            missing.append(bch.label)
            output_exists = False
    if output_exists == False:
        logging.debug(f'missing BC labels: {missing}')
    else:
        logging.info('all output exists.')
    return output_exists
   

def process_fastq_pair(config, read1file, read2file, bclist, outdir, force=False):

    # if all the output files for bclist exist, don't recalc unless force=True. 
    if outdir is None:
        outdir = "."
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
    
        if read1file.endswith('.gz'):
             read1f = gzip.open(read1file, "rt")
        if read2file.endswith('.gz'):
             read2f = gzip.open(read2file, "rt")         
            
        recs1 = SeqIO.parse(read1f, "fastq")
        recs2 = SeqIO.parse(read2f, "fastq")
    
        seqshandled = 0
        unmatched = 0
        didmatch = 0
        
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
                    logging.info(f'handled {seqshandled} reads. matched={didmatch} unmatched={unmatched}')
            
            except StopIteration as e:
                logging.debug(f'iteration stopped')
                break
        
        umf.close()
        pf.close()
        for bch in bclist:
            bch.finalize()    
        # close possible gzip filehandles??
        max_mismatch = bclist[0].max_mismatch
        logging.info(f'handled {seqshandled} sequences. {didmatch} matched. {unmatched} unmatched')
    else:
        logging.info('all output exists and force=False. Not recalculating.')


def make_summaries(config, bcolist):
    # make jobset
    pass


def find_connected_components():
    pass


#%find connected components
#graph=[];
#for i=1:length(bcnfilt)
#    %find the connected graph components
#    [graph(i).S, graph(i).G]= graphconncomp( clustermatrix1(i).C, 'Directed', 'false'); 
#end
#% save('graph1.mat','graph'); 
# https://stackoverflow.com/questions/53277121/dulmage-mendelsohn-matrix-decomposition-in-python

#
# /grid/zador/data_nlsas_norepl/MAPseq/processing/MATLAB common/MATLAB common/charlienash-nricp-5d1cb79/dependencies/geom3d-2017.12.01/geom3d/meshes3d
#
#function [S,C] = conncomp(G)
#% CONNCOMP Drop in replacement for graphconncomp.m from the bioinformatics
#% toobox. G is an n by n adjacency matrix, then this identifies the S
#% connected components C. This is also an order of magnitude faster.
#%
#% [S,C] = conncomp(G)
#%
#% Inputs:
#%   G  n by n adjacency matrix
#% Outputs:
#%   S  scalar number of connected components
#%   C

#% Transpose to match graphconncomp
#G = G';

#[p,~,r] = dmperm(G+speye(size(G)));
#S = numel(r)-1;
#C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
#C(p) = C;
#end

#%collapse barcodes to most abundant member of the connected graph component
#  
#for i=1:length(bcnfilt)
#    x=1:graph(i).S;
#    [tf,loc]=ismember(x,graph(i).G,'R2012a');
#    collapsedreads=data(i).reads(loc,:);
#    collapsedcounts=accumarray(graph(i).G',data(i).counts);%'
#    [corrected(i).counts2u,ix]=sort(collapsedcounts,'descend');
#    corrected(i).reads2u=collapsedreads(ix,:);
#end




