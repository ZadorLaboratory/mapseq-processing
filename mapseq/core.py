import gzip
import logging
import os
import sys
import traceback

from configparser import ConfigParser

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from cshlwork.utils import dataframe_to_seqlist, run_command_shell, NonZeroReturnException, setup_logging
from alignment.bowtie import run_bowtie, make_bowtie_df


def get_default_config():
    dc = os.path.expanduser('~/git/mapseq-processing/etc/mapseq.conf')
    cp = ConfigParser()
    cp.read(dc)
    return cp


def process_bcfasta(config, infile, outdir=None):
    '''
    by default, outdir will be same dir as infile
    '''
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    logging.debug(f'handling {filepath}')
    logging.info('calc counts...')
    df = make_counts_df(config, infile)
    logging.info('save initial counts df...')
    of = os.path.join(dirname , f'{base}.counts.tsv')
    df.to_csv(of, sep='\t') 
    
    logging.info('threshold...')
    df = do_threshold(config, df)
    bctool = config.get('bcfasta','tool')
      
    
    # process regular barcodes
    logging.info('remove spikeins...')
    bcdf = get_barcodes(config, df)
    of = os.path.join(dirname , f'{base}.bc.seq.fasta')
    logging.debug(f'fasta for {bctool} = {of}') 
    logging.info(f'make fasta seqfile for {bctool}...')
    seqfasta = write_fasta_for_bowtie(config, bcdf, outfile=of)
    of = os.path.join(dirname , f'{base}.bc.seq.{bctool}')
    logging.info(f'running {bctool}...')
    afile = run_bowtie(config, seqfasta, of, tool=bctool )
    logging.info(f'handle bowtie align file: {afile}')
    btdf = make_bowtie_df(afile)
    
    # make sparse matrix 
    btmdf = matrix_df_from_btdf(btdf)


    return btmdf


    
    # process spikeins
    #logging.info('getting spike ins...')
    #sidf = get_spikeins(config,df)
    #of = os.path.join(dirname , f'{base}.si.seq.fasta')
    #logging.debug(f'fasta for {bctool} = {of}') 
    #logging.info(f'make fasta seqfile for {bctool}...')
    #seqfasta = write_fasta_for_bowtie(config, sidf, outfile=of)
    #of = os.path.join(dirname , f'{base}.si.seq.{bctool}')
    #logging.info(f'running {bctool}...')
    #afile = run_bowtie(config, seqfasta, of, tool=bctool  )
    #logging.info(f'handle bowtie align file: {afile}')
    
    
    # process L1 spikeins
    #logging.info('getting L1 spike ins...')
    #l1df = get_L1barcodes(config,df)
    #of = os.path.join(dirname , f'{base}.l1.seq.fasta')
    #logging.debug(f'fasta for {bctool} = {of}') 
    #logging.info(f'make fasta seqfile for {bctool}...')
    #seqfasta = write_fasta_for_bowtie(config, l1df, outfile=of)
    #of = os.path.join(dirname , f'{base}.l1.seq.{bctool}')
    #logging.info(f'running {bctool}...')
    #afile = run_bowtie(config, seqfasta, of, tool=bctool  )
    #logging.info(f'handle bowtie align file: {afile}')


def make_counts_df(config, infile):
    stt = int(config.get('bcfasta', 'start'))
    end = int(config.get('bcfasta', 'end')) 
    
    slist = []
    
    rcs = SeqIO.parse(infile, "fasta")
    handled = 0
    for sr in rcs:
        s = sr.seq
        if 'N' in sr:
            pass
        else:
            logging.debug(f'{s}')
            s = s[stt:end]
            logging.debug(f'{s}')
            slist.append(str(s))
        handled += 1    
    logging.debug(f"kept {len(slist)} non-'N' sequences out of {handled}")
    
    df = pd.DataFrame(slist, columns=['sequence'] )
    ser = df.sequence.value_counts()
    df = pd.DataFrame()
    df['sequence'] = ser.index
    df['counts'] = ser.values
    logging.debug(f'counts df = \n{df}')
    return df


def do_threshold(config, df):
    cthr = int(config.get('bcfasta', 'countthreshold'))
    spend = int(config.get('bcfasta', 'spikend'))
    logging.debug(f'thresh = {cthr}')
    df = df[df.counts > cthr]
    df['sequence'] = df.sequence.str[:spend]
    #df.drop(['nsequence'],axis=1, inplace=True)
    return df


def get_barcodes(config, df):
    # 
    # contains    '[TC][TC]$'
    #  df[df["col"].str.contains("this string")==False]
    si = config.get('bcfasta', 'spikein')
    logging.debug(f'before filter {len(df)}')
    df = df[df['sequence'].str.contains(si, regex=False) == False ]
    logging.debug(f'after first filter {len(df)}')    
    pat = '[TC][TC]$'
    df = df[df['sequence'].str.contains(pat) == True ]
    logging.debug(f'after second filter {len(df)}')    
    return df

def get_L1barcodes(config, df):
    #
    # contains   '[AG][AG]$'
    #  df[df["col"].str.contains("this string")==False]
    si = config.get('bcfasta', 'spikein')
    logging.debug(f'before filter {len(df)}')
    df = df[df['sequence'].str.contains(si, regex=False) == False ]
    logging.debug(f'after first filter {len(df)}')    
    pat = '[AG][AG]$'
    df = df[df['sequence'].str.contains(pat) == True ]
    logging.debug(f'after second filter {len(df)}')      
    return df


def get_spikeins(config, df):
    #  df[df["col"].str.contains("this string")==False]
    si = config.get('bcfasta', 'spikein')
    logging.debug(f'before filter {len(df)}')
    df = df[df['sequence'].str.contains(si, regex=False) == True ]
    logging.debug(f'after filter {len(df)}')    
    return df


def write_fasta_for_bowtie(config, df, outfile=None):
    logging.debug(f'creating bowtie input')
    srlist = dataframe_to_seqlist(df)
    logging.debug(f'len srlist={len(srlist)}')
    if outfile is not None:
        SeqIO.write(srlist, outfile, 'fasta')
    else:
        logging.error('outfile is None, not implemented.')
    return outfile


def matrix_df_from_btdf(df):
    '''
    takes bowtie read df
    emits boolean adjacency matrix of 'name_read','name_align'
    
    '''
    labels = np.unique(df[['name_read','name_align']])
    sdf = df.filter(['name_read','name_align'], axis=1)
    sdf['val'] = True
    mdf = sdf.pivot(index = 'name_read', columns='name_align', values='val').reindex(columns=labels, index=labels, fill_value=False)
    return mdf

class BarCodeHandler(object):
    
    def __init__(self, label, barcode, outdir):
        self.barcode = barcode
        self.label = label
        
        if outdir is None:
            outdir = "."
        self.ofname = os.path.abspath(f'{outdir}/{label}.fasta')
        self.of = open(self.ofname, 'w')
        logging.debug(f'open file for {self.label}')
        
    def do_match(self, id, seq):
        r = False
        if self.barcode in seq:
            id = str(id)
            sr = SeqRecord( seq, id=id, name=id, description=id)
            try:
                SeqIO.write([sr], self.of, 'fasta')
            except:
                logging.warning(f'problem with {self.label}...')
            
            r = True
        return r

    def finalize(self):
        logging.debug(f'closing file for {self.label}')
        self.of.close()
            

def load_barcodes(config, bcfile, labels = None, outdir=None):
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
                bch = BarCodeHandler(label, bcseq, outdir)
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


def process_fastq_pair(config, read1file, read2file, bclist, outdir):

    if outdir is None:
        outdir = "."
    outfile = os.path.abspath(f'{outdir}/unmatched.fasta')
    
    umf = open(outfile, 'w')

    r1s = int(config.get('fastq','r1start'))
    r1e = int(config.get('fastq','r1end'))
    r2s = int(config.get('fastq','r2start'))
    r2e = int(config.get('fastq','r2end'))   

    if read1file.endswith('.gz'):
         read1f = gzip.open(read1file, "rt")
    if read2file.endswith('.gz'):
         read2f = gzip.open(read1file, "rt")         
        
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
            matched = False
            for bch in bclist:
                r = bch.do_match(seqshandled, fullread)
                if r:
                    didmatch += 1
                    if didmatch % 10000 == 0:
                        logging.debug(f'match {didmatch}: found bc {bch.label} in {fullread}!!')
                    matched = True
                    break
            if not matched:
                unmatched += 1
                if unmatched % 10000 == 0:
                    logging.debug(f'{unmatched} unmatched so far.')
                id = str(seqshandled)
                sr = SeqRecord( fullread, id=id, name=id, description=id)
                SeqIO.write([sr], umf, 'fasta')
            
            seqshandled += 1
            if seqshandled % 500000 == 0: 
                logging.info(f'handled {seqshandled} reads. matched={didmatch} unmatched={unmatched}')
        
        except StopIteration as e:
            logging.debug(f'iteration stopped')
            break
    
    umf.close()
    for bch in bclist:
        bch.finalize()    
    # close possible gzip filehandles??
    logging.info(f'handled {seqshandled} read sequences. {unmatched} unmatched')


def make_summaries(config, bcolist):
    # make jobset
    pass






