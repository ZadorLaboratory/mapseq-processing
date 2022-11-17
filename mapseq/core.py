import logging
import os
import sys

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from cshlwork.utils import dataframe_to_seqlist

def run_bowtie(config, infile, outfile):
    '''
    bowtie-build -q BC1_seq.fasta indexes/BC1_bt 
    bowtie -v 3 -p 10 -f --best -a indexes/BC1_bt BC1_seq.fasta BC1_bowtie.txt
    
    
    '''
    logging.info(f'running bowtie on {infile} -> {outfile}')
    r = "bowtiefile"
    return r

def process_bcfasta(config, infile):
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    logging.debug(f'handling {filepath}')
    
    df = make_counts_df(config, infile)
    df = do_threshold(config, df)
    df = remove_spikeins(config,df)
    of = os.path.join(dirname , f'{base}_seq.fasta')
    logging.debug(f'fasta for bowtie = {of}') 
    seqfasta = write_fasta_for_bowtie(config, df, outfile=of)
    of = os.path.join(dirname , f'{base}_seq.bowtie')
    a = run_bowtie(config, seqfasta, of )
    return df

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


def remove_spikeins(config, df):
    #  df[df["col"].str.contains("this string")==False]
    si = config.get('bcfasta', 'spikein')
    logging.debug(f'before filter {len(df)}')
    df = df[df['sequence'].str.contains(si) == False ]
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









