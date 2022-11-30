import logging
import os
import sys
import traceback

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from cshlwork.utils import dataframe_to_seqlist, run_command_shell, NonZeroReturnException, setup_logging

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
    logging.info('threshold...')
    df = do_threshold(config, df)
    logging.info('remove spikeins...')
    df = remove_spikeins(config,df)
    of = os.path.join(dirname , f'{base}.seq.fasta')
    logging.debug(f'fasta for bowtie = {of}') 
    logging.info('make fasta seqfile for bowtie...')
    seqfasta = write_fasta_for_bowtie(config, df, outfile=of)
    of = os.path.join(dirname , f'{base}.seq.bowtie')
    logging.info('run bowtie...')
    afile = run_bowtie(config, seqfasta, of )
    logging.info(f'handle bowtie align file: {afile}')
    df = handle_alignment(afile)
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


def run_bowtie(config, infile, outfile):
    '''
    bowtie-build -q BC1.seq.fasta indexes/BC1.bt 
    bowtie -v 3 -p 10 -f --best -a indexes/BC1.bt BC1_seq.fasta BC1.bt.algn
    
    
    '''
    logging.info(f'running allxall bowtie on {infile} -> {outfile}')
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    logging.debug(f'handling {filepath}')
    
    idxdir = os.path.abspath(f'{dirname}/indexes')
    os.makedirs(idxdir, exist_ok = True )
    idxpfx = f'{idxdir}/{base}'
    cmd = ['bowtie2-build',
           #'-q',
           infile,
           idxpfx, 
           ]
    logging.debug(f'running bowtie-build...')
    try:
        run_command_shell(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with infile {infile}')
        logging.error(traceback.format_exc(None))
        raise     

    logging.info(f'bowtie-build done.')
    #  bowtie -v 3 -p 10 -f --best -a indexes/BC1.bt BC1_seq.fasta BC1.bt.algn
    #
    #  bowtie2 -v 3 -p 10 -f --best -a 
    #     /Users/jhover/project/mapseq/M205testout/indexes/BC1.seq 
    #     /Users/jhover/project/mapseq/M205testout/BC1.seq.fasta 
    #     /Users/jhover/project/mapseq/M205testout/BC1.seq.bowtie
    #
    #
       
    cmd1 = ['bowtie',
           '-v', '3',
           '-p','10', # # threads
           '-f',      # -f query input files are (multi-)FASTA .fa/.mfa
           '--best',
           '-a',      # -a/--all report all alignments; very slow, MAPQ not meaningful
           idxpfx,
           infile,
           outfile
           ]
    
    cmd2 = ['bowtie2',
           #'-N', '3',
           '-p','10',   # # threads
           '-f',       # -f query input files are (multi-)FASTA .fa/.mfa
           #'--best',
           '--all',   # -a/--all report all alignments; very slow, MAPQ not meaningful
           '-x', idxpfx,
           infile,
           outfile
           ]
    
    
    logging.debug(f'running bowtie...')
    try:
        run_command_shell(cmd2)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with infile {infile}')
        logging.error(traceback.format_exc(None))
        raise         
    logging.info(f'bowtie done.')
    return outfile

def handle_alignment(infile):
    pass







