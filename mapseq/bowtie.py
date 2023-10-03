import logging
import os
import sys
import traceback

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

import numpy as np
import pandas as pd

from cshlwork.utils import run_command_shell, NonZeroReturnException, setup_logging

#
# Possibly set up to do pipline internally?
# See:  https://www.youtube.com/watch?v=_cTbADrGLCQ  
#
#    ConfigParser entries:
#
#    [bowtie]
#    threads = 10
#
#    [bowtie2]
#    threads = 10
#

BOWTIE_1_COLS=['name_read', 'strand','name_align','offset','seq','quals','ceil','mm_desc']

def run_bowtie(config, infile, outfile, tool='bowtie'):
    '''
    bowtie-build -q BC1.seq.fasta indexes/BC1.bt 
    bowtie -v 3 -p 10 -f --best -a indexes/BC1.bt BC1_seq.fasta BC1.bt.algn
   
    '''
    logging.debug(f'running allxall bowtie on {infile} -> {outfile}')
     
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    logging.debug(f'handling {filepath}')
    
    idxdir = os.path.abspath(f'{dirname}/indexes')
    os.makedirs(idxdir, exist_ok = True )
    idxpfx = f'{idxdir}/{base}'

    # don't run seconds if first fails. 
    build_ok = False
    

    if tool == 'bowtie':
        cmd = ['bowtie-build',
               #'-q',
               infile,
               idxpfx, 
               ]
    elif tool == 'bowtie2':
        cmd = ['bowtie2-build',
           #'-q',
           infile,
           idxpfx, 
           ]
    logging.debug(f'running {cmd}')
    try:
        run_command_shell(cmd)
        build_ok = True
        
    except NonZeroReturnException as nzre:
        logging.error(f'problem with infile {infile}')
        logging.error(traceback.format_exc(None))
        raise     

    logging.debug(f'bowtie-build done.')

    threads = config.get(tool, 'threads')
    max_mismatch = config.get(tool, 'max_mismatch')
    if tool == 'bowtie':    
        cmd = ['bowtie',
               '-v', max_mismatch,
               '-p', threads, # # threads
               '-f',      # -f query input files are (multi-)FASTA .fa/.mfa
               '--best',
               '-a',      # -a/--all report all alignments; very slow, MAPQ not meaningful
               idxpfx,
               infile,
               outfile
               ]
    elif tool == 'bowtie2':
        
        cmd = ['bowtie2',
               #'-N', '3',
               '-p',threads,   # # threads
               '-f',       # -f query input files are (multi-)FASTA .fa/.mfa
               #'--best',
               '--all',   # -a/--all report all alignments; very slow, MAPQ not meaningful
               '-x', idxpfx,
               '-S', outfile, 
               infile,
               ]
    
    
    logging.debug(f'running bowtie...')
    try:
        run_command_shell(cmd)
    
    except NonZeroReturnException as nzre:
        logging.error(f'problem with infile {infile}')
        logging.error(traceback.format_exc(None))
        raise         
    logging.debug(f'bowtie done.')
    return outfile

def make_bowtie_df(infile):
    '''
    parse standard bowtie output and create pandas dataframe with standard columns
    
    '''
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)    
    bdf = pd.read_csv(filepath, sep='\t',header=None, names=BOWTIE_1_COLS)
    return bdf


def make_adjacency_df(bowtiedf):
    '''
    consume custom bowtie dataframe (from all x all alignment) and
    create adjacency matrix dataframe
    '''
    labels = np.unique(btdf[['name_read','name_align']])
    sdf = btdf.filter( ['name_read','name_align'], axis=1 )
    sdf['val'] = 1
    mdf = sdf.pivot(index = 'name_read', 
                    columns='name_align', 
                    values='val').reindex(columns=labels, index=labels, fill_value=0)
    mdf.fillna(0, inplace=True)
    return mdf


    
def make_bowtie2_df(filepath):
    logging.warning('not implemented')