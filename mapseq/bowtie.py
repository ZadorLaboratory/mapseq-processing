#!/usr/bin/env python
# 

import argparse
import logging
import os
import sys
import traceback

import numpy as np
import pandas as pd

from collections import defaultdict

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

# Stop using cshlwork utils. Copy any needed directly into mapseq package. 
#from cshlwork.utils import *

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

#
#  first row is mandatory fields. 
#   second row are optional expected fields. 
BOWTIE_2_COLS=['name_read', 'flagsum', 'name_align','offset','qual', 'cigar', 'mate', 'mate_offset', 'fraglen', 'seq', 'quals', 
                  'score', 'next', 'n_amb', 'n_mismatch', 'n_gap', 'n_gapext', 'distance', 'md', 'yt' ]
BOWTIE_2_INT_COLS= ['flagsum', 'offset', 'qual', 'mate_offset', 'fraglen', 'score', 'n_amb', 'n_mismatch', 'n_gap', 'n_gapext', 'distance', ]


BOWTIE_OPT_COLS= ['score', 'next', 'n_amb', 'n_mismatch', 'n_gap', 'n_gapext', 'distance', 'md', 'yt' ]
OPT_MAP = { 'score'     : 'AS',
            'next'      : 'XS',
            'n_amb'     : 'XN',
            'n_mismatch': 'XM',
            'n_gap'     : 'XO',
            'n_gapext'  : 'XG',
            'distance'  : 'NM', 
            'md'        : 'MD',
            'yt'        : 'YT'    
            }


def run_bowtie(config, infile, outfile, tool='bowtie'):
    '''
    bowtie-build -q BC1.seq.fasta indexes/BC1.bt 
    bowtie -v 3 -p 10 -f --best -a indexes/BC1.bt BC1_seq.fasta BC1.bt.algn
   
    
    mkdir indexes.bt2
    bowtie2-build -f BC4.seq.fasta indexes/BC4
    bowtie2 -p 10 -f --end-to-end --fast --all -x indexes/BC4 -S BC4.bowtie2 BC4.real.seq.fasta  
   
    '''
    logging.debug(f'running allxall {tool} on {infile} -> {outfile}')
     
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
           '-f',
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

    if tool == 'bowtie':    
        max_mismatch = config.get(tool, 'max_mismatch')
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



def make_bowtie_df(infile, max_mismatch=3, ignore_self=False):
    with open(infile) as f:
        line=f.readline()
    if line.startswith('@HD'):
        logging.debug('Detected bowtie2 input.')
        df = make_bowtie2_df(infile)
        logging.debug(f'df before max_mismatch =< {max_mismatch}')
        df = df[df['n_mismatch'] <= max_mismatch]
        # alignments to *other* sequences only
        if ignore_self:
            df = df[df['n_mismatch'] > 0]
        logging.debug(f'df after max_mismatch < {max_mismatch} =\n{df}')
    else:
        logging.debug('Detected bowtie1 input.')
        df = make_bowtie1_df(infile)
        df['n_mismatch'] = max_mismatch
    return df


def make_bowtie1_df( infile):
    '''
    parse standard bowtie output and create pandas dataframe with standard columns
    
    '''
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)    
    bdf = pd.read_csv(filepath, sep='\t',header=None, names=BOWTIE_1_COLS)
    return bdf


def make_bowtie2_df(infile):
    '''
    bt1: name_read  strand  name_align  offset  seq  quals   ceil    mm_desc
    
    bt2: name_read         name_align  offset  seq  quals   
            algn_score  next_score  n_amb  n_mm  n_gaps n_gapext  distance      
    
    parse bowtie2 output format.
    produce dataframe.     
    '''
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)
    outfile = os.path.join(dirname, f'{base}.btdf')
    logging.debug(f'handling base={base}  ->  {outfile}')    

    try:
        logging.debug(f" attempting to open '{infile}'")
        filehandle = open(infile, 'r')
    except FileNotFoundError:
        logging.error(f"No such file {infile}")
        raise                
    
    current = 0
    sumreport = 1
    suminterval = 1000000
    repthresh = sumreport * suminterval
    
    # list of lists to hold data
    lol = []
    
    #
    def def_value(): 
        return "None"
    
    # BOWTIE_2_COLS=['name_read', 'flagsum', 'name_align','offset', 'qual', 'cigar', 'mate', 'mate_offset', 'fraglen', 'seq', 'quals', 
    # 'score', 'next_score', 'n_amb', 'n_mismatch', 'n_gaps', 'n_gapext', 'distance','md','yt' ]
    try:
        while True:
            line = filehandle.readline()
            if line == '':
                break
            if line.startswith("@HD"):
                pass
            elif line.startswith("@SQ"):
                pass                   
            elif line.startswith("@PG"):
                pass
            else:
                allfields = [x.strip() for x in line.split('\t')]
                mands = allfields[0:11]
                flist = mands
                optfields = allfields[11:]
                #logging.debug(f'allfields length={len(allfields)} mands len={len(mands)} opts len={len(optfields)}')
                #logging.debug(f'mands=\n{mands}')
                #logging.debug(f'opts=\n{optfields}')
                optdict = defaultdict(def_value)
                for of in optfields:
                    ofields = of.split(':')
                    key = ofields[0]
                    val = ofields[2]
                    optdict[key] = val 
                    #logging.debug(f'handling optkey {key} = {val}')
                for colname in BOWTIE_OPT_COLS:
                    mapid = OPT_MAP[colname]
                    val = optdict[mapid]
                    #logging.debug(f'for colname={colname} map={mapid} val={val}')
                    flist.append(val)
                
                lol.append(flist)
                current += 1                        
                if current >= repthresh:
                    logging.info(f"Processed {current} entries... ")
                    sumreport +=1
                    repthresh = sumreport * suminterval
    except Exception as e:
        traceback.print_exc(file=sys.stdout)                
    
    if filehandle is not None:
        filehandle.close()          

    df = pd.DataFrame(data=lol, columns=BOWTIE_2_COLS )
    df = fix_columns_int(df, BOWTIE_2_INT_COLS)
    return df


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

    



