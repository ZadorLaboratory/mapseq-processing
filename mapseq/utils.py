import glob
import gzip
import io
import itertools
import os
import logging
import pprint
import shutil
import subprocess
import sys
import tempfile
import threading
import traceback
import urllib

from configparser import ConfigParser

import datetime as dt

import numpy as np
import pandas as pd

from scipy import sparse, stats
from ftplib import FTP

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)

#
# Multiprocessing  using explicity command running. 
#            jstack = JobStack()
#            cmd = ['program','-a','arg1','-b','arg2','arg3']
#            jstack.addjob(cmd)
#            jset = JobSet(max_processes = threads, jobstack = jstack)
#            jset.runjobs()
#
#            will block until all jobs in jobstack are done, using <max_processes> jobrunners that
#            pull from the stack.
#

class JobSet(object):
    def __init__(self, max_processes, jobstack):
        self.max_processes = max_processes
        self.jobstack = jobstack
        self.threadlist = []
        
        for x in range(0, self.max_processes):
            jr = JobRunner(self.jobstack, label=f'{x}')
            self.threadlist.append(jr)
        logging.debug(f'made {len(self.threadlist)} jobrunners. ')


    def runjobs(self):
        logging.debug(f'starting {len(self.threadlist)} threads...')
        for th in self.threadlist:
            th.start()
            
        logging.debug(f'joining threads...')    
        for th in self.threadlist:
            th.join()

        logging.debug(f'all threads joined. returning...')


class JobStack(object):
    def __init__(self):
        self.stack = []

    def addjob(self, cmdlist):
        '''
        List is tokens appropriate for 
        e.g. cmd list :  [ '/usr/bin/x','-x','xarg','-y','yarg']
        '''
        self.stack.append(cmdlist)
    
    def pop(self):
        return self.stack.pop()


class JobRunner(threading.Thread):

    def __init__(self, jobstack, label=None):
        super(JobRunner, self).__init__()
        self.jobstack = jobstack
        self.label = label
        
    def run(self):
        while True:
            try:
                cmdlist = self.jobstack.pop()
                cmdstr = ' '.join(cmdlist)
                logging.info(f'[{self.label}] running {cmdstr}')
                logging.debug(f'[{self.label}] got command: {cmdlist}')
                run_command_shell(cmdlist)
                logging.debug(f'[{self.label}] completed command: {cmdlist}')
                logging.info(f'[{self.label}] completed {cmdstr} ')
            
            except NonZeroReturnException:
                logging.warning(f'[{self.label}] NonZeroReturn Exception job: {cmdlist}') 
            
            except IndexError:
                logging.info(f'[{self.label}] Command list empty. Ending...')
                break
                

class NonZeroReturnException(Exception):
    """
    Thrown when a command has non-zero return code. 
    """

def setup_logging(level):
    """ 
    Setup logging using e.g. level=logging.DEBUG
    """
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
    logging.basicConfig()
    logger = logging.getLogger()
    streamHandler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(FORMAT)
    streamHandler.setFormatter(formatter)
    logger.addHandler(streamHandler)
    logger.setLevel(level)



def findrepeats(s, n=2):
    '''
    Finds strings of repeated characters longer than <num>
    Typically used to find homopolymers in sequences for exclusion from results
    as repeats interfere with amplification/replication.
    '''
    rcount = 1
    last = None
    res = False
    for c in s:
        if last is None:
            last = c
        else:
            if c == last:
                rcount += 1
            else:
                last = c
                rcount = 1
        if rcount >= n:
            # short-circuit if found. 
            return True
        else:
            pass
    return res

def remove_base_repeats(df, col='sequence', n=7):
    '''
    removes rows where column string repeats C,G,A, or T <n> times or more.    
    '''
    startlen=len(df)
    logging.debug(f'removing n>{n} repeats from df len={startlen}')
    bases = ['C','G','A','T']
    for b in bases:
        pat = b*n
        logging.debug(f'searching for {pat}')
        df = df[ ~df[col].str.contains(pat)]
    
    endlen=(len(df))
    logging.debug(f'df without {n} repeats len={endlen}')
    df.reset_index(drop=True, inplace=True)
    return df

def remove_singletons(listoflists):
    '''
    Assumes input of a list of lists. 
    Return only lists with more than one element.
    
    E.g. used to remove single-element components from align_and_collapse
    
    '''
    logging.debug(f'len before = {len(listoflists)}')
    outlist = [ x for x in listoflists if len(x) > 1 ]
    logging.debug(f'len after = {len(outlist)}')
    return outlist
        


def has_base_repeats(seqstring, n=7):
    '''
    if string repeats C,G,A, or T <n> times or more.    
    '''
    bases = ['C','G','A','T']
    for b in bases:
        pat = b*n
        if pat in seqstring:
            return True
    return False

def flatten_list(listoflists):
    flat_list = []
    for row in listoflists:
        flat_list += row
    return flat_list


def fix_columns_float(df, columns):
    for col in columns:
        try:
            logging.debug(f'trying to fix col {col}')
            fixed = np.array(df[col], np.float64)
            df[col] = fixed
        except ValueError:
            logging.debug(f'invalid literal in {col}')
    # np cast sets nans to 0, change them back:
    df.replace(0, np.nan, inplace=True)
    return df

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

                   


def add_rowlist_column(rowlist, colval):
    """
    For use during dataframe construction. Adds col to list of rows with specified.
       
    """
    for row in rowlist:
        row.append(colval)
    return rowlist
    


def gzip_decompress(filename):
    """
    default for copyfileobj is 16384
    https://blogs.blumetech.com/blumetechs-tech-blog/2011/05/faster-python-file-copy.html

    """
    log = logging.getLogger('utils')
    if filename.endswith('.gz'):
        targetname = filename[:-3]
        bufferlength = 10 * 1024 * 1024  # 10 MB
        with gzip.open(filename, 'rb') as f_in:
            with open(targetname, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out, length=bufferlength)
    else:
        log.warn(
            f'tried to gunzip file without .gz extension {filename}. doing nothing.')

 
def run_command(cmd):
    """
    cmd should be standard list of tokens...  ['cmd','arg1','arg2'] with cmd on shell PATH.
    @return 
    
    
    """
    cmdstr = " ".join(cmd)
    logging.debug(f"command: {cmdstr} running...")
    start = dt.datetime.now()
    cp = subprocess.run( cmd, 
                    text=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT
                    )
    end = dt.datetime.now()
    elapsed =  end - start
    logging.debug(f"ran cmd='{cmdstr}' return={cp.returncode} {elapsed.seconds} seconds.")
    
    if cp.stderr is not None:
        logging.warn(f"got stderr: {cp.stderr}")
    if cp.stdout is not None:
        logging.debug(f"got stdout: {cp.stdout}")   
    if str(cp.returncode) == '0':
        #logging.debug(f'successfully ran {cmdstr}')
        logging.debug(f'got rc={cp.returncode} command= {cmdstr}')
    else:
        logging.warn(f'got rc={cp.returncode} command= {cmdstr}')
       
        #raise NonZeroReturnException(f'For cmd {cmdstr}')
    return cp

def run_command_shell(cmd):
    """
    maybe subprocess.run(" ".join(cmd), shell=True)
    cmd should be standard list of tokens...  ['cmd','arg1','arg2'] with cmd on shell PATH.
    
    """
    cmdstr = " ".join(cmd)
    logging.debug(f"running command: {cmdstr} ")
    start = dt.datetime.now()
    cp = subprocess.run(" ".join(cmd), 
                    shell=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT)

    end = dt.datetime.now()
    elapsed =  end - start
    logging.debug(f"ran cmd='{cmdstr}' return={cp.returncode} {elapsed.seconds} seconds.")
    
    if cp.stderr is not None:
        logging.warn(f"got stderr: {cp.stderr}")
        pass
    if cp.stdout is not None:
        #logging.debug(f"got stdout: {cp.stdout}")
        pass
    if str(cp.returncode) == '0':
        #logging.debug(f'successfully ran {cmdstr}')
        logging.debug(f'got rc={cp.returncode} command= {cmdstr}')
    else:
        logging.warn(f'got rc={cp.returncode} command= {cmdstr}')
        raise NonZeroReturnException(f'For cmd {cmdstr}')
    return cp



def string_modulo(instring, divisor):
    """
    Takes instring. Converts to bytes. Takes hex() value of bytes. converts to integer. 
    returns final integer % modbase
    
    """
    encoded = instring.encode('utf-8')
    hstring = encoded.hex()
    intval = int(hstring, 16)
    return intval % divisor


def modulo_filter(inlist, divisor, remainder):
    """
    Takes a list, returns list containing items in inlist that 
    have the given remainder modulo divisor. 
    """
    newlist = []
    for e in inlist:
        if string_modulo(e, divisor) == remainder:
            newlist.append(e)
    logging.debug(f'inlist len={len(inlist)}, {divisor} servers, {remainder} server idx. outlist len={len(newlist)}')
    return newlist
    


def get_configstr(cp):
    with io.StringIO() as ss:
        cp.write(ss)
        ss.seek(0)  # rewind
        return ss.read()

def format_config(cp):
    cdict = {section: dict(cp[section]) for section in cp.sections()}
    s = pprint.pformat(cdict, indent=4)
    return s

def get_mainbase(filepath):
    '''
    for any full or relative filename XXXXX.y.z.w.ext 
    returns just XXXXX (first portion before dot). 
    
    '''
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)
    base = base.split('.')[0]
    return base


def write_fasta_from_df(df, outfile=None):
    '''
    Assumes df has 'sequence' column
    
    '''
    logging.debug(f'writing {len(df)} sequence as fasta from DF...')
    idx = 0
    if outfile is not None:
        with open(outfile, 'w') as of:
            for s in df['sequence']:
                of.write(f'>{idx}\n{s}\n')
                idx += 1  
    else:
        logging.error('outfile is None, not implemented.')
    return outfile


def write_fasta_from_df_bioconda(df, outfile=None):
    '''
    Assumes df has 'sequence' column
    
    '''
    logging.debug(f'writing fasta from DF..')
    srlist = dataframe_to_seqlist(df)
    logging.debug(f'len srlist={len(srlist)}')
    if outfile is not None:
        SeqIO.write(srlist, outfile, 'fasta')
    else:
        logging.error('outfile is None, not implemented.')
    return outfile


def read_fasta_to_df_bioconda(infile, seqlen=None):
    '''
    input fasta 
    Optionally trim sequence to seqlen, returning a two-column dataframe 'sequence' 'tail'
    None means keep all. 
    '''   
    slist = []
    tlist = []
    rcs = SeqIO.parse(infile, "fasta")
    handled = 0
    df = None
    if seqlen is None:
        for sr in rcs:
            s = sr.seq
            slist.append(str(s))
            handled += 1    
        df = pd.DataFrame(slist, columns=['sequence'] )
    else:
        seqlen = int(seqlen)
        for sr in rcs:
            s = sr.seq[:seqlen]
            t = sr.seq[seqlen:]
            slist.append( [str(s),str(t)] )
            handled += 1
        df = pd.DataFrame(slist, columns=['sequence','tail'] )
    logging.debug(f"handled {handled}  sequences. df=\n{df}")    
    return df


def read_fasta_to_df(infile, seqlen=None):
    '''
    input fasta 
    Optionally trim sequence to seqlen, returning a two-column dataframe 'sequence' 'tail'
    None means keep all.
    
    use string[pyarrow] ?
    https://pythonspeed.com/articles/pandas-string-dtype-memory/
    
     
    '''
    handled_interval = 10000000   
    slist = []
    tlist = []
    handled = 0
    df = None
    if seqlen is None:
        # whole sequence goes in sequence, no tail
        with open(infile) as fa:
            for line in fa:
                if line.startswith(">"):
                    pass
                else:
                    slist.append(line[:-1])
                    handled += 1 
                    if handled % handled_interval == 0:
                        logging.info(f'handled={handled}')  
            logging.info(f'making df from {len(slist)} sequences...')    
            df = pd.DataFrame(slist, columns=['sequence'] )
            
    else:
        with open(infile) as fa:
            for line in fa:
                if line.startswith(">"):
                    pass
                else:
                    slist.append(line[:seqlen])
                    tlist.append(line[seqlen:-1])
                    handled += 1 
                    if handled % handled_interval == 0:
                        logging.info(f'handled={handled}') 
            logging.info(f'making df from {len(slist)} sequences and tails...')    
            ss = pd.Series(slist, dtype="string[pyarrow]")
            ts = pd.Series(tlist, dtype="string[pyarrow]")
            df = pd.DataFrame(columns=['sequence','tail'] )     
            df['sequence'] = ss
            df['tail'] = ts
    logging.debug(f"handled {handled}  sequences. df=\n{df}")    
    return df


def dataframe_to_seqlist(df, seqcol='sequence',idcol=None, desccols=None, sep=':'):
    '''
    reads df and produces python list of BioPython SeqRecord objects. 
    
    seqcol    column to be output as sequence
    desccols  list of column names to be included, in order after ">"
    sep       character to separate fields of idcols/desccols 
    '''
    logging.debug('converting dataframe to SeqRecord list...')
    srlist = []
    for index, row in df.iterrows():
        s = row[seqcol]
        seq = Seq(s)
        if idcol is not None:
            id = str(row[idcol])
        else:
            id = str(index)
        if desccols is not None:
            slist = []
            for coln in desccols:
                slist += str(row[coln])
                desc = sep.join(slist)
        else:
            desc = str('')
        sr = SeqRecord(seq, id = id, name=id, description=desc )
        srlist.append(sr)
    logging.debug(f'made list of {len(srlist)} SeqRecords')
    return srlist    
            

def convert_numeric(df):
    '''
    Convert all columns that can be successfully cast to numeric. 
    '''
    for c in df.columns:
        try:
            newc = pd.to_numeric(df[c])
            df[c] = newc
        except:
            logging.debug(f'column {c} not convertible to numeric.')
    # changes made in place, but also return so df = convert_numeric(df) works. 
    return df


def load_df(filepath):
    """
    Convenience method to load DF consistently across modules. 
    """
    logging.debug(f'loading {filepath}')
    filepath = os.path.expanduser(filepath)
    df = pd.read_csv(filepath, sep='\t', index_col=0, keep_default_na=False, dtype=str, comment="#")
    #df.fillna(value='', inplace=True)

    df = df.convert_dtypes(convert_integer=False)
    for col in df.columns:
        try:
            df[col] = df[col].astype('uint32')
        except ValueError:
            logging.debug(f'skipping column {col}')
    logging.debug(f'{df.dtypes}')
    return df


def merge_dfs( dflist):
    newdf = None
    for df in dflist:
        if newdf is None:
            newdf = df
        else:
            newdf = pd.concat([df, newdf], ignore_index=True, copy=False)
    logging.debug(f'merged {len(dflist)} dataframes newdf len={len(newdf)}')
    return newdf


def merge_tsvs( filelist):
    logging.debug(f'merge list {filelist}')
    newdf = None  
    for f in filelist:
        try:
            df = pd.read_csv(f, sep='\t', index_col=0)
        
            if newdf is None:
                newdf = df
            else:
                newdf = pd.concat([df, newdf], ignore_index=True, copy=False)
        except FileNotFoundError as fnfe:
            logging.warning(f'expected file {f} not found. Ignoring...')
    logging.debug(f'merged {len(filelist)} tsv files new df len={len(newdf)}')
    newdf.reset_index(drop=True, inplace=True)
    return newdf
    

def merge_write_df(newdf, filepath,  mode=0o644):
    """
    Reads existing, merges new, drops duplicates, writes to temp, renames temp. 
    """
    log = logging.getLogger('utils')
    log.debug(f'inbound new df:\n{newdf}')
    filepath = os.path.abspath( os.path.expanduser(filepath) )
    if os.path.isfile(filepath):
        df = load_df(filepath)
        log.debug(f'read df:\n{df}')
        df = pd.concat([df, newdf], ignore_index=True, copy=False)
        df.fillna(value='', inplace=True)
        df = df.astype('str', copy=False)
        log.debug(f'appended df:\n{df}')
    else:
        df = newdf
        df.fillna(value='', inplace=True)
        df = df.astype('str', copy=False)
    logging.debug(f"df length before dropping dupes is {len(df)}")
    df.drop_duplicates(inplace=True, ignore_index=True, keep='first')
    logging.debug(f"df length after dropping dupes is {len(df)}")
    df = df.reset_index(drop=True)
    rootpath = os.path.dirname(filepath)
    basename = os.path.basename(filepath)
    try:
        (tfd, tfname) = tempfile.mkstemp(suffix=None,
                                         prefix=f"{basename}.",
                                         dir=f"{rootpath}/",
                                         text=True)
        logging.debug(f"made temp {tfname}")
        df.to_csv(tfname, sep='\t')
        os.rename(tfname, filepath)
        os.chmod(filepath, mode)
        logging.debug(f"wrote df to {filepath}")

    except Exception as ex:
        logging.error(traceback.format_exc(None))
        raise ex


def write_config(config, filename, timestamp=True, datestring=None):
    '''
    writes config file to relevant name,
    if timestamp=True, puts date/time code dot-separated before extension. e.g.
    filename = /path/to/some.file.string.txt  ->  /path/to/some.file.string.202303081433.txt
    date is YEAR/MO/DY/HR/MI
    if datestring is not None, uses that timestamp
    
    '''
    filepath = os.path.abspath(filename)    
    dirname = os.path.dirname(filepath)
    basefilename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(basefilename) 
    
    if timestamp:
        if datestring is None:
            datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
        else:
            datestr = datestring
        filename = f'{dirname}/{base}.{datestr}{ext}'

    os.makedirs(dirname, exist_ok=True)
        
    with open(filename, 'w') as configfile:
        config.write(configfile)
    logging.debug(f'wrote current config to {filename}')
    return os.path.abspath(filename)
    


def write_df(newdf, filepath,  mode=0o644):
    """
    Writes df in standard format.  
    """
    logging.debug(f'inbound df:\n{newdf}')
    try:
        df = newdf.reset_index(drop=True)
        rootpath = os.path.dirname(filepath)
        basename = os.path.basename(filepath)
        df.to_csv(filepath, sep='\t')
        os.chmod(filepath, mode)
        logging.debug(f"wrote df to {filepath}")

    except Exception as ex:
        logging.error(traceback.format_exc(None))
        raise ex

def write_tsv(df, outfile=None):
    if outfile is None:       
        outfile = sys.stdout
    logging.debug(f'writing {len(df)} lines output to {outfile}')      
    df.to_csv(outfile, sep='\t')

def readlist(filepath):
    '''
    Assumes file is a list of strings, one per line. 
    Ignores lines beginning with a has '#'
    Ignores characters in a line afeter a '#'
    '''

    if filepath is not None:
        logging.debug(f'reading file: {filepath}')
        flist = []
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if len(line) > 0:
                        idx = line.find('#')
                        if idx == -1:
                            flist.append(line.strip())
                        elif idx > 0:
                            flist.append(line[:idx].strip())
                    else:
                        pass   # empty line
                        
            logging.debug(f'got list with {len(flist)} items.')
            return flist
        except:
            return []
    else:
        logging.debug('no file. return [].')
        return []



def writelist(filepath, dlist, mode=0o644):
    logging.debug(f"writing list length={len(dlist)} to file='{filepath}'")
    rootpath = os.path.dirname(filepath)
    basename = os.path.basename(filepath)
    try:
        (tfd, tfname) = tempfile.mkstemp(suffix=None,
                                         prefix=f"{basename}.",
                                         dir=f"{rootpath}/",
                                         text=True)
        logging.debug(f"made temp {tfname}")
        with os.fdopen(tfd, 'w') as f:
            nlines = 0
            for item in dlist:
                f.write(f"{item}\n")
                nlines += 1
        os.rename(tfname, filepath)
        os.chmod(filepath, mode)
        logging.debug(f"wrote {nlines} to {filepath}")
    except Exception as ex:
        logging.error(traceback.format_exc(None))

    finally:
        pass



