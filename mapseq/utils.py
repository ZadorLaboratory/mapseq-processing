import glob
import gzip
import itertools
import os
import logging
import shutil
import subprocess
import sys
import tempfile
import threading
import traceback
import urllib

from configparser import ConfigParser
import datetime as dt

import io
import numpy as np
import pandas as pd

from scipy import sparse, stats
from ftplib import FTP

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)


class JobSet(object):
    def __init__(self, max_processes, jobstack):
        self.max_processes = max_processes
        self.jobstack = jobstack
        self.threadlist = []
        
        for x in range(0,self.max_processes):
            jr = JobRunner(self.jobstack)
            self.threadlist.append(jr)
        logging.debug(f'made {len(self.threadlist)} jobrunners. ')


    def runjobs(self):
        logging.debug(f'starting threads...')
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

    def __init__(self, jobstack):
        super(JobRunner, self).__init__()
        self.jobstack = jobstack
        

    def run(self):
        while True:
            try:
                cmdlist = self.jobstack.pop()
                logging.debug(f'got command: {cmdlist}')
                run_command_shell(cmdlist)
                logging.debug(f'completed command: {cmdlist}')
            except NonZeroReturnException:
                logging.warning(f'NonZeroReturn Exception job: {cmdlist}') 



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


def remove_parentpath(parent, fullpath):
    '''
      remove parent from fullpath
    '''
    logging.debug(f'fullpath ={fullpath}')
    compath = os.path.commonpath([parent, fullpath])
    logging.debug(f'commonpath ={compath}')
    comlist = splitall(compath)
    logging.debug(f'comlist ={comlist}')
    logging.debug(f'len comlist = {len(comlist)}')
    comlength = len(comlist)
    plist = splitall(parent)
    logging.debug(f'plist ={plist}')
    fulllist = splitall(fullpath)
    logging.debug(f'fulllist ={fulllist}')
    sublist = fulllist[comlength:]
    logging.debug(f'sublist ={sublist}')
    subpath = os.path.join(*sublist)
    logging.debug(f'subpath is {subpath}')
    return subpath


def splitall(path):
    '''
     Takes standard unix path and splits into list of components 
    '''
    if os.path.isdir(path):
        head = path
        plist = []
        while head != os.path.sep:
            head, tail = os.path.split(head)
            plist.append(tail)
        plist.reverse() 
        logging.debug(f'plist={plist}')
        return plist
    else:
        logging.error('this method is for directories only')
        return ['']

def flatten_tree(indir, outdir, mapfile):
    '''
     takes directory hierachy and puts all files into one outdir
     creates a mapfile of filenames (without extension) so that other output 
     can be restored to original hierarchy
    '''
    indir = os.path.abspath(indir)
    outdir = os.path.abspath(outdir)
    logging.debug(f'flattening {indir} to {outdir} with mapfile {mapfile}')
    
    logging.debug(f'confirming outdir {outdir} exists...')
    os.makedirs(outdir, exist_ok = True)    
    
    mapstring = ""
     
    for root, dirs, files in os.walk(indir):
        for d in dirs:
            ddir = os.path.join(root, d)    
            #logging.debug(f'{ddir}')
        for f in files:
            ffile = os.path.join(root, f) 
            logging.debug(f'{ffile}')
            dpath = os.path.dirname(ffile)
            subdir = remove_parentpath(indir, dpath)
            fname = os.path.basename(ffile)
            sample, ext = os.path.splitext(fname)
            mapstring += f'{subdir} {sample}\n'
            outfile = os.path.join(outdir, fname)
            logging.debug(f'copying {ffile} -> {outfile}')
            shutil.copyfile(ffile, outfile)
    #logging.debug(f'{mapstring}')   

    with open(mapfile, 'w') as mf:
        mf.write(mapstring)
     
def findmatches(dirpath, prefix, ext='*'):
    '''
    Return list of all absolute filepaths in dirpath that match prefix (minus file extension)
    '''
    logging.debug(f'looking for matches to {prefix} in {dirpath}')
    rlist = glob.glob(f'{dirpath}/{prefix}.{ext}')
    return rlist


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
    return df


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

                   
def unflatten_tree(indir, rootdir, mapfile, ext):
    '''
    takes files (ignoring extension) in indir and restores (copies) them to 
    hierarchy as stored in mapfile, based at rootdir. 

    '''
    indir = os.path.abspath(indir)
    rootdir = os.path.abspath(rootdir)
    logging.debug(f'unflattening {indir} to rootdir {rootdir} from mapfile {mapfile}')
    
    with open(mapfile, 'r') as mf:
        lines = mf.readlines()
        for line in lines:
            (subdir, sample ) = line.split()
            logging.debug(f'subdir = {subdir} sample={sample} ')
            flist = findmatches(indir, sample, ext)
            logging.debug(f'found matches={flist}')
            for fname in flist:
                bname = os.path.basename(fname)
                outfile = f'{rootdir}/{subdir}/{bname}'
                outdir = f'{rootdir}/{subdir}'
                os.makedirs(outdir, exist_ok = True)
                logging.debug(f'copying {fname} -> {outfile}')
                shutil.copyfile(fname, outfile)
            
            


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



def add_rowlist_column(rowlist, colval):
    """
    For use during dataframe construction. Adds col to list of rows with specified.
       
    """
    for row in rowlist:
        row.append(colval)
    return rowlist
    

def chmod_recurse(path, perms=0o755):
    """
    Recursively set permissions...
    0o755 is world read+execute. 
    
    """
    for root, dirs, files in os.walk(path):
        for d in dirs:
            os.chmod(os.path.join(root, d), perms)
        for f in files:
            os.chmod(os.path.join(root, f), perms)


def download_wget(srcurl, destpath, finalname=None, overwrite=True, decompress=True, rate='1M'):
    """
    GNU Wget 1.20.1, a non-interactive network retriever.
    Usage: wget [OPTION]... [URL]...
    
    Startup:
      -V,  --version                   display the version of Wget and exit
      -h,  --help                      print this help
      -v,  --verbose                   be verbose (this is the default)
      -nv, --no-verbose                turn off verboseness, without being quiet
           --report-speed=TYPE         output bandwidth as TYPE.  TYPE can be bits
      -t,  --tries=NUMBER              set number of retries to NUMBER (0 unlimits)
           --retry-connrefused         retry even if connection is refused
           --retry-on-http-error=ERRORS    comma-separated list of HTTP errors to retry
      -O,  --output-document=FILE      write documents to FILE
      -nc, --no-clobber                skip downloads that would download to
                                         existing files (overwriting them)
    
      -c,  --continue                  resume getting a partially-downloaded file
           --progress=TYPE             select progress gauge type
           --show-progress             display the progress bar in any verbosity mode
      -N,  --timestamping              don't re-retrieve files unless newer than
                                         local
           --no-if-modified-since      don't use conditional if-modified-since get
                                         requests in timestamping mode
           --no-use-server-timestamps  don't set the local file's timestamp by
                                         the one on the server
       -T,  --timeout=SECONDS           set all timeout values to SECONDS
           --dns-timeout=SECS          set the DNS lookup timeout to SECS
           --connect-timeout=SECS      set the connect timeout to SECS
           --read-timeout=SECS         set the read timeout to SECS
      -w,  --wait=SECONDS              wait SECONDS between retrievals
           --waitretry=SECONDS         wait 1..SECONDS between retries of a retrieval
           --random-wait               wait from 0.5*WAIT...1.5*WAIT secs between retrievals
    
           --limit-rate=RATE           limit download rate e.g. 1M  1 MB/s      
    """
    logging.debug(f'wget file {srcurl}')
    cmd = ['wget',
           '--no-verbose',
           '--no-use-server-timestamps',
           '--limit-rate', rate,
           '--continue', 
           '-O', f'{destpath}',
           f'{srcurl}']
    cmdstr = " ".join(cmd)
    logging.debug(f"wget command: {cmdstr} running...")
    
    start = dt.datetime.now()
    cp = subprocess.run(cmd, 
                        universal_newlines=True, 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
    end = dt.datetime.now()
    elapsed =  end - start
    logging.debug(f"ran cmd='{cmdstr}' return={cp.returncode} {type(cp.returncode)} ")
    if str(cp.returncode) == '0':
        logging.debug(f"got stderr: {cp.stderr}")
        logging.debug(f"got stdout: {cp.stdout}")
        if len(cp.stderr) > 10:
            dlbytes = parse_wget_output_bytes(cp.stderr)
            logging.info(f'downloaded {dlbytes} bytes {destpath} successfully, in {elapsed.seconds} seconds. ')
        else:
            logging.warn(f'file already downloaded.')
    else:
        logging.error(f'non-zero return code for src {srcurl}')
    return cp.returncode


def parse_wget_output_bytes(outstr):
    """
    E.g. 2021-07-20 14:33:09 URL:https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5529542/SRR5529542 [17019750/17019750] -> "SRR5529542.sra" [1]
    """
    logging.debug(f'handling stderr string {outstr}')    
    fields = outstr.split()
    bstr = fields[3][1:-1]
    dlbytes = int(bstr.split('/')[0])
    return dlbytes


def download_ftpurl(srcurl, destpath, finalname=None, overwrite=True, decompress=True):
    """
    destpath is directory
    
    Downloads via FTP from ftp src url to local destpath, 
    If finalname is specified, renames final output. 
    overwrite: won't re-download if filename already exists. 
    decompress: if filename ends with .gz , will gunzip  
    """
    log = logging.getLogger('star')
    # source FTP info
    (scheme, host, fullpath, p, q, f) = urllib.parse.urlparse(srcurl)
    filename = os.path.basename(fullpath)
    dirname = os.path.dirname(fullpath)
    
    # local files?
    localfile = f'{destpath}/{filename}'
    localfinal = f'{destpath}/{finalname}'
    destexists = os.path.exists(localfile) or os.path.exists(localfinal)
    log.debug(f'checking if {localfile} or {localfinal} exist -> {destexists}')
    
    if destexists and not overwrite:
        log.info(f"Destination files already exist and overwrite=false. Skipping.")
    else:
        log.info(
            f"Downloading file {filename} at path {dirname}/ on host {host} via FTP.")
        ftp = FTP(host)
        ftp.login('anonymous', 'hover@cshl.edu')
        ftp.cwd(dirname)
        log.debug(f'opening file {destpath}/{filename}. transferring...')
        with open(f'{destpath}/{filename}', 'wb') as fp:
            ftp.retrbinary(f'RETR {filename}', fp.write)
        log.debug(f"done retrieving {destpath}/{filename}")
        ftp.quit()
    
        if decompress and filename.endswith('.gz'):
            log.debug(f'decompressing gzip file {destpath}/{filename}')
            gzip_decompress(f'{destpath}/{filename}')
            os.remove(f'{destpath}/{filename}')
            filename = filename[:-3]
    
        if finalname is not None:
            src = "/".join([destpath, filename])
            dest = "/".join([destpath, finalname])
            log.info(f'renaming {src} -> {dest}')
            os.rename(src, dest)


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

 
def peek_tarball(tfile, subfile, numlines):
    """
    reads start of subfile in tarball without extracting.
    """
    
    cmd = ['tar', 
           '-xOf', 
           tfile,
           subfile, 
           '|',
           'zcat',
           '-d',
           '|',
           'head',
           f'-{numlines}'
           ]
    cmdstr = " ".join(cmd)
    #logging.debug(f"command: {cmdstr} running...")
    
    try:
        err, out, returncode = run_command_shell(cmd)
        out = out.decode()
        #logging.debug(f"got output:\n{out}")
        return err, out, returncode
        
    except Exception as e:
        logging.error(f'problem with {tfile}')
        logging.error(traceback.format_exc(None))


def remove_pathlist(pathlist):
    """
    recursively removes everything in pathlist
    if file, removes, 
    if directory, removes recursively. 
    
    """
    for fp in pathlist:
        if os.path.exists(fp):
            try:
                if os.path.isfile(fp):
                    os.remove(fp)
                elif os.path.isdir(fp):
                    shutil.rmtree(fp)
                logging.debug(f'removed {fp}')
            except Exception as ex:
                logging.error(f'problem removing {fp}')
                logging.error(traceback.format_exc(None))

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
    

def get_default_config():
    cp = ConfigParser()
    cp.read(os.path.expanduser("~/git/cshlwork/etc/utils.conf"))
    return cp


def get_configstr(cp):
    with io.StringIO() as ss:
        cp.write(ss)
        ss.seek(0)  # rewind
        return ss.read()


def write_fasta_from_df(config, df, outfile=None):
    '''
    Assumes df has 'sequence' column
    
    
    '''
    logging.debug(f'creating bowtie input')
    srlist = dataframe_to_seqlist(df)
    logging.debug(f'len srlist={len(srlist)}')
    if outfile is not None:
        SeqIO.write(srlist, outfile, 'fasta')
    else:
        logging.error('outfile is None, not implemented.')
    return outfile


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
            


def run_egad(go, nw, **kwargs):
    """EGAD running function
    
    Wrapper to lower level functions for EGAD

    EGAD measures modularity of gene lists in co-expression networks. 

    This was translated from the MATLAB version, which does tiled Cross Validation
    
    The useful kwargs are:
    int - nFold : Number of CV folds to do, default is 3, 
    int - {min,max}_count : limits for number of terms in each gene list, these are exclusive values


    Arguments:
        go {pd.DataFrame} -- dataframe of genes x terms of values [0,1], where 1 is included in gene lists
        nw {pd.DataFrame} -- dataframe of co-expression network, genes x genes
        **kwargs 
    
    Returns:
        pd.DataFrame -- dataframe of terms x metrics where the metrics are 
        ['AUC', 'AVG_NODE_DEGREE', 'DEGREE_NULL_AUC', 'P_Value']
    """
    assert nw.shape[0] == nw.shape[1] , 'Network is not square'
    assert np.all(nw.index == nw.columns) , 'Network index and columns are not in the same order'
    nw_mask = nw.isna().sum(axis=1) != nw.shape[1]
    nw = nw.loc[nw_mask, nw_mask].astype(float)
    np.fill_diagonal(nw.values, 1)
    return _runNV(go, nw, **kwargs)


def _runNV(go, nw, nFold=3, min_count=20, max_count=1000):

    #Make sure genes are same in go and nw
    genes_intersect = go.index.intersection(nw.index)

    go = go.loc[genes_intersect, :]
    nw = nw.loc[genes_intersect, genes_intersect]

    #Make sure there aren't duplicates
    duplicates = nw.index.duplicated(keep='first')
    nw = nw.loc[~duplicates, ~duplicates]

    go = go.loc[:, (go.sum(axis=0) > min_count) & (go.sum(axis=0) < max_count)]
    go = go.loc[~go.index.duplicated(keep='first'), :]

    roc = _new_egad(go.values, nw.values, nFold)

    col_names = ['AUC', 'AVG_NODE_DEGREE', 'DEGREE_NULL_AUC', 'P_Value']
    #Put output in dataframe
    return pd.DataFrame(dict(zip(col_names, roc)), index=go.columns)


def _new_egad(go, nw, nFold):

    #Build Cross validated Positive
    x, y = np.where(go)
    cvgo = {}
    for i in np.arange(nFold):
        a = x[i::nFold]
        b = y[i::nFold]
        dat = np.ones_like(a)
        mask = sparse.coo_matrix((dat, (a, b)), shape=go.shape)
        cvgo[i] = go - mask.toarray()
        
    CVgo = np.concatenate(list(cvgo.values()), axis=1)

    sumin = np.matmul(nw.T, CVgo)

    degree = np.sum(nw, axis=0)

    predicts = sumin / degree[:, None]

    np.place(predicts, CVgo > 0, np.nan)

    #Calculate ranks of positives
    rank_abs = lambda x: stats.rankdata(np.abs(x))
    predicts2 = np.apply_along_axis(rank_abs, 0, predicts)

    #Masking Nans that were ranked (how tiedrank works in matlab)
    predicts2[np.isnan(predicts)] = np.nan

    filtering = np.tile(go, nFold)

    #negatives :filtering == 0
    #Sets Ranks of negatives to 0
    np.place(predicts2, filtering == 0, 0)

    #Sum of ranks for each prediction
    p = np.nansum(predicts2, axis=0)

    #Number of predictions
    #Number of 1's masked for each GO term for each CV
    n_p = np.sum(filtering, axis=0) - np.sum(CVgo, axis=0)

    #Number of negatives
    #Number of GO terms - number of postiive
    n_n = filtering.shape[0] - np.sum(filtering, axis=0)
    roc = (p / n_p - (n_p + 1) / 2) / n_n
    U = roc * n_p * n_n
    Z = (np.abs(U - (n_p * n_n / 2))) / np.sqrt(n_p * n_n *
                                                (n_p + n_n + 1) / 12)
    roc = roc.reshape(nFold, go.shape[1])
    Z = Z.reshape(nFold, go.shape[1])
    #Stouffer Z method
    Z = np.nansum(Z, axis=0) / np.sqrt(nFold)
    #Calc ROC of Neighbor Voting
    roc = np.nanmean(roc, axis=0)
    P = stats.norm.sf(Z)

    #Average degree for nodes in each go term
    avg_degree = degree.dot(go) / np.sum(go, axis=0)

    #Calc null auc for degree
    ranks = np.tile(stats.rankdata(degree), (go.shape[1], 1)).T

    np.place(ranks, go == 0, 0)

    n_p = np.nansum(go, axis=0)
    nn = go.shape[0] - n_p
    p = np.nansum(ranks, axis=0)

    roc_null = (p / n_p - ((n_p + 1) / 2)) / nn

    return roc, avg_degree, roc_null, P


def MetaMarkers_PR(enrichment, class_pred = None):
    '''
    enrichment should be a dataframe of cells by cell type - from MetaMarkers
    '''
    # copy - otherwise overwrites the adata object if one is passed in....
    enr = enrichment.copy() 

    if class_pred is not None:
        groups = class_pred.predicted.unique()
        if len(groups) == 1 :
            cols = ~ enr.columns.str.contains(groups[0])
            enr.loc[:,cols] = 0
        else:
            for group, df in class_pred.groupby('predicted') :
                cols = ~ enr.columns.str.contains(group)
                enr.loc[df.index,cols]  = 0

    enr = enr.astype(float)
    pr = enr.T.melt().sort_values('value',ascending=False)
    pr = pr.reset_index(drop=True)
    pr['first_occurence'] = ~ pr.variable.duplicated()
    pr['TP'] = np.cumsum(pr['first_occurence'])
    pr['P'] = pr.index +1
    pr['Recall'] = pr.TP / enr.shape[0]
    pr['Precision'] = pr.TP / pr.P
    # print(np.trapz(pr.Precision,pr.Recall))
   
    return(pr)


#
#  SCQC/Bionformatic-specific functions. 
#  Assumes knowledge of dataframe formats. 
#


def taxon_to_spec(taxid= '10090'):
    d = {   '10090': "mouse",
            '9606':"human"}
    return(d[taxid])


def compare_barcode_to_whitelists(barcodes, 
        whitelistpaths = ['resource/whitelist_10xv1.txt','resource/whitelist_10xv2.txt','resource/whitelist_10xv3.txt']):
    f = open(whitelistpaths[0],'r').split('/n')
    # set(barcodes) & 
    pass



def read_identifier_file(filepath, flatten=True):
    """
    read and parse several formats for protein file
        
    1 per line
    multi per line: 
        if flatten=True just add to overall list. 
        otherwise each item in return list represents one line. multi-items in sub-list
    comma-separated
    ignore empty line
    ignore comment lines/extensions
    remove duplicates
    
    return list of items. 
    
    ABCD   EFGH
    IJK
    LMN
    
    -> [ ['ABCD', 'EFGH'],
    
    
    
    
    
    """
    idlist = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                idx = line.find('#')
                # get rid of comments
                if idx > -1:
                    line = line[:idx].strip()  
                if len(line) > 0:
                    if ',' in line:
                        line = line.replace(',',' ')
                    fields = line.split()
                    if flatten:
                        for f in fields:
                            f = f.strip()
                            if len(f) >0:
                                idlist.append(f)
                    else:
                        idlist.append(fields)
        idlist = list(set(idlist))
        logging.debug(f'got list with {len(idlist)} items.')
        return idlist
    except:
        return []    


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
    filepath = os.path.expanduser(filepath)
    df = pd.read_csv(filepath, sep='\t', index_col=0, keep_default_na=False, dtype =str, comment="#")
    df.fillna(value='', inplace=True)
    df = df.astype('str', copy=False)
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
        df = pd.read_csv(f, sep='\t', index_col=0)
        if newdf is None:
            newdf = df
        else:
            newdf = pd.concat([df, newdf], ignore_index=True, copy=False)
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


def write_config(config, filename, timestamp=True):
    '''
    writes config file to relevant name,
    if adddate, puts date/time code dot-separated before extension. e.g.
    filename = /path/to/some.file.string.txt  ->  /path/to/some.file.string.202303081433.txt
    date is YEAR/MO/DY/HR/MI
    
    '''
    filepath = os.path.abspath(filename)    
    dirname = os.path.dirname(filepath)
    basefilename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(basefilename) 
    
    if timestamp:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
        filename = f'{dirname}/{base}.{datestr}{ext}'

    os.makedirs(dirname, exist_ok=True)
        
    with open(filename, 'w') as configfile:
        config.write(configfile)
    logging.debug(f'wrote current config to {filename}')
        
    


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
    


def matrix2table(df, symmetric=True, merge_label=True, combo_char='x'):
    '''
    Takes  A x B  square matrix of permutations and converts to table of AxB with values.
    Symmetric means A-> B is same as B->A and need appear only once. 
    output is alphabetized on 'comb' column. 
    merge_labels means produce 'comb' column with combo_chr between label names.  

        A   B   C
    A   .8  .2  .1
    B  .2   .6  .3
    C  .1   .3  .5

    first  second  comb   val    
    A      A       AxA    .8
    A      B       AxB    .2
    A      C       AxC    .1
    B      B       BxB    .6
    B      C       BxC    .3
    C      C       CxC    .5
    
    '''
    #tdf = df.stack().reset_index()
    #tdf.columns = ['first','second','val']
    #tdf['combo'] = tdf['first'].astype(str) + combo_char + tdf['second'].astype(str)
    #items = list(set(list(tdf['first'])))
    items = list(df.columns)
    items.sort()

    lol = []
    
    for (f,s) in itertools.combinations_with_replacement(items, 2):
        r = [ f'{f}{combo_char}{s}', df.at[f,s] ]
        logging.debug(f'row is {r}')
        lol.append( r )

    tdf = pd.DataFrame(lol, columns=['pair','val'])
    return tdf


def listdiff(list1, list2):
    logging.debug(f"got list1: {list1} list2: {list2}")
    s1 = set(list1)
    s2 = set(list2)
    sd = s1 - s2
    dl = list(sd)
    dl.sort()
    logging.debug(f"diff has length {len(dl)}")
    return dl


def listmerge(list1, list2):
    logging.debug(f"got list1: {list1} list2: {list2}")
    s1 = set(list1)
    s2 = set(list2)
    sd = s1 | s2
    dl = list(sd)
    dl.sort()
    logging.debug(f"merged has length {len(dl)}")
    return dl

