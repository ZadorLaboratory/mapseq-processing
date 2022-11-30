#!/usr/bin/env python
#
#


import argparse
import gzip
import logging
import os
import sys

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from configparser import ConfigParser

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from cshlwork.utils import JobRunner, JobStack, JobSet

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
            SeqIO.write([sr], self.of, 'fasta')
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
            if seqshandled % 50000 == 0: 
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




if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARN)
    
    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')
    
    parser.add_argument('-c','--config', 
                        metavar='config',
                        required=False,
                        default=os.path.expanduser('~/git/mapseq-processing/etc/mapseq.conf'),
                        type=str, 
                        help='out file.')    

    parser.add_argument('-b','--barcodes', 
                        metavar='barcodes',
                        required=False,
                        default=os.path.expanduser('~/git/mapseq-processing/etc/barcode_v2.txt'),
                        type=str, 
                        help='barcode file space separated.')

    parser.add_argument('-s','--sampleinfo', 
                        metavar='sampleinfo',
                        required=True,
                        default=None,
                        type=str, 
                        help='XLS sampleinfo file. ')

    parser.add_argument('-o','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. cwd if not given.') 


    parser.add_argument('infiles' ,
                        metavar='infiles', 
                        type=str,
                        nargs=2,
                        default=None, 
                        help='Read1 and Read2 fastq files')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    cp.read(args.config)
    cdict = {section: dict(cp[section]) for section in cp.sections()}
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'infiles={args.infiles}')
    
    sampdf = load_sample_info(cp,args.sampleinfo)
    logging.debug(sampdf)
    rtlist = list(sampdf.rtprimer.dropna())
        
    bcolist = load_barcodes(cp, args.barcodes, labels=rtlist, outdir=args.outdir)
    logging.debug(bcolist)
    process_fastq_pair(cp, args.infiles[0], args.infiles[1], bcolist, outdir=args.outdir)
    make_summaries(cp, bcolist)
    
    
    
    
    