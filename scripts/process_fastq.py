#!/usr/bin/env python


import argparse
import logging
import os

from configparser import ConfigParser

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# https://github.com/chapmanb/bcbb/tree/master/gff
#from BCBio import GFF
#from BCBio.GFF import *

class BarCodeHandler(object):
    
    def __init__(self, label, barcode):
        self.barcode = barcode
        self.label = label
        self.of = open(f'{label}.fasta', 'w')
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
            

def load_barcodes(config, bcfile):
    bclist = []
    with open(bcfile) as fh:
        while True:
            ln = fh.readline()
            if len(ln) < 2:
                break
            (label, bcseq) = ln.split()
            bch = BarCodeHandler(label, bcseq)
            bclist.append(bch)
    logging.debug(f'made list of {len(bclist)} barcode handlers.')
    return bclist
        

def process_fastq_pair(config, read1file, read2file, bclist):

    umf = open('unmatched.fasta', 'w')

    r1s = int(config.get('fastq','r1start'))
    r1e = int(config.get('fastq','r1end'))
    r2s = int(config.get('fastq','r2start'))
    r2e = int(config.get('fastq','r2end'))   

    recs1 =  SeqIO.parse(read1file, "fastq")
    recs2 = SeqIO.parse(read2file, "fastq")

    seqshandled = 0
    unmatched = 0
    
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
                    logging.debug(f'found bc {bch.label} in {fullread}!!')
                    matched = True
                    break
            if not matched:
                unmatched += 1
                logging.debug(f'found unmatched. writing...')
                id = str(seqshandled)
                sr = SeqRecord( fullread, id=id, name=id, description=id)
                SeqIO.write([sr], umf, 'fasta')
            
            seqshandled += 1
            if seqshandled % 50000 == 0: 
                logging.info(f'handled {seqshandled} read sequences. {unmatched} unmatched so far.')
        
        except StopIteration as e:
            logging.debug(f'iteration stopped')
            break
    
    umf.close()
    for bch in bclist:
        bch.finalize()    
    
    logging.info(f'handled {seqshandled} read sequences. {unmatched} unmatched')




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
                        default=os.path.expanduser('barcode_v2.txt'),
                        type=str, 
                        help='barcode file space separated.')

    parser.add_argument('infiles' ,
                        metavar='infiles', 
                        type=str,
                        nargs=2,
                        default=None, 
                        help='Read1file  Read2file')
       
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
    bclist = load_barcodes(cp, args.barcodes)
    process_fastq_pair(cp, args.infiles[0], args.infiles[1], bclist)
    
    
    
    