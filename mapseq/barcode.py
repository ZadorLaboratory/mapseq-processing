#!/usr/bin/env python
#
#  Module for generic barcode/SSI handling, not specific to MAPseq or BARseq. 
#
#

import logging
import os
import traceback

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


class BarcodeHandler(object):
    '''
    Basically implements fastx_barcode_handler.pl. 
    Matches against given barcode/SSI sequence, and
    writes out target fasta (minus SSI) to SSI-specific file.
    
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

    def __str__(self):
        s = f"BarCodeHandler: label={self.label} barcode={self.barcode} df=\n{self.dataframe}"
        return s

    def __repr__(self):
        s = f"BarCodeHandler: label={self.label} barcode={self.barcode} df=\n{self.dataframe}"
        return s    

    @classmethod
    def make_all_bch_counts(cls, config, ssilist, outdir=None):
        '''
        Makes counts dataframes for all barcode handler objects. 
        Saves 
        '''
        alldf = pd.DataFrame({'label': pd.Series(dtype='str'),
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
            cdf['label'] = bch.label       
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
        alldf = pd.DataFrame({'label': pd.Series(dtype='str'),
                   'sequence': pd.Series(dtype='str'),
                   'counts': pd.Series(dtype='int')})
        for df in dflist:
            alldf = pd.concat([alldf, df], ignore_index = True, copy=False)
        return alldf


def load_barcodes(config, bcfile, labels = None, outdir=None, eol=True, max_mismatch=0):
    '''
     labellist is a python list of ids. barcode labels will be checked against
     <id> from labellist. 
    '''
    codelist = [ f'{x}' for x in labels]
    logging.debug(f'codelist={codelist}')
    bclist = []
    with open(bcfile) as fh:
        logging.debug(f"opened barcode file {bcfile}")
        while True:
            ln = fh.readline()
            if len(ln) < 2:
                break
            (label, bcseq) = ln.split()
            logging.debug(f'handling label={label} and bcseq={bcseq}')
            
            if labels is None or label in codelist:
                bch = BarcodeHandler(label, bcseq, outdir, eol, max_mismatch)
                bclist.append(bch)
    logging.debug(f'made list of {len(bclist)} barcode handlers.')
    return bclist

def check_output(bclist):
    '''
    check to see if valid output exists for all input barcode handlers, in 
    order to facility short-circuiting. 
    '''
    output_exists = True
    missing = []
    for bch in bclist:
        logging.debug(f'checking path {bch.filename}...')
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