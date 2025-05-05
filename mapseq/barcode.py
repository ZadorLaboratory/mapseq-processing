#!/usr/bin/env python
#
#  Module for generic barcode/SSI handling, not specific to MAPseq or BARseq. 
#
#
import atexit
import datetime as dt
import logging
import os
import pprint 
import traceback

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from collections import defaultdict


from mapseq.stats import *
from mapseq.utils import *


def write_config(cp, filename, timestamp=True, datestring=None):
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
        cp.write(configfile)
    logging.debug(f'wrote current config to {filename}')
    return os.path.abspath(filename)


def split_fasta(cp, infile, barcodefile, outdir, force=False, datestr=None):
    '''
    Take paired FASTA file as input rather than raw FASTQ
    Use combined barcode sorter structure to minimize character comparisons. 

    '''
    if outdir is None:
        outdir = "."
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logging.debug(f'made outdir={outdir}')
    
    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    
    bclist = load_barcodes( barcodefile, 
                            labels=None, 
                            outdir=outdir, 
                            eol=True, 
                            max_mismatch=0,
                            cp=config)
    logging.info(f'made list of barcode handlers, length={len(bclist)}')
    logging.debug(bclist)
    
    output_exists = check_output(bclist)
    cfilename = f'{outdir}/splitter.config.txt'
    bc_length=len(bclist[0].barcode)
    unmatched = os.path.abspath(f'{outdir}/unmatched.fasta')
    
    seqhandled_interval = int(cp.get('splitter','seqhandled_interval')) 
    matched_interval = int(cp.get('splitter','matched_interval'))
    unmatched_interval = int(cp.get('splitter','unmatched_interval'))
    
    logging.info(f'performing split. outdir={outdir} output_exists={output_exists} force={force} ')
    
    if ( not output_exists ) or force:
        write_config(cp, cfilename, timestamp=True, datestring=datestr)
        sh = StatsHandler(cp, outdir=outdir, datestr=datestr)        
        # create structure to search for all barcode simultaneously
        labeldict = {}
        for bco in bclist:
            labeldict[bco.barcode] = bco.label
        matchdict, seqdict, unmatched_file = build_bcmatcher(bclist) 

        pairshandled = 0
        num_handled_total = 0
        num_matched_total = 0
        num_unmatched_total = 0

        # handle all sequences in input 
        with open(infile) as f:
            num_handled = 0
            num_matched = 0
            num_unmatched = 0
            
            logging.debug(f'handling file {infile}')
            while True:
                try:
                    line = f.readline()
                    if line.startswith('>'):
                        pass
                    else:
                        if len(line) == 0:
                            raise StopIteration
                        
                        fullread = line.strip()
                        # handle sequence
                        matched = False
                        seq = fullread[-bc_length:]
                        # id, seq, matchdict, fullseq, unmatched
                        matched = do_match(num_handled, seq, matchdict, fullread, unmatched_file)
                        num_handled += 1
                        num_handled_total += 1
    
                        # report progress...                    
                        if not matched:
                            num_unmatched += 1
                            num_unmatched_total += 1                        
                            if num_unmatched % unmatched_interval == 0:
                                logging.debug(f'{num_unmatched} unmatched so far.')
                        else:
                            num_matched += 1
                            num_matched_total += 1
                            if num_matched % matched_interval == 0:
                                logging.debug(f'match {num_matched}: found SSI in {fullread}!')                        
    
                        if num_handled % seqhandled_interval == 0: 
                            logging.info(f'handled {num_handled} matched={num_matched} unmatched={num_unmatched}')
                    
                except StopIteration as e:
                    logging.debug('iteration stopped')    
                    break
            
        logging.debug(f'finished with {infile}')
        sh.add_value('/fasta','reads_handled', num_handled_total )
        sh.add_value('/fasta','reads_unmatched', num_unmatched_total )
        f.close()
        
        matchrate = 0.0
        if num_matched_total > 0: 
            unmatchrate = num_unmatched_total / num_matched_total              
            matchrate = 1.0 - unmatchrate
        logging.info(f'handled {num_handled_total} sequences. {num_matched_total} matched. {num_unmatched_total} unmatched matchrate={matchrate}')
    else:
        logging.warn('All FASTA output exists and force=False. Not recalculating.')


def set_ssi_df(df, column='ssi'):
    '''
    assumed ssi 
    sets barcode id and label from sequence. 
 
    '''
    initlen = len(seqdf)
    logging.debug(f'setting read counts for sequence DF len={len(seqdf)}')
    cdf = make_read_counts_df(cp, seqdf, label=None)
    fmap = cdf.set_index('sequence')['read_count'].to_dict()
    seqdf['read_count'] = seqdf['sequence'].map(fmap)
    # change made to inbound df, but return anyway
    return seqdf
    

def build_bcmatcher(bclist):
    '''
    Take list of BarcodeHandler objects, and make matching dict of dicts:
    
    AGAT
    ACGA
    CGGT
                      a         c 
                   c    g       g
                    g  a        g
                   a    t       t 
    
    
    seqdict    { 'AGAT' : [],
                  'ACGA' : [],
                  'CGGT': []
                 }
                 
    unmatched  list      
    '''
    matchdict = {}
    seqdict = {}

    #unmatched = []    
    # calculate directory
    afile = bclist[0].filename
    filepath = os.path.abspath(afile)    
    dirname = os.path.dirname(filepath)
    unmatchedfile = f'{dirname}/unmatched.fasta'
    unmatched = open(unmatchedfile,'w')
    logging.debug(f'opened {unmatchedfile} for writing...')


    for bco in bclist:
        seq = bco.barcode
        seqlen = len(seq)
        stoplen = seqlen - 1

        subdict = matchdict
        for i, nt in enumerate(seq):
            logging.debug(f'i={i} nt={nt}')
            if i < stoplen:
                try:
                    subdict = subdict[nt]
                    logging.debug(f'subdict = {subdict}')
                except KeyError:
                    subdict[nt] = {}
                    subdict = subdict[nt]
            else:
                logging.debug(f'end condition reached')
                #subdict[nt] = [ ]
                subdict[nt] = bco
                seqdict[seq] = subdict[nt]
            
        #print(matchdict)
    return (matchdict, seqdict, unmatched)


def do_match(id, seq, matchdict, fullseq, unmatched):
    '''
    given a target sequence, walk the matchdict and either place a (id, fullseq) tuple in the list 
    at the end of the seq or add the (id, fullseq) tuple to the unmatched list. 
    
    return whether matched or not True | False
    
    '''
    subdict = matchdict
    seqlen = len(seq)
    stoplen = len(seq) - 1
    
    for i, nt in enumerate(seq):
        #logging.debug(f'i={i} nt={nt}')
        if i < stoplen:
            try:
                subdict = subdict[nt]
            except KeyError:
                #unmatched.append((id, fullseq))
                unmatched.write(f'>{id}\n{fullseq}\n')
                return False
        else:
            try:
                #seqlist = subdict[nt]
                # seqlist.append((id, fullseq[:-seqlen]))
                bco = subdict[nt]
                bco.writeseq(id, fullseq[:-seqlen] )
                #bco.of.write(f'>{id}\n{fullseq[:-seqlen]}\n')
                return True
            
            except KeyError:
                #unmatched.append((id, fullseq))
                # unmatched is a file object. 
                unmatched.write(f'>{id}\n{fullseq}\n')
                return False
        
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
        self.max_mismatch = int(max_mismatch)
        self.of = None  # do not open file until necessary. 
        self.dataframe = None
        if outdir is None:
            outdir = "."
        atexit.register(self.cleanup)

    def cleanup(self):
        logging.debug(f'running cleanup(). Closing outfile {self.filename} ')
        if self.of is not None:
            self.of.close()

    def do_match(self, id, seq ):
        '''
        test sequence exactly against this barcode.
        just test end of sequence the same length as barcode. 
        only write part of sequence up to barcode/SSI.
        EOL testing only now, add BOL...
         
        '''
        if self.of is None:
            logging.debug(f'opening new file for {self.label}')
            self.of = open(self.filename, 'w')
            self.of.write(f'sequences for  barcode {self.label}\n')
        
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
    
    
    def writeseq(self, header, sequence ):
        if self.of is None:
            logging.debug(f'opening new file for {self.label}')
            self.of = open(self.filename, 'w')
        self.of.write(f'>{header}\n{sequence}\n')        
        
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


def load_barcodes(bcfile, 
                  labels=None, 
                  outdir=None, 
                  eol=True, 
                  max_mismatch=0,
                  cp=None):
    '''
     labellist is a python list of ids. barcode labels will be checked against
     <id> from labellist. 
    '''
    if cp is None:
        cp = get_default_config()
    if labels is not None:
        codelist = [ f'{x}' for x in labels]
        logging.debug(f'codelist={codelist}')
    bcdict = get_barcode_dict(bcfile)
    bclist = []
    #codelist = list(bcdict.values())
    for bcseq in bcdict.keys():
        label = bcdict[bcseq]
        logging.debug(f'handling label={label} and barcode={bcseq}')        
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
        #if os.path.exists(bch.filename) and ( os.path.getsize(bch.filename) >= 1 ):
        nonzero = False
        pathexists = False
        pathexists = os.path.exists(bch.filename)
        if pathexists:
            nonzero = os.path.getsize(bch.filename) >= 1
        
        logging.info(f'file={bch.filename} pathexists={pathexists} nonzero={nonzero}')
        if pathexists and nonzero : 
            logging.info(f'Non-empty BC{bch.label}.fasta exists.')
        else:
            logging.info(f"{bch.filename} doesn't exist. output_exists=False")
            missing.append(bch.label)
            output_exists = False
    
    if not output_exists:
        logging.debug(f'missing BC labels: {missing}')
    else:
        logging.info('all output exists.')
    return output_exists


def get_rtprimer_dict(bcfile, labels=None):
    '''
    read barcode file, parse into dict. 
    optionally include only labels in provided list. 
    '''
    rtdict = {}  # map from seq to primer number
    with open(bcfile) as fh:
        logging.debug(f"opened barcode file {bcfile}")
        while True:
            ln = fh.readline()
            if len(ln) < 2:
                break
            (label, bcseq) = ln.split()
            rtn = label.replace('BC','')
            if labels is None or label in labels: 
                rtdict[bcseq] = rtn
                logging.debug(f'label={label} rtn={rtn} barcode={bcseq}')    
    logging.debug(f'got dict len={len(rtdict)}')
    return rtdict


def get_barcode_dict(bcfile, labels=None):
    '''
    read barcode file, parse into dict. 
    optionally include only labels in provided list. 
    '''
    bcdict = {}  # map from seq to barcode label
    with open(bcfile) as fh:
        logging.debug(f"opened barcode file {bcfile}")
        while True:
            ln = fh.readline()
            if len(ln) < 2:
                break
            (label, bcseq) = ln.split()
            rtn = label.replace('BC','')
            if labels is None or label in labels: 
                bcdict[bcseq] = label
                logging.debug(f'label={label} rtn={rtn} barcode={bcseq}')    
    logging.debug(f'got dict len={len(bcdict)}')
    return bcdict



