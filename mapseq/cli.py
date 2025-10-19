import argparse
import logging
import os
import sys

from configparser import ConfigParser

import pandas as pd

from mapseq.core import *
from mapseq.barcode import *
from mapseq.utils import *
from mapseq.stats import *


def process_fastq_pairs():
    
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

    parser.add_argument('-s','--sampleinfo', 
                        metavar='sampleinfo',
                        required=False,
                        default=None,
                        type=str, 
                        help='XLS or TSV sampleinfo file. ')

    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Combined read, read_count TSV')   

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. output file base dir if not given.')

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')

    parser.add_argument('-f','--force', 
                    action="store_true", 
                    default=False, 
                    help='Recalculate even if output exists.') 

    parser.add_argument('-w','--write_each', 
                    action="store_true", 
                    default=None, 
                    help='Write output for each input pair. ')

    parser.add_argument('-F','--filter_non_dominant', 
                    action="store_true", 
                    default=None, 
                    help='Remove non-dominant values.')

    parser.add_argument('-L','--logfile', 
                    metavar='logfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Logfile for subprocess.')
   
    parser.add_argument('infiles' ,
                        metavar='infiles', 
                        type=str,
                        nargs='+',
                        default=None, 
                        help='Read1 and Read2 [Read1B  Read2B ... ] fastq files')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   

    cp = ConfigParser()
    o = cp.read(args.config)
    if len(o) < 1:
        logging.error(f'No valid configuration. {args.config}')
        sys.exit(1)
    cdict = format_config(cp)
    logging.debug(f'Running with config. {args.config}: {cdict}')
    logging.debug(f'force={args.force} infiles={args.infiles}')

    # override via config if needed. 
    if args.filter_non_dominant is not None:
        cp.set('fastq','filter_non_dominant', str(args.filter_non_dominant))

    if args.write_each is not None:
        cp.set('fastq','write_each', str( args.write_each ))

    # check nargs. 
    if (len(args.infiles) < 2)  or (len(args.infiles) % 2 != 0 ):
        parser.print_help()
        print('error: the following arguments are required: 2 or multiple of 2 infiles')
        sys.exit(1)
       
    # set outdir / outfile
    outdir = os.path.abspath('./')
    outfile = f'{outdir}/reads.tsv'
    if args.outdir is None:
        if args.outfile is not None:
            logging.debug(f'outdir not specified. outfile specified.')
            outfile = os.path.abspath(args.outfile)
            filepath = os.path.abspath(outfile)    
            dirname = os.path.dirname(filepath)
            filename = os.path.basename(filepath)
            (base, ext) = os.path.splitext(filename)   
            head = base.split('.')[0]
            outdir = dirname
            logging.debug(f'outdir set to {outdir}')
        else:
            logging.debug(f'outdir/file not specified.')        
    else:
        # outdir specified
        outdir = os.path.abspath(args.outdir)
        if args.outfile is not None:
            logging.debug(f'outdir specified. outfile specified.')
            outfile = os.path.abspath(args.outfile)
        else:
            logging.debug(f'outdir specified. outfile not specified.')
            outfile = f'{outdir}/reads.tsv'

    logging.debug(f'making missing outdir: {outdir} ')
    os.makedirs(outdir, exist_ok=True)
    logging.info(f'outdir={outdir} outfile={outfile}')

    logging.info(f'handling {args.infiles[0]} and {args.infiles[1]} to outdir {args.outfile}')
    infilelist = package_pairs(args.infiles)   
    
    logging.debug(f'infilelist = {infilelist}')
    
    if args.logfile is not None:
        log = logging.getLogger()
        FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
        formatter = logging.Formatter(FORMAT)
        logStream = logging.FileHandler(filename=args.logfile)
        logStream.setFormatter(formatter)
        log.addHandler(logStream)

    if args.sampleinfo is not None:
        logging.debug(f'loading sample DF...')
        sampdf = load_sample_info(args.sampleinfo, cp=cp)
        logging.debug(f'\n{sampdf}')
        sampdf.to_csv(f'{outdir}/sampleinfo.tsv', sep='\t')
    else:
        sampdf = None
        logging.warning('No sampleinfo. Cannot guarantee SSI per FASTQ pair.')
    
        
    if args.datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    else:
        datestr = args.datestr
    sh = StatsHandler(outdir=outdir, datestr=datestr) 

    df = process_fastq_grouped( infilelist, 
                                outdir,
                                sampdf = sampdf,                          
                                force=args.force,
                                cp = cp)
    
    write_mapseq_df(df, outfile)    
    logging.info('Done process_fastq_pairs')  