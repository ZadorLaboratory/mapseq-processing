import argparse
import logging
import os
import sys

from configparser import ConfigParser

from mapseq.core import *
from mapseq.barcode import *
from mapseq.utils import *
from mapseq.stats import *

# default is single infile
INFILES_COMMANDS = ['aggregate_reads_byfile',
                    'merge_mapseq_dfs',
                    'merge_stats',
                    'process_fastq_pairs',
                    'process_fastq_pairs_concat']

# default is single outfile
OUTDIR_COMMANDS = ['aggregate_reads_byfile',
                   'make_matrices',
                   'process_fastq_pairs']

SAMPLEINFO_COMMANDS = ['process_fastq_pairs',
                       'process_fastq_pairs_concat',
                       'make_readtable',
                       'make_vbctable']

#
#   Command line helper functions 
# 
def initialize_logging():
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARN)    

def set_logging(args):
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    if args.logfile is not None:
        log = logging.getLogger()
        FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
        formatter = logging.Formatter(FORMAT)
        logStream = logging.FileHandler(filename=args.logfile)
        logStream.setFormatter(formatter)
        log.addHandler(logStream)

def handle_sampleinfo(args, cp):
    sampdf = None
    if args.sampleinfo is not None:
        logging.debug(f'loading sample DF...')
        sampdf = load_sample_info(args.sampleinfo, cp=cp)
        logging.debug(f'\n{sampdf}')
        # sampdf.to_csv(f'{outdir}/sampleinfo.tsv', sep='\t')
    return sampdf
 
def build_argparser(cmd):
    logging.debug(f'Building args for {cmd}')

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

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')

    parser.add_argument('-L','--logfile', 
                    metavar='logfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Logfile for subprocess.')

    parser.add_argument('-f','--force', 
                    action="store_true", 
                    default=False, 
                    help='Recalculate even if output exists.')

    if cmd in SAMPLEINFO_COMMANDS:
        parser.add_argument('-s','--sampleinfo', 
                        metavar='sampleinfo',
                        required=False,
                        default=None,
                        type=str, 
                        help='XLS or TSV sampleinfo file. ')    

        parser.add_argument('-S','--samplesheet', 
                        metavar='samplesheet',
                        required=False,
                        default='Sample information',
                        type=str, 
                        help='XLS sheet tab name.')

    if cmd in INFILES_COMMANDS:
        parser.add_argument('infiles' ,
                        metavar='infiles', 
                        type=str,
                        nargs='+',
                        default=None, 
                        help='One or more input files.')
    else:
        parser.add_argument('infile',
                        metavar='infile',
                        type=str,
                        help='Single input file.')               

    if cmd in OUTDIR_COMMANDS:
        parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=os.path.abspath('./'), 
                    type=str, 
                    help='Output directory. Current directory if omitted.')        
    else:
        parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=True,
                    default=None, 
                    type=str, 
                    help='Output file name.')        

    return parser

def handle_config(args):
    cp = ConfigParser()
    o = cp.read(args.config)
    if len(o) < 1:
        logging.error(f'No valid configuration. {args.config}')
        sys.exit(1)
    cdict = format_config(cp)
    logging.debug(f'Running with config. {args.config}: {cdict}')
    return cp

def parse_input_output(args):
    '''
    Check args for infiles, outfile, outdir
    outfile and outdir are mutually exclusive. 
    @return infile, infiles, outfile, outdir
    
    '''
    infile = None
    infiles = None
    outfile = None
    outdir = None

    if 'outfile' in args:
        outfile = os.path.abspath( os.path.expanduser( args.outfile ) )
        filepath = os.path.abspath(outfile)    
        dirname = os.path.dirname(filepath)
        filename = os.path.basename(filepath)
        outdir = dirname
        logging.debug(f'outdir implicitly set to {outdir}')
    elif 'outdir' in args:
        outdir = os.path.abspath(  os.path.expanduser( args.outdir ) )
        logging.debug(f'outdir explicitly set to {outdir}')
    else:
        logging.warning('Neither outfile or outdir in args namespace.')

    if outdir is not None:
        logging.debug(f'making missing outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)

    if 'infile' in args:
        infile = args.infile
    elif 'infiles' in args:
        infiles = args.infiles

    return (infile, infiles, outfile, outdir)


#
# Pipeline command line functions.
#
def aggregate_reads_byfile():
    '''
    Aggregate reads within each input file. 
    '''
    cmd = 'aggregate_reads_byfile'
    initialize_logging()
    parser = build_argparser(cmd=cmd)

    # CUSTOM ARG(S)
    parser.add_argument('-k', '--use_dask', 
                        action="store_true", 
                        dest='use_dask',
                        help='Use Dask subsystem.')

    parser.add_argument('-t','--dask_temp', 
                    metavar='dask_temp',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Fast, roomy storage for DASK temp files.')

    parser.add_argument('-m','--min_reads', 
                        metavar='min_reads',
                        required=False,
                        default=None,
                        type=int, 
                        help='Min reads to retain initial full read.')

    parser.add_argument('-C','--column', 
                    metavar='column',
                    required=False,
                    default=['sequence','source'], 
                    type=str, 
                    help='Read column to aggregate.')

    args = parser.parse_args()
    set_logging(args)
    logging.debug('Handling command {cmd}')
    cp = handle_config(args)

    # COMMAND-SPECIFIC CUSTOMIZATIONS
    # override via config if needed. 

    # STANDARDIZED STUFF sampleinfo, infile(s), outfile(s)/outdir
    (infile, infiles, outfile, outdir) = parse_input_output(args)

    for infile in infiles:
        dirpath, filename = os.path.split(infile)
        base, ext = os.path.splitext(filename)
        outbase = base.replace('reads','aggregated')
        outfile = os.path.join(outdir, f'{outbase}.tsv')
        logging.info(f'calculated outfile = {outfile}')

        if (not os.path.exists(outfile) or args.force):
            df = load_mapseq_df( infile, fformat='reads', use_dask=args.use_dask, use_arrow=True)
            logging.debug(f'loaded. len={len(df)} dtypes =\n{df.dtypes}') 
            df = aggregate_reads( df, 
                                  column=args.column,
                                  outdir=outdir,
                                  min_reads = args.min_reads,
                                  use_dask = args.use_dask, 
                                  dask_temp = args.dask_temp,
                                  cp=cp 
                                )
            logging.debug(f'got aggregated df len={len(df)}:\n{df} ')
            write_mapseq_df(df, outfile)
        else:
            logging.info(f'Output {outfile} exists and force={args.force} Skipping.')
    logging.info('Done aggregate_reads')

def process_fastq_pairs():
    '''
    Read in FASTQ pairs and join them.
    Write to outdir.
    Do not concat. 
    
    '''
    cmd = 'process_fastq_pairs'
    initialize_logging()
    parser = build_argparser(cmd=cmd)

    # CUSTOM ARG(S)
    parser.add_argument('-F','--filter_non_dominant', 
                    action="store_true", 
                    default=None, 
                    help='Remove non-dominant values.')

    args = parser.parse_args()
    set_logging(args)
    logging.debug('Handling command {cmd}')
    cp = handle_config(args)

    # COMMAND-SPECIFIC CUSTOMIZATIONS
    # override via config if needed. 
    if args.filter_non_dominant is not None:
        cp.set('fastq','filter_non_dominant', str(args.filter_non_dominant))

    cp.set('fastq','write_each', str( True ))

    # check nargs. 
    if (len(args.infiles) < 2) or (len(args.infiles) % 2 != 0 ):
        parser.print_help()
        print('Error: the following arguments are required: 2 or multiple of 2 infiles')
        sys.exit(1)

    # STANDARDIZED STUFF sampleinfo, infile(s), outfile(s)/outdir
    sampdf = handle_sampleinfo(args, cp)
    (infile, infiles, outfile, outdir) = parse_input_output(args)

    process_fastq_grouped_noconcat( infiles, 
                                    outdir, 
                                    sampdf = sampdf,
                                    force = args.force,
                                    cp = cp, 
                                    )
    logging.info('Done')

def process_fastq_pairs_concat():
    '''
    Read in FASTQ pairs and join them.
    Write concatenated output to outfile.

    '''
    cmd = 'process_fastq_pairs_concat'
    initialize_logging()
    parser = build_argparser(cmd=cmd)

    parser.add_argument('-w','--write_each', 
                    action="store_true", 
                    default=None, 
                    help='Write output for each input pair. ')

    parser.add_argument('-F','--filter_non_dominant', 
                    action="store_true", 
                    default=None, 
                    help='Remove non-dominant values.')

    args = parser.parse_args()
    set_logging(args)
    logging.debug('Handling command {cmd}')
    cp = handle_config(args)

    # COMMAND-SPECIFIC CUSTOMIZATIONS
    # override via config if needed. 
    if args.filter_non_dominant is not None:
        cp.set('fastq','filter_non_dominant', str(args.filter_non_dominant))

    if args.write_each is not None:
        cp.set('fastq','write_each', str( args.write_each ))

    # check nargs. 
    if (len(args.infiles) < 2) or (len(args.infiles) % 2 != 0 ):
        parser.print_help()
        print('Error: the following arguments are required: 2 or multiple of 2 infiles')
        sys.exit(1)

    # STANDARDIZED STUFF sampleinfo, infile(s), outfile(s)/outdir
    sampdf = handle_sampleinfo(args, cp)
    (infile, infiles, outfile, outdir) = parse_input_output(args)

    process_fastq_grouped_concat(   args.infiles, 
                                    outfile,
                                    sampdf = sampdf,                         
                                    force = args.force, 
                                    cp=cp)
    logging.info('Done')

    

