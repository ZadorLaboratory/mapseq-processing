#!/usr/bin/env python
#
import argparse
import logging
import os
import sys

from configparser import ConfigParser

import pypandoc

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.utils import *
from mapseq.stats import *

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

    parser.add_argument('-L','--logfile', 
                    metavar='logfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Logfile for subprocess.')


    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Generated report file.') 

    parser.add_argument('-D','--datestr', 
                    metavar='datestr',
                    required=False,
                    default=None, 
                    type=str, 
                    help='Include datestr in relevant files.')

    parser.add_argument('-t','--template', 
                    metavar='template',
                    required=False,
                    default=os.path.expanduser('~/git/mapseq-processing/etc/user_report_template.md'), 
                    type=str, 
                    help='Report template.') 

    parser.add_argument('-m','--metadata', 
                    metavar='metadata',
                    required=False,
                    default=os.path.expanduser('~/git/mapseq-processing/etc/report_metadata.yaml'), 
                    type=str, 
                    help='Report template.')
   
    parser.add_argument('infiles',
                        metavar='infiles',
                        type=str,
                        nargs ="+",
                        help='All stats files to be consumed.')
        
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
    logging.debug(f'Running with config. {args.config}:\n{cdict}')
    logging.debug(f'infiles={args.infiles}')
    
      
    # set outdir / outfile
    if args.outfile is not None:
        logging.debug(f'outfile specified.')
        outfile = os.path.abspath(args.outfile)
        filepath = os.path.abspath(outfile)    
        dirname = os.path.dirname(filepath)
        filename = os.path.basename(filepath)
        (base, ext) = os.path.splitext(filename)   
        head = base.split('.')[0]
        outdir = dirname
        logging.debug(f'making missing outdir: {outdir} ')
        os.makedirs(outdir, exist_ok=True)
        logging.info(f'outdir={outdir} outfile={outfile}')
    else:
        outfile = sys.stdout
         
    if args.logfile is not None:
        log = logging.getLogger()
        FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
        formatter = logging.Formatter(FORMAT)
        logStream = logging.FileHandler(filename=args.logfile)
        logStream.setFormatter(formatter)
        log.addHandler(logStream)

    if args.datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    else:
        datestr = args.datestr
    
    logging.info(f'loading {args.infiles}') 
 
    rmd = generate_report(args.template, args.infiles, cp=cp)
    
    if outfile is None:
        fh = sys.stdout
        fh.write(rmd)
    else:
        logging.info(f'Saving report as Markdown doc to {outfile}...')
        with open(args.outfile, 'w') as fh:
            fh.write(rmd)
            
    # pandoc --metadata-file ~/git/mapseq-processing/etc/user_report_metadata.yaml  
    # statsreport.md -o statsreport.pdf
    if args.outfile is not None:
        logging.info(f'converting Markdown to PDF...')
        outfile = os.path.abspath(args.outfile)
        filepath = os.path.abspath(outfile)    
        dirname = os.path.dirname(filepath)
        filename = os.path.basename(filepath)
        (base, ext) = os.path.splitext(filename)
        pdfoutfile = os.path.join(dirname, f'{base}.pdf')           
        #rv = pypandoc.convert_file(args.outfile, 
        #                           'pdf', 
        #                           outputfile=pdfoutfile,
        #                           extra_args=[f'--metadata-file', f'{args.metadata}']
        #                           )
        #logging.info(f'writing PDF to {pdfoutfile} ')
        
        docxoutfile = os.path.join(dirname, f'{base}.docx') 
        rv = pypandoc.convert_file( args.outfile, 
                                    'docx',
                                    outputfile=docxoutfile,
                                    extra_args=[f'--metadata-file', f'{args.metadata}']
                                    )
    
