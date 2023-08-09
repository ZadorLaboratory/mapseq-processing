#!/usr/bin/env python
# 
#  Parses current standard NGS CSHL notification email, and
#  retrieves data on GPFS by running cp remotely over ssh. 
#  
#  Requires password-less login with SSH key active on Mac:
#     https://apple.stackexchange.com/questions/48502/how-can-i-permanently-add-my-ssh-private-key-to-keychain-so-it-is-automatically
#  Requires Python
#  Requires pandas
#  Recommend setting up  within a Conda environement.  
#
#  Produces path table:
# 
# cat paths.20221207.tsv | while read line; do echo $line | xargs bash -c 'mkdir $0/$1' ; done
# cat paths.20221207.tsv | while read line; do echo $line | xargs bash -c 'cp -vr $2/* $0/$1/' ; done
#
#
#
#
import argparse
import email
import logging
import os
import sys
import traceback

import paramiko
import pandas as pd

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)
gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

#
# If Yield (MBases) > 50000000 then nextseq, else miseq
#
NEXTSEQ_THRESHOLD = 50000000
PATHTABLE_COLUMNS = ['seqtype','expid','srcdir']
# OUTROOT = '/grid/mbseq/data_norepl/mapseq/library_data'


def email_to_table(file):
    '''
    take list of .eml files ->  .tsv file   (no index, no headers)
    <studyid>  <path>  <studytype>
    
    path grep is 'Zador_Lab', following field is studyid
    
    '''
    lol = []
    if file.endswith('.eml'):
        emailtext = eml_to_text(file)
        items = parse_ngs_emailtxt(emailtext)
        if items is not None:
            lol.append(items)
    else:
        logging.warning(f'file {file} not .eml file extension. ignoring.')
    logging.debug(f'made list of {len(lol)} row(s)...')
    df = pd.DataFrame(data=lol, index=None, columns=PATHTABLE_COLUMNS)
    return df

def parse_ngs_emailtxt(emailtext):
    '''
    *ASSUMES*
    -- Input in .eml format with that extension
    -- Body of email with info is last plain/text Content-Type part of multipart email 
    -- 'Here are' line with expid in parentheses
    -- 'Zador_Lab' in remote file path. 
    -- last directory in remote path is libid
    -- collapsed line following line with 'MBases' is number. 
    
    returns studyid, path, studytype
    
    studytype = nextseq  miseq
    '''
    #logging.debug(emailtext)
    lines = emailtext.split('\n')
    logging.debug(f'e.g.lines {lines[6:10]}')
    expid = None
    studyid = None
    path = None
    libid = None
    studytype = 'miseq'
    for idx, line in enumerate(lines):
        if line.startswith('Here are'):
            logging.debug('Found "Here are" line...')
            #print(line)
            # get experiment id from note. 
            a = line.find('(')
            b = line.find(')')
            expid = line[a+1:b]
            expid = expid.strip()
            logging.debug(f'found exp id {expid}')
                      
        elif 'Zador_Lab' in line:
            logging.debug(f'Found "Zador_Lab" line: {line}')
            path = line.strip()
            fields = path.split('/')
            for i,val in enumerate(fields):
                if val.strip() == 'Zador_Lab':
                    #print(f'with Zador_Lab = {fields}')
                    studyid = fields[i+1]
                    libid = fields[-1]
                    logging.debug(f'found libid={libid}')
        
        elif 'MBases' in line:
            logging.debug(f'Found nbases line: {line}')
            nbases = lines[idx + 1]
            logging.debug(f'nbases line = {nbases}')
            fields = nbases.split()
            nbases = int(fields[-1].replace(',',''))
            #print(f'nbases = {nbases}')            
            if nbases > NEXTSEQ_THRESHOLD:
                studytype = 'nextseq'
        #print(line)
    
    if libid is not None:
        studyid = f'LID{libid}_{studyid}'
    if expid is not None:
        studyid = f'{expid}_{studyid}'
    logging.debug(f'final studyid={studyid}')
    return [studytype, studyid, path ]

def add_lines(part):
    previous = None
    new = []
    part = part.get_payload(decode=True)
    if type(part) == bytes:
        part = part.decode(errors='replace')


    lines = str(part).split('\n')
    for line in lines:
        if len(line) > 2:
            #print(f'{line}')                
            if line.endswith('='):
                line = line[:-1]
                if previous is not None:
                    previous = f'{previous}{line}'
                else:
                    previous = line
                    
            else:
                if previous is not None:
                    new.append(f'{previous}{line}')
                    previous = None
                else:
                    new.append(line)
    emailtext  = '\n'.join(new)
    return emailtext



def eml_to_text(infile):
    '''
    pulls out body of *last* plain/text part of email. 
    fixes line wraps. 
    removes empty lines. 
    returns collapsed string.  
    
    '''
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    logging.debug(f'handling {filepath}')

    msg = email.message_from_file(open(infile))
    emailtext = ""
    for part in msg.walk():
        logging.debug(f'part content_type= {part.get_content_type()}')
        if part.get_content_type() == 'multipart/related':
            for subpart in part.walk():
                logging.debug(f'subpart content_type= {subpart.get_content_type()}')
                if subpart.get_content_type() == 'text/plain':
                    emailtext += add_lines(subpart)
                else:
                    logging.debug(f'not handling {subpart.get_content_type()}')
        elif part.get_content_type() == 'text/plain':
            emailtext += add_lines(part)
    logging.debug(emailtext)
    return emailtext


def handle_df(df, user, server, outdir, dryrun=True):
    '''
    perform server copy with values from df. 
    
    '''
    do_dryrun = dryrun
    for index, row in df.iterrows():
        expid = row['expid']
        seqtype = row['seqtype']
        srcdir = row['srcdir']
        srcdir = os.path.normpath(srcdir)
        logging.debug(f"\nuser={user}\nserver={server}\nseq={seqtype}\nexpid={expid}\nsrcdir={srcdir}")
        servercopy(srcdir, expid, seqtype, user, server, outdir, dryrun)


def servercopy(srcdir, expid, seqtype, user, server, outdir, dryrun=True):
    '''
    performs remote server copy 
    logs in as user @ server (assumes ssh agent/keys) 
    Creates target directory. 
    Copies over all files in 'basecalls' subdir to outdir.   
    https://www.geeksforgeeks.org/how-to-execute-shell-commands-in-a-remote-machine-using-python-paramiko/  
    '''
    cmd1 = f"mkdir -vp {outdir}/{expid}"
    cmd2 = f"cp -vr {srcdir}/* {outdir}/{expid}/"
    logging.info(f"ssh {user}@{server} '{cmd1}'")
    logging.info(f"ssh {user}@{server} '{cmd2}'")
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, port=22, username=user, timeout=3)
    if dryrun:
        cmd = '/usr/bin/pwd'
        (stdin, stdout, stderr) = client.exec_command(cmd)
        cmd_output = stdout.read().strip()
        logging.info(f'cmd output= {cmd_output}')
    else:
        try:
            (stdin, stdout, stderr) = client.exec_command(cmd1)
            cmd_output = stdout.read().strip()
            logging.info(f'cmd output= {cmd_output}')            
            (stdin, stdout, stderr) = client.exec_command(cmd2)
            cmd_output = stdout.read().strip()
            logging.info(f'cmd output= {cmd_output}')
        except Exception as ex:
            logging.warning('problem with ssh connection.')
            logging.warning(traceback.format_exc(None))    

                    
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

    parser.add_argument('-c', '--copyfiles', 
                        action="store_true", 
                        dest='copyfiles', 
                        default=False,
                        help='perform remote copy?')

    parser.add_argument('-o','--outdir', 
                    metavar='outdir',
                    required=False,
                    default='/grid/mbseq/data_norepl/mapseq/raw_data/nextseq', 
                    type=str, 
                    help='outdir on server. ') 

    parser.add_argument('-s','--server', 
                    metavar='server',
                    default='bamdev1.cshl.edu',
                    required=False, 
                    type=str, 
                    help='server to copy on.') 

    parser.add_argument('-e','--expid', 
                    metavar='expid',
                    required=False,
                    default=None,
                    type=str, 
                    help='explicitly provided experiment id')

    parser.add_argument('-u','--user', 
                    metavar='user',
                    default=os.getlogin(),
                    required=False, 
                    type=str, 
                    help='username to copy as.') 
    
    parser.add_argument('infile' ,
                        metavar='infile', 
                        type=str,
                        nargs='?',
                        default=None, 
                        help='.eml files from CSHL NGS')
       
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    logging.debug(f'handling {args.infile} ...')
    df = email_to_table(args.infile)
    if args.expid is not None:
        df['expid'] = args.expid    
    df.to_csv(sys.stdout, sep='\t', index=False, header=True)
    
    if args.copyfiles: 
        logging.info(f'performing copy on {args.user}@{args.server} to {args.outdir}')
        handle_df(df, args.user, args.server, args.outdir, dryrun=False) 
        
    
    