#!/usr/bin/env python
#  quick and dirty script to remove spaces, caps, and special characters
#  from filenames 

import sys 
import os 
import getopt
import logging
import string

# Characters to be removed, and their substitutions
FILTERMAP = {
             "&" : "and",
             ":" : ".",
             ")" : "",
             "(" : "",
             ";" : ".",
             "," : "",
             "`" : "",
             '"' : "",
             "'" : "",
             ' ' : '_',
             '[' : "",
             ']' : "",
             '-' : ".",
             }

def filterName(name):
    #for k in FILTERMAP.keys():
    #    name = replace(name, k, FILTERMAP[k])
    #return name
    return name.translate(str.maketrans(FILTERMAP))


if __name__ == "__main__":  
    usage = '''normalizenames.py [FILES]
     -h --help       This helpful message.
     -d --debug      Debug messages.
     -v --verbose    Verbose output
'''
    debug = 0
    verbose = 0
    
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "hdv", ["help", "debug", "verbose"])
    
    except getopt.GetoptError:
        print("Unknown option...")
        print(usage)                          
        sys.exit(1)        
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(usage)                     
            sys.exit()            
        elif opt in ("-d", "--debug"):
            debug = 1
            verbose = 1
        elif opt in ("-v", "--verbose"):
            verbose = 1
            
    files = args

    logging.basicConfig()
    log = logging.getLogger()

    if debug:
        log.setLevel(logging.DEBUG)
    elif verbose :
        log.setLevel(logging.INFO)
    else:
        log.setLevel(logging.WARN)

    for f in files:
        log.debug("File: %s" % f)
        fixed = filterName(f)
        log.debug("Fixed name: %s" % fixed )
        log.info("%s -> %s" % (f,fixed))
        os.rename(f,fixed)
