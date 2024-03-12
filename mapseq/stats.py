#!/usr/bin/env python
# mechanism to collect statistics/information during processing to be written out
# to file for later reporting. 
# default:  stats.json
#  fastq
#  ssifasta
#  merged
#
#
#

import os
import json
import logging
import sys

import datetime as dt

from pprint import pprint
from collections import defaultdict


MYSTATS = None

def def_dict_value(): 
    return defaultdict(def_dict_value)


class StatsHandler(object):
    
    def __init__(self, config, outdir=None, datestr = None ):
        global MYSTATS
        self.config = config
        if outdir is None:
            outdir = "."
        self.outdir = outdir
        self.statsdict = {}
        if datestr is None:
            self.datestr = dt.datetime.now().strftime("%Y%m%d%H%M")    
        else:
            self.datestr = datestr
        self.filename = os.path.abspath(f'{outdir}/stats.{self.datestr}.json')
        MYSTATS = self
        self.write_stats()
        
        
    def __repr__(self):
        return str(self.statsdict)

    def add_value(self, path, key, value):
        '''
        add_value'/a/b/c', key , 'myvalue') ->
        { 'a': {'b': {'c': { 'key' : 'myvalue'}}}}
        add_value('/a/b/c', key2, 'myvalue2') ->
        { 'a': {'b': {'c': { 'key' : 'myvalue', key2 : 'myvalue2'}}}}
        
        call creates path as needed. 
        '''
        pathlist = path.split('/')
        pathlist = [x for x in pathlist if x]
        logging.debug(f'pathlist = {pathlist}')
        cur_dict = self.statsdict
        for k in pathlist:
            try:
                cur_dict = cur_dict[k]
            except KeyError:
                cur_dict[k] = {}
                cur_dict = cur_dict[k]
        
        cur_dict[key] = value
        self.write_stats()
        # User has a reference to dict, but return it anyway
        return self.statsdict
    
    def get_statsdict(self):
        '''
        utility function, do not call externally for normal usage.  
        '''    
        return self.statsdict
    
    def write_stats(self):
        with open(self.filename ,'w', encoding='utf-8') as json_data:
            json.dump(self.statsdict, json_data, ensure_ascii=False, indent=4)
        logging.debug(f'wrote {self.statsdict}')


def get_default_stats():
    global MYSTATS
    if MYSTATS is not None:
        return MYSTATS
    else:
        raise Exception('Stats not initialized yet. Do so. ')
        #cp = get_default_config()
        #sh = StatsHandler(cp)
        #MYSTATS = sh
        #return sh

def make_stats_handler(config, outdir=None, datestr = None ):
    sh = StatsHandler(config, outdir, datestr )
    return sh



def old_get_stats_object():    
    try:
        with open('stats.json') as json_data:
            statsobj = json.load(json_data)
        pprint(statsobj)    
        
    except FileNotFoundError:
        logging.warning(f'no stats.json file found. creating new one.')
        statsobj = {}
        with open('stats.json','w', encoding='utf-8') as json_data:
            json.dump(statsobj, json_data, ensure_ascii=False, indent=4)
            logging.debug(f'wrote {statsobj}')
    return statsobj 

def old_save_stats_object(statsobj):
    with open('stats.json','w', encoding='utf-8') as json_data:
        json.dump(statsobj, json_data, ensure_ascii=False, indent=4)
    logging.debug(f'wrote {statsobj}')


def runtests():
    o = get_stats_object()
    pprint(o)
    o['stats'] = {}
    save_stats_object(o)
   
