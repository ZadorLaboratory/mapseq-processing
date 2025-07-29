import os
import json
import logging
import sys

import datetime as dt

from pprint import pprint
from collections import defaultdict

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import *
from mapseq.barcode import *
from mapseq.utils import *

from jinja2 import Template
import codecs
from mergedeep import merge


MYSTATS = None

def def_dict_value(): 
    return defaultdict(def_dict_value)

class StatsHandler(object):
    
    def __init__(self, outdir=None, datestr = None, outfile=None ):
        global MYSTATS
        if outdir is None:
            outdir = "./"
        self.outdir = outdir
        self.statsdict = {}
        if datestr is None:
            self.datestr = dt.datetime.now().strftime("%Y%m%d%H%M")    
        else:
            self.datestr = datestr
        if outfile is None:
            self.filename = os.path.abspath(f'{outdir}/stats.{self.datestr}.json')
        else:
            self.filename = os.path.abspath(outfile)
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
        #raise Exception('Stats not initialized yet. Do so. ')
        logging.debug(f'creating new StatsHandler')
        sh = StatsHandler()
        MYSTATS = sh
        return sh

def make_stats_handler( outdir=None, datestr = None ):
    sh = StatsHandler( outdir, datestr )
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


def config_asdict(cp):
    '''
    make configparser object into dict of dicts. 
    '''
    cdict = {}
    for s in cp.sections():
        sdict = {}
        for k in list(cp[s].keys()):
            sdict[k] = cp.get(s, k)
        cdict[s] = sdict
    return cdict
        

def generate_report( template, json_list, cp=None):
    '''
    @arg template    Markdown text file template for output doc. 
    @arg json_list    Python list of json files containing keys/values for doc.
    
    @return string of Markdown text with values interpolated.
    '''

    conf = config_asdict(cp)
    logging.debug(f'cdict = {conf}')

    stats = {}        
    for sfile in json_list:
        with open(sfile) as json_data:
            logging.debug(f'importing {sfile} ...')
            obj = json.load(json_data)
            stats = merge(stats, obj)
    
    logging.debug(f'stats={stats}')        
    
    with open(template, 'r') as file:
        template = Template(file.read(), trim_blocks=True)
        rendered_file = template.render(stats=stats, conf=conf)

        print(rendered_file)
        #output the file
        #output_file = codecs.open("report.md", "w", "utf-8")
        #output_file.write(rendered_file)
        #output_file.close()
    return rendered_file
                   
def calc_template_switch(df, cp=None):
    '''
    from full VBCtable, calculate template switch rate. 
    
    only relevant for nextseq, not novaseq

    ONLY target areas. because normally injecton are not pooled, and so swapping doesn't matter. 
    normalize by spikes by target
    
        
    load('barcodematrixL1M253_CR.mat');    % Load L1 matrix
    load('barcodematrixM253_CR.mat');      % load L2 matrix
    load('spikesM253_CR.mat');  % load spike in
    target_L2 = 6:70; %change target site info
    target_L1 = 74; %change L1 target site info
    num_spikein = zeros(1,length(spikes)); %for too many spike-in counts
    for i=1:length(spikes)
        num_spikein(i) = length(spikes(i).counts2u);
    end
    
    
    num_target_L2 = sum(sum(barcodematrix(:,target_L2)));       % total num of L2 molecules in L2 targets
    num_target_L1 = sum(sum(barcodematrixL1(:,target_L1)));         % total num of L1 molecules in L1 targets
    num_spikes_L1 = sum(num_spikein(target_L1));           % total num of spike in molecules in L1 targets
    num_spikes_L2 = sum(num_spikein(target_L2));           % total num of spike in molecules in L2 targets
    num_templateswitching = sum(sum(barcodematrix(:,target_L1)));   % num of L2 molecules detected in L1 targets
    c = num_templateswitching / (num_target_L2 * (num_spikes_L1 + num_target_L1) );        % template switching coefficient
    
    ratio_ts = 0.5 * c * num_target_L2+c * num_spikes_L2 ;         % ratio of template swtiching molecules in all L2 targets

    '''


def calc_false_positive(df, cp=None):
    '''
    from full VBCtable, calculate false positive rate. 
    
    from UMI_threshold.m 

    threshold_UMI = 0:10;       % test a variety of UMI threshold_UMI
    threshold_injection = 30;   % threshold for injection sites
    idx_injection = 1:2;        % SSI for injection sites
    idx_target = [6:70 73];     % SSI for target sites (including negative control)
    idx_negative_ctrl = 73;     % SSI for ctrl

    error_rate_false_positive = zeros(1,length(threshold_UMI));

    for i=1:length(threshold_UMI)       % calculate the error rate for each UMI threshold
        % # of neurons with false positive projection: max injection UMI>50 AND max ctrl UMI > threshold
        num_false_positive = sum( max( barcodematrix(:,idx_injection), [],2) > threshold_injection & max(barcodematrix(:,idx_negative_ctrl),[],2)>threshold_UMI(i) ); 
    
        % # of projection neurons: max injection UMI>50 AND max target UMI > threshold
        num_total = sum( max(barcodematrix(:,idx_injection),[],2) > threshold_injection & max(barcodematrix(:,idx_target),[],2)>threshold_UMI(i) );
    
        error_rate_false_positive(i) = num_false_positive/num_total;
    end    
    
    '''    
    if cp is None:
        cp = get_default_config()
        
#
#        QC/ Assessment routines. 
#
        
def assess_readinfo(df,
                    outdir=None, 
                    cp=None):
    ''' 
    Gather info/classify quality of sequencing read info/ fields. 
    Input is split sequence from aggregated reads. 
    
    -- L1/L2 vs spike-in
    -- distribution of rtag/rrtag variation.  
    
    '''
    # parameters/inputs
    spikeseq = cp.get('readtable','spikeseq')
    realregex = cp.get('readtable', 'realregex' )
    loneregex = cp.get('readtable', 'loneregex' )
    use_libtag = cp.getboolean('readtable','use_libtag')
    bcfile = os.path.expanduser( cp.get('barcodes','ssifile') )
    inj_min_reads = int(cp.get('vbctable','inj_min_reads'))
    target_min_reads = int(cp.get('vbctable','target_min_reads'))

    use_libtag = True    

    sh = get_default_stats()

    #sampdf = load_sample_info(sinfo, cp=cp)
    #sinfo = 'M297.sampleinfo.xlsx'
    #sampdf = load_sample_info(sinfo, cp=cp)
    #sampdf
    #labels = get_rtlist(sampdf)
    #    bcdict = get_barcode_dict(bcfile, labels)
    #    rtdict = get_rtprimer_dict(bcfile, labels)
    #    rtbcdict = get_rtbc_dict(bcfile, labels)
    #labels = get_rtlist(sampdf)
    #bcdict = get_barcode_dict(bcfile, labels)
    #rtdict = get_rtprimer_dict(bcfile, labels)
    #rtbcdict = get_rtbc_dict(bcfile, labels)
    
    # Calculate matches/non-matches for 
    smap = df['spikeseq'] == spikeseq
    df.loc[smap, 'type'] = 'spike'
    #df['type'].value_counts()
    sh.add_value('/readinfo/all', 'n_spikes', str(smap.sum() ) )
    rmap = df['libtag'].str.match(realregex)
    sh.add_value('/readinfo/all', 'n_real', str( rmap.sum() ) )
    lmap = df['libtag'].str.match(loneregex)
    sh.add_value('/readinfo/all', 'n_lones', str( lmap.sum() ) )
    df.loc[rmap, 'type'] = 'real'
    df.loc[lmap, 'type'] = 'lone'
    namap = df['type'].isna()
    sh.add_value('/readinfo/all', 'n_lones', str(  namap.sum() ) )    
    sh.add_value('/readinfo', 'inj_min_reads', str(inj_min_reads ) )
    sh.add_value('/readinfo', 'target_min_reads', str(target_min_reads ) )

    tdf = df[df['site'].str.startswith('target')]
    tdf = tdf[tdf['read_count'] >= int(target_min_reads)]
    idf = df[df['site'].str.startswith('injection')]
    idf = idf[idf['read_count'] >= int(inj_min_reads)]    
    df = pd.concat([tdf, idf])
    df.reset_index(drop=True, inplace=True)
    
    
    
    # Calculate all again after read thresholding.     
    sh.add_value('/readinfo/readthresholded', 'n_total', str(len(df) ) )
    smap = df['spikeseq'] == spikeseq
    df.loc[smap, 'type'] = 'spike'
    #df['type'].value_counts()
    sh.add_value('/readinfo/readthresholded', 'n_spikes', str(smap.sum() ) )
    rmap = df['libtag'].str.match(realregex)
    sh.add_value('/readinfo/readthresholded', 'n_real', str( rmap.sum() ) )
    lmap = df['libtag'].str.match(loneregex)
    sh.add_value('/readinfo/readthresholded', 'n_lones', str( lmap.sum() ) )
    df.loc[rmap, 'type'] = 'real'
    df.loc[lmap, 'type'] = 'lone'
    namap = df['type'].isna()
    sh.add_value('/readinfo/all', 'n_lones', str(  namap.sum() ) )    
    print(sh)


def merge(a: dict, b: dict, path=[]):
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                merge(a[key], b[key], path + [str(key)])
            elif a[key] != b[key]:
                raise Exception('Conflict at ' + '.'.join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a

def merge_stats(infiles, outfile):
    '''
    combine multiple stats files into one unified (additive) JSON tree.
    needed because some processes put entries into the same top-level dict, 
    e.g. 'fastq' from both process_fastq_pairs, and filter_reads
     
    '''
    merged = None
    for infile in infiles:
        logging.debug(f'reading {infile}')
        with open(infile) as jfh:
            jo = json.load(jfh)
        if merged is None:
            logging.debug(f'merged is None. setting.')
            merged = jo
        else:
            logging.debug(f'merged exists. merging...')
            merged = merge(merged, jo)
    logging.debug(f'finished merging. writing to {outfile}')
    with open(outfile ,'w', encoding='utf-8') as jofh:
        json.dump(merged, jofh, ensure_ascii=False, indent=4)
    logging.debug('done.')


def qc_make_readmatrix( df, sampdf=None, outdir='./', cp=None):
    '''
    make matrix of SSI vs. FASTQ  
    
    '''
    logging.debug(f'inbound df len={len(df)} columns={df.columns}')
    sh = StatsHandler(outdir=outdir)

    if cp is None:
        cp = get_default_config()
    bcfile = os.path.expanduser( cp.get('barcodes','ssifile') ) 

    # Map label, rtprimer to SSIs    
    logging.debug(f'getting rt labels...')
    labels = None
    if sampdf is not None:
        labels = get_rtlist(sampdf)
    n_reads = len(df)
    sh.add_value('/qcreads/sources', 'n_reads', str(n_reads ) )
    logging.debug(f'rtlabels={labels}')
    bcdict = get_barcode_dict(bcfile, labels)
    n_ssis =  len(bcdict)
    sh.add_value('/qcreads/sources', 'n_ssis', str(n_ssis) )
    rtdict = get_rtprimer_dict(bcfile, labels)
    rtbcdict = get_rtbc_dict(bcfile, labels)
    logging.debug(f'got {len(bcdict)} barcodes with labels and primer number.')    

    logging.info('filling in rtprimer number by SSI sequence...')    
    df['rtprimer'] = df['ssi'].map(rtdict)
    df['rtprimer'] = df['rtprimer'].astype('category')
    
    logging.info('filling in label by rtprimer...')
    df['label'] = df['rtprimer'].map(rtbcdict)
    df['label'] = df['label'].astype('category')
    
    
    gdf = df.dropna()
    gdf.reset_index(inplace=True, drop=True)
    n_readswssi =  len(gdf)
    sh.add_value('/qcreads/sources', 'n_readswssi', str(n_readswssi) )
    good_pct = ( n_readswssi / n_reads ) * 100
    good_spct = f'{good_pct:.2f}'
    sh.add_value('/qcreads/sources', 'n_reads_valid_ssi', str(good_spct) )
    
    gadf = gdf.groupby(by=['label','source'], observed=True).agg( {'source':'count'})
    gadf.columns = ['count']
    gadf.reset_index(inplace=True,drop=False)
    
    sldf = gadf.pivot(index='source',columns='label',values='count')
    sldf.fillna(0.0, inplace=True)
    scol = natsorted(list(sldf.columns))
    sldf = sldf[scol]
    nsindex = natsorted(sldf.index)
    sldf = sldf.loc[nsindex]
    
    outfile = os.path.join(outdir, 'EXP.source_reads.xlsx')
    with pd.ExcelWriter(outfile) as writer:
        sldf.to_excel(writer, sheet_name='source by SSI')
    return sldf
    
    


