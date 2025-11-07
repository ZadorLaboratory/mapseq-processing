#!/usr/bin/env python 



import itertools
import logging
import os
import subprocess
import sys
import traceback

import datetime as dt

from configparser import ConfigParser
from collections import defaultdict

import networkx as nx
import numpy as np
import pandas as pd

from Bio import SeqIO

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import *
from mapseq.stats import *


def align_collapse_pd(df,
                      column='vbc_read',
                      pcolumn='read_count', 
                      max_mismatch=None,
                      max_recursion=None, 
                      outdir=None, 
                      datestr=None,
                      force=False,
                      min_reads = None,
                      drop=True, 
                      cp=None):
    '''
    Assumes dataframe with sequence and read_count columns
    Use read_count to choose parent sequence.
    max_recursion appears to be important, rather than memory reqs. 40000 or 50000 
        novaseq (1.5B reads) may be needed. 
    By default, picking up wherever the algorithm left off, depending on
        file availability. Force overrides and recalculates from beginning. 
      
    TODO:
    Speed up use of pandas, map() function rather than apply()
    relationship between read_count on full_read sequences (52nt), and 
        value_counts() on unique VBCs (30nts)
    
    '''
    # housekeeping...
    if cp is None:
        cp = get_default_config()
    
    aligner = cp.get('collapse','tool')    
    if max_mismatch is None:
        max_mismatch = int(cp.get('collapse', 'max_mismatch'))
    else:
        max_mismatch = int(max_mismatch)

    if min_reads is None:
        min_reads = int(cp.get('collapse', 'min_reads'))
    else:
        min_reads = int(min_reads)

    if max_recursion is not None:
        rlimit = int(max_recursion)
        logging.info(f'set new recursionlimit={rlimit}')
        sys.setrecursionlimit(rlimit)
    else:
        rlimit = int(cp.get('collapse','max_recursion'))
        logging.info(f'set new recursionlimit={rlimit}')
        sys.setrecursionlimit(rlimit)        

    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    if outdir is None:
        outdir = './'
    
    outdir = os.path.abspath(outdir)    
    os.makedirs(outdir, exist_ok=True)
    logging.debug(f'collapse: aligner={aligner} max_mismatch={max_mismatch} outdir={outdir}')    
    
    sh = get_default_stats()      
    sh.add_value('/collapse','n_full_sequences', len(df) )
    sh.add_value('/collapse','max_mismatch_call', str(max_mismatch) )
    logging.info('Getting unique DF...')    
    udf = df[column].value_counts().reset_index() 
    sh.add_value('/collapse','n_unique_sequences', len(udf) )    

    # handle writing unique fasta
    of = os.path.join( outdir , f'{column}.unique.tsv')  
    if os.path.exists(of) and not force:
        logging.debug(f'Output {of} exists and not force.')
    else:
        logging.debug(f'Writing uniques to {of}')
        udf.to_csv(of, sep='\t')
        
    of = os.path.join( outdir , f'{column}.unique.fasta')      
    if os.path.exists(of) and not force:
        logging.debug(f'Output {of} exists and not force.')
        seqfasta = of
    else: 
        logging.info(f'Writing uniques as FASTA to {of}')
        seqfasta = write_fasta_from_df(udf, outfile=of, sequence=[column])    

    # run allXall bowtiex
    of = os.path.join( outdir , f'unique_sequences.bt2.sam')
    if os.path.exists(of) and not force:
        logging.debug(f'Output {of} exists and not force.')
        btfile = of
    else:
        logging.info(f'Running {aligner} on {seqfasta} file to {of}')
        btfile = run_bowtie(cp, seqfasta, of, tool=aligner)
    
    logging.info(f'Bowtie done. Produced output {btfile}. Checking DF.')
    of = os.path.join( outdir , f'unique_sequences.btdf.tsv') 
    if os.path.exists(of) and not force:
        logging.debug(f'Output {of} exists and not force.')
        btdf = load_bowtie_df(of)
    else:
        logging.info(f'Creating btdf dataframe from {btfile} max_mismatch={max_mismatch}')    
        btdf = make_bowtie_df(btfile, max_mismatch=max_mismatch, ignore_self=True)
        logging.debug(f'Writing output to {of}')
        btdf.to_csv(of, sep='\t') 
    sh.add_value('/collapse','n_bowtie_entries', len(btdf) )

    # handle edgelist 
    of = os.path.join( outdir , f'edgelist.txt')     
    if os.path.exists(of) and not force:
        logging.debug(f'Output {of} exists and not force.')
        edgelist = read_listlist(of)
    else:
        edgelist = edges_from_btdf(btdf)
        writelist(of, edgelist)   
    btdf = None  # help memory usage
    logging.debug(f'edgelist len={len(edgelist)}')
    sh.add_value('/collapse','n_edges', len(edgelist) )

    # handle components.     
    of = os.path.join( outdir , f'components.txt')
    if os.path.exists(of) and not force:
        logging.debug(f'Output {of} exists and not force.')
        components = read_listlist(of)
    else:
        logging.info('Calculating Hamming components...')    
        components = get_components(edgelist)
        writelist(of, components)
    edgelist = None  # help memory usage
    logging.debug(f'all components len={len(components)}')
    sh.add_value('/collapse','n_components', len(components) )
    
    # assess components...
    # components is list of lists.
    data = [ len(c) for c in components]
    data.sort(reverse=True)
    ccount = pd.Series(data)
    of = os.path.join( outdir , f'component_count.tsv')
    ccount.to_csv(of, sep='\t')            

    mcomponents = remove_singletons(components)
    logging.debug(f'multi-element components len={len(mcomponents)}')
    sh.add_value('/collapse','n_multi_components', len(mcomponents) )
    of = os.path.join( outdir , f'multi_components.json')
    logging.debug(f'writing components len={len(components)} t {of}')
    with open(of, 'w') as fp:
        json.dump(mcomponents, fp, indent=4)

    logging.info(f'Collapsing {len(components)} components...')
    newdf = collapse_by_components_pd(df, udf, mcomponents, column=column, pcolumn=pcolumn, outdir=outdir)
    newcol = f'{column}_col'        
    newdf.rename( {'newsequence': newcol}, inplace=True, axis=1)    
    logging.info(f'Got collapsed DF. len={len(newdf)}')

    rdf = newdf[[ column, newcol, pcolumn ]]
    of = os.path.join( outdir , f'read_collapsed.tsv')
    logging.info(f'Writing reduced mapping TSV to {of}')
    rdf.to_csv(of, sep='\t')
    if drop:
        newdf.drop(column, inplace=True, axis=1)
        newdf.rename( { newcol : column }, inplace=True, axis=1)       
    return newdf        


def align_collapse_pd_grouped(df,
                              column='vbc_read',
                              pcolumn='read_count',
                              gcolumn='brain', 
                              max_mismatch = None,
                              max_recursion = None, 
                              outdir= None, 
                              datestr= None,
                              force= False,
                              min_reads = None, 
                              cp=None):
    '''
    Groups alignment and collapse by gcolumn value. [brain]
    Uses sub-directories for standard intermediate output/scratch. 
    Assumes dataframe with sequence (vbc_read) and read_count columns
    Use read_count to choose parent sequence.
        
    max_recursion appears to be important, rather than memory reqs. 40000 or 50000 
        novaseq (1.5B reads) may be needed. 
    By default, picking up wherever the algorithm left off, depending on
        file availability. Force overrides and recalculates from beginning. 
      
    TODO:
    Speed up use of pandas, map() function rather than apply()
    relationship between read_count on full_read sequences (52nt), and 
        value_counts() on unique VBCs (30nts)
    
    '''
    # housekeeping...
    if cp is None:
        cp = get_default_config()
        
    if max_mismatch is None:
        max_mismatch = int(cp.get('collapse', 'max_mismatch'))
    else:
        max_mismatch = int(max_mismatch)

    if min_reads is None:
        min_reads = int(cp.get('collapse', 'min_reads', fallback=1))
    else:
        min_reads = int(min_reads)

    if max_recursion is not None:
        rlimit = int(max_recursion)
        logging.info(f'set new recursionlimit={rlimit}')
        sys.setrecursionlimit(rlimit)
    else:
        rlimit = int(cp.get('collapse','max_recursion'))
        logging.info(f'(from config) set new recursionlimit={rlimit}')
        sys.setrecursionlimit(rlimit)        

    aligner = cp.get('collapse','tool')

    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    
    if outdir is None:
        outdir = './'
    outdir = os.path.abspath(outdir)    
    os.makedirs(outdir, exist_ok=True)

    logging.debug(f'collapse: aligner={aligner} max_mismatch={max_mismatch} outdir={outdir}')    
    sh = get_default_stats()
    
    sh.add_value('/collapse','n_initial_sequences', len(df) )
    if min_reads > 1:
        logging.info(f'thresholding by min_reads={min_reads} before collapse. len={len(df)} ')
        df = df[ df['read_count'] >= min_reads  ]
        df.reset_index(inplace=True, drop=True)
        logging.info(f'after min_reads={min_reads} filter. len={len(df)} ')
    sh.add_value('/collapse','n_thresholded_sequences', len(df) )    
    sh.add_value('/collapse','n_full_sequences', len(df) )
    sh.add_value(f'/collapse','api_max_mismatch', str(max_mismatch) )

    # Get list of ids to group collapse by...
    df[gcolumn] = df[gcolumn].astype('string')
    gidlist = list( df[gcolumn].dropna().unique() )
    gidlist = [ x for x in gidlist if len(x) > 0 ]
    gidlist.sort()
    logging.debug(f'handling group list: {gidlist}')

    logging.info(f'pulling out no group for merging at end...')
    nogroup_df = df[ df[gcolumn] == '' ]
    sh.add_value('/collapse','n_no_group', len(nogroup_df) )

    gdflist = []

    for gid in gidlist:
        logging.info(f"collapsing '{column}' by {gcolumn} = '{gid}' ")
        gdir = os.path.join( outdir, f'{gcolumn}.{gid}' )
        os.makedirs(gdir, exist_ok=True)
        
        gdf = df[df[gcolumn] == gid]
        gdf.reset_index(inplace=True, drop=True)
        initial_len = len(gdf)  
        logging.info(f'[{gcolumn}:{gid}] initial len={len(gdf)} subdir={gdir}')        
        
        # get reduced dataframe of unique head sequences
        logging.info('Getting unique DF...')    
        udf = gdf[column].value_counts().reset_index() 
        sh.add_value(f'/collapse/{gcolumn}_{gid}','n_unique_sequences', len(udf) )    
    
        of = os.path.join( gdir , f'{column}.unique.tsv')
        logging.info(f'Writing unique DF to {of}')
        udf.to_csv(of, sep='\t') 
        
        # handle writing unique fasta    
        of = os.path.join( gdir , f'{column}.unique.fasta')      
        if os.path.exists(of) and not force:
            logging.debug(f'Output {of} exists and not force.')
            seqfasta = of
        else: 
            logging.info(f'Writing uniques as FASTA to {of}')
            seqfasta = write_fasta_from_df(udf, outfile=of, sequence=[column])        
     
        # run allXall bowtiex
        of = os.path.join( gdir , f'unique_sequences.bt2.sam')
        if os.path.exists(of) and not force:
            logging.debug(f'Output {of} exists and not force.')
            btfile = of
        else:
            logging.info(f'Running {aligner} on {seqfasta} file to {of}')
            btfile = run_bowtie(cp, seqfasta, of, tool=aligner) 
            
        logging.info(f'Bowtie done. Produced output {btfile}. Checking DF.')
        of = os.path.join( gdir , f'unique_sequences.btdf.tsv')
        if os.path.exists(of) and not force:
            logging.debug(f'Output {of} exists and not force.')
            btdf = load_bowtie_df(of)
        else:
            logging.info(f'Creating btdf dataframe from {btfile} max_mismatch={max_mismatch}')    
            btdf = make_bowtie_df(btfile, 
                                  max_distance=max_mismatch, 
                                  ignore_self=True)
            logging.debug(f'Writing output to {of}')
            btdf.to_csv(of, sep='\t')
        sh.add_value(f'/collapse/{gcolumn}_{gid}','n_bowtie_entries', len(btdf) )
    
        # handle edgelist 
        of = os.path.join( gdir , f'edgelist.txt')     
        if os.path.exists(of) and not force:
            logging.debug(f'Output {of} exists and not force.')
            edgelist = read_listlist(of)
        else:
            edgelist = edges_from_btdf(btdf)
            writelist(of, edgelist)   
        btdf = None  # help memory usage
        logging.debug(f'edgelist len={len(edgelist)}')
        sh.add_value(f'/collapse/{gcolumn}_{gid}','n_edges', len(edgelist) )

        # handle components.     
        of = os.path.join( gdir , f'components.txt')
        if os.path.exists(of) and not force:
            logging.debug(f'Output {of} exists and not force.')
            components = read_listlist(of)
        else:
            logging.info('Calculating Hamming components...')    
            components = get_components(edgelist)
            writelist(of, components)
        edgelist = None  # help memory usage
        logging.debug(f'all components len={len(components)}')
        
        # assess components...
        # components is list of lists.
        data = [ len(c) for c in components]
        data.sort(reverse=True)
        ccount = pd.Series(data)
        of = os.path.join( gdir , f'component_count.tsv')
        ccount.to_csv(of, sep='\t')            
    
        mcomponents = remove_singletons(components)
        logging.debug(f'multi-element components len={len(mcomponents)}')
        sh.add_value(f'/collapse/{gcolumn}_{gid}','n_multi_components', len(mcomponents) )
        of = os.path.join( gdir , f'multi_components.json')
        logging.debug(f'writing components len={len(components)} t {of}')
        with open(of, 'w') as fp:
            json.dump(mcomponents, fp, indent=4)
    
        logging.info(f'Collapsing {len(components)} components...')
        newdf = collapse_by_components_pd(gdf, udf, mcomponents, column=column, pcolumn=pcolumn, outdir=gdir)
        # newdf has sequence and newsequence columns, rename to orig_seq and sequence
        newcol = f'{column}_col'        
        newdf.rename( {'newsequence': newcol}, inplace=True, axis=1)    
        logging.info(f'Got collapsed DF. len={len(newdf)}')
    
        rdf = newdf[[ column, newcol, pcolumn ]]
        of = os.path.join( gdir , f'read_collapsed.tsv')
        logging.info(f'Writing reduced mapping TSV to {of}')
        rdf.to_csv(of, sep='\t')
        
        newdf.drop(column, inplace=True, axis=1)
        newdf.rename( { newcol : column }, inplace=True, axis=1) 
        gdflist.append(newdf)       
    
    # merge all brains into one dataframe...
    logging.debug(f'sizes={[ len(x) for x in gdflist ]} adding nogroup len={len(nogroup_df)}')    
    gdflist.append(nogroup_df)
    outdf = pd.concat(gdflist, ignore_index = True)
    outdf.reset_index(inplace=True, drop=True)
    logging.info(f'All groups. Final DF len={len(outdf)}')
    return outdf

#
#  BOWTIE/TARJAN FUNCTIONS
#
def edges_from_btdf(btdf):
    readlist = btdf.name_read.values.tolist()
    alignlist = btdf.name_align.values.tolist()  
    edgelist = [ list(t) for t in zip(readlist, alignlist)]
    return edgelist


def get_components(edgelist, integers=True):
    '''
    returns strongly connected components from list of edges via Tarjan's algorithm. 
    assumes labels are integers (for later use as indices in dataframes. 
    '''
    complist = []
    logging.debug(f'getting connected components from edgelist len={len(edgelist)}')
    if len(edgelist) < 100:
        logging.debug(f'{edgelist}')
    for g in tarjan(from_edges(edgelist)):
        #logging.debug(f'g={g}')
        complist.append(g)
    logging.debug(f'{len(complist)} components.')
    if len(edgelist) < 100:
        logging.debug(f'{complist}')
    if integers:
        outlist = []
        for g in complist:
            outlist.append( [int(x) for x in g])
        complist = outlist    
    return complist

#
#  Tarjan's algorithm, same as Matlab graphconncomp()   
#  https://rosettacode.org/wiki/Tarjan#Python:_As_function
#
def from_edges(edges):        
    class Node:
        def __init__(self):
            # root is one of:
            #   None: not yet visited
            #   -1: already processed
            #   non-negative integer: what Wikipedia pseudo code calls 'lowlink'
            self.root = None
            self.succ = []

    nodes = defaultdict(Node)
    for v,w in edges:
        nodes[v].succ.append(nodes[w])

    for i,v in nodes.items(): # name the nodes for final output
        v.id = i

    return nodes.values()
    
    
def tarjan(V):
    '''
    May get recursion limit errors if input is large. 
    https://stackoverflow.com/questions/5061582/setting-stacksize-in-a-python-script/16248113#16248113
    
    import resource, sys
    resource.setrlimit(resource.RLIMIT_STACK, (2**29,-1))
    sys.setrecursionlimit(10**6)
    
    Same algorithm as used in MATLAB pipeline. 
    https://www.mathworks.com/help/matlab/ref/graph.conncomp.html 
    
    Non-recursive version? 
    https://stackoverflow.com/questions/46511682/non-recursive-version-of-tarjans-algorithm
    
    https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm 
    
    '''
    def strongconnect(v, S): 
        v.root = pos = len(S)
        S.append(v)
        for w in v.succ:
            if w.root is None:  # not yet visited
                yield from strongconnect(w, S)

            if w.root >= 0:  # still on stack
                v.root = min(v.root, w.root)
        if v.root == pos:  # v is the root, return everything above
            res, S[pos:] = S[pos:], []
            for w in res:
                w.root = -1
            yield [r.id for r in res]
    
    for v in V:
        if v.root is None:
            yield from strongconnect(v, [])

       

def build_seqmapdict(udf, components, column='vbc_read'):
    '''
    Create mappings from all unique sequences to component sequence
    '''
    seqmapdict = {}
    comphandled = 0
    comphandled_interval = 100000
    comp_len = len(components)
    
    for comp in components:
        ser = list(udf[column].iloc[comp])
        t = ser[0]
        for s in ser: 
            seqmapdict[s]= t    
        if comphandled % comphandled_interval == 0:
            logging.debug(f'[{comphandled}/{comp_len}] t = {t}')
        comphandled += 1
    return seqmapdict


def build_seqmapdict_pd(udf, components, column='vbc_read', pcolumn='count'):
    '''
    Create mappings from all unique sequences to component sequence
    Can we do this faster? Choose most abundant variant?
    dict should be oldsequence -> newsequence
    
    '''
    seqmapdict = {}
    comphandled = 0
    pcolumn='count'
    logging.debug(f'udf len={len(udf)} components len={len(components)} column={column} pcolumn={pcolumn} ')
    comphandled_interval = 1000
    comp_len = len(components)    
    for i, indexlist in enumerate( components):
        cdf = udf[[column, pcolumn]].iloc[indexlist]
        cdf.reset_index(inplace=True, drop=True)
        #logging.debug(f'component [{i}/{comp_len}]: len={len(cdf)}')
        maxid = int(cdf[pcolumn].idxmax())
        t = cdf[column].iloc[maxid]
        for compseq in list(cdf[column]): 
            #logging.debug(f'compseq={compseq} -> {t}')
            seqmapdict[compseq]= t    
        if comphandled % comphandled_interval == 0:
            logging.debug(f'[{i}/{comp_len}]: len={len(cdf)} seq = {t} ')
        comphandled += 1
    return seqmapdict


def collapse_by_components_pd(fulldf, 
                              uniqdf, 
                              components, 
                              column, 
                              pcolumn, 
                              outdir=None):
    #
    # *** Assumes components are multi-element components only ***
    #
    # components consist of indices within the uniqdf. 
    # assumes component elements are integers, as they are used as dataframe indices. 
    # create map of indices in original full DF that correspond to all members of a (multi-)component
    # from alignment
    # hash key is index of first appearance of sequence in full DF (i.e. its index in uniqdf.) 
    # 
    # returns copy of input fulldf with sequence column collapsed to most-common sequence in component. 
    #
    #  SETUP
    #  infile is all reads...
    #      fdf = read_fasta_to_df(infile, seq_length=32)
    #      udf = pd.DataFrame(fdf['sequence'].unique(), columns=['sequence'])
    #      components = remove_singletons(components)
    
    logging.debug(f'multi-element components len={len(components)}')
    logging.debug(f'fulldf length={len(fulldf)} uniqdf length={len(uniqdf)} {len(components)} components.')
    logging.info(f'building seqmapdict {len(uniqdf)} unique seqs, {len(components)} components, for {len(fulldf)} raw sequences. ')
    smd = build_seqmapdict_pd(uniqdf, components, column, pcolumn)
      
    # Make new full df:
    logging.info('seqmapdict built. Applying.')
    if outdir is not None:
        outfile = f'{outdir}/seqmapdict.json'
        logging.debug(f'writing seqmapdict len={len(smd)} tp {outfile}')
        with open(outfile, 'w') as fp:
            json.dump(smd, fp, indent=4)
    else:
        logging.debug(f'no outdir given.')
    logging.info(f'applying seqmapdict...')
    # make deep copy of original sequence column
    fulldf.loc[:, f'{column}_col'] = fulldf.loc[:, column]    
    # map old to new
    fulldf[f'{column}_col'] = fulldf[column].map(smd, na_action='ignore')
    # fill in NaN with original values. 
    
    fulldf.fillna( { 'vbc_read_col' : fulldf[column] }, inplace=True)    
    logging.info(f'New collapsed df = \n{fulldf}')
    log_objectinfo(fulldf, 'fulldf')
    return fulldf
      

def collapse_only_pd(fdf, udf, mcomponents, column, pcolumn, outdir, cp ):
    '''
    Handle partially completed align-collapse. 
    Inputs:
       fulldf TSV
       uniqudf TSV
       multi-components list. 
    
    Output:
        collapsed DF. 
        df = collapse_pd(fulldf,
                     uniquedf,
                     mcomponents,   
                     column=args.column,
                     pcolumn=args.parent_column,
                     outdir=outdir, 
                     cp=cp)
    '''
  
    # assess components...
    sh = get_default_stats()
    logging.debug(f'multi-element components len={len(mcomponents)}')
    sh.add_value('/collapse','n_multi_components', len(mcomponents) )

    logging.info(f'Collapsing {len(mcomponents)} mcomponents...')
    newdf = collapse_by_components_pd(fdf, udf, mcomponents, column=column, pcolumn=pcolumn, outdir=outdir)

    # newdf has sequence and newsequence columns, rename to orig_seq and sequence
    newcol = f'{column}_col'        
    newdf.rename( {'newsequence': newcol}, inplace=True, axis=1)    
    logging.info(f'Got collapsed DF. len={len(newdf)}')

    rdf = newdf[[ column, newcol, 'read_count' ]]
    of = os.path.join( outdir , f'read_collapsed.tsv')
    logging.info(f'Writing reduced mapping TSV to {of}')
    rdf.to_csv(of, sep='\t')
    logging.debug(f'dropping {column}')
    newdf.drop(column, inplace=True, axis=1)   
    return newdf        



#
#  Non-recursive version--just replaces stack with explicit structure.  
#
#  https://stackoverflow.com/questions/46511682/non-recursive-version-of-tarjans-algorithm 
#
class Node:
    def __init__(self, name):
        self.name = name
        self.index = None
        self.lowlink = None
        self.adj = []
        self.on_stack = False

def tarjans_nostack(vertices, edgelist):
    N = len(vertices)  # number of vertices
    es = edgelist      # list of edges, [(0,1), (2,4), ...]
    
    vs = [Node(i) for i in range(N)]
    for v, w in es:
        vs[v].adj.append(vs[w])
    
    i = 0
    stack = []
    call_stack = []
    comps = []
    for v in vs:
        if v.index is None:
            call_stack.append((v,0))
            while call_stack:
                v, pi = call_stack.pop()
                # If this is first time we see v
                if pi == 0:
                    v.index = i
                    v.lowlink = i
                    i += 1
                    stack.append(v)
                    v.on_stack = True
                # If we just recursed on something
                if pi > 0:
                    prev = v.adj[pi-1]
                    v.lowlink = min(v.lowlink, prev.lowlink)
                # Find the next thing to recurse on
                while pi < len(v.adj) and v.adj[pi].index is not None:
                    w = v.adj[pi]
                    if w.on_stack:
                        v.lowlink = min(v.lowlink, w.index)
                    pi += 1
                # If we found something with index=None, recurse
                if pi < len(v.adj):
                    w = v.adj[pi]
                    call_stack.append((v,pi+1))
                    call_stack.append((w,0))
                    continue
                # If v is the root of a connected component
                if v.lowlink == v.index:
                    comp = []
                    while True:
                        w = stack.pop()
                        w.on_stack = False
                        comp.append(w.name)
                        if w is v:
                            break
                    comps.append(comp)



"""
   Tarjan's algorithm and topological sorting implementation in Python
   
   by Paul Harrison
   Public domain, do with it as you will
   
   https://www.logarithmic.net/pfh-files/blog/01208083168/sort.py
   Tarjan's considered a sort algorithm because the order in which the strongly 
   connected components are identified constitutes a reverse topological sort of 
   the DAG formed by the components.   
   
"""

def strongly_connected_components(graph):
    """ Find the strongly connected components in a graph using
        Tarjan's algorithm.
        
        graph should be a dictionary mapping node names to
        lists of successor nodes.
        
        Also recursive:  visit()
    """
    
    result = [ ]
    stack = [ ]
    low = { }
        
    def visit(node):
        if node in low: 
            return
    
        num = len(low)
        low[node] = num
        stack_pos = len(stack)
        stack.append(node)
        for successor in graph[node]:
            visit(successor)
            low[node] = min(low[node], low[successor])
        
        if num == low[node]:
            component = tuple(stack[stack_pos:])
            del stack[stack_pos:]
            result.append(component)
            for item in component:
                low[item] = len(graph)
    
        for node in graph:
            visit(node)
        
        return result


def topological_sort(graph):
    count = { }
    for node in graph:
        count[node] = 0
    for node in graph:
        for successor in graph[node]:
            count[successor] += 1

    ready = [ node for node in graph if count[node] == 0 ]
    
    result = [ ]
    while ready:
        node = ready.pop(-1)
        result.append(node)
        
        for successor in graph[node]:
            count[successor] -= 1
            if count[successor] == 0:
                ready.append(successor)
    return result


def robust_topological_sort(graph):
    """ First identify strongly connected components,
        then perform a topological sort on these components. 
    """
    components = strongly_connected_components(graph)

    node_component = { }
    for component in components:
        for node in component:
            node_component[node] = component

    component_graph = { }
    for component in components:
        component_graph[component] = [ ]
    
    for node in graph:
        node_c = node_component[node]
        for successor in graph[node]:
            successor_c = node_component[successor]
            if node_c != successor_c:
                component_graph[node_c].append(successor_c) 

    return topological_sort(component_graph)


def read_listlist(infile):
    '''
    Reads a file, assumes each line is a Python list. 
    Returns list of lists. 
    
    '''
    logging.debug(f'reading file: {infile}')
    flist = []
    try:
        with open(infile, 'r') as f:
            for line in f:
                line = line.strip()
                li = eval(line)
                flist.append(li)
            logging.debug(f'got list with {len(flist)} items.')
        return flist
    except:
        logging.warning(f'error reading{infile} ')
        return []

def max_hamming(seq_list, max_ok=3):
    '''
    Calculate naive max hamming for set of sequences. 
    Keep track of how many exceed a given max. 
    
    '''
    hmax = 0
    n_pairs = 0
    n_exceed = 0
    for a,b in itertools.combinations(seq_list, 2):
        n_pairs += 1
        hd = calc_hamming(a,b, use_rc=False)
        if hd > max_ok:
            n_exceed += 1
        if hd > hmax:
            hmax = hd
    logging.debug(f'considered {n_pairs} pairs, {n_exceed} were > {max_ok}')
    return ( hmax, n_pairs, n_exceed, max_ok)       
    
    
def calc_hamming( a, b, use_rc=False):
    '''
    simple naive Hamming calculation.  
    Optionally check against reverse complement, and return minimum.
    use_rc is very slow. do not use in production. MAPseq only uses forward
    reads anyway--so no purpose. 
    '''
    mismatch = 0
    for i,c in enumerate(a):
        if c != b[i]:
            mismatch += 1

    if use_rc:
        b = Seq(b)
        rb = str( b.reverse_complement())
                
        mismatch_rc = 0
        for i,c in enumerate(a):
            if c != rb[i]:
                mismatch_rc += 1
        mismatch = min(mismatch, mismatch_rc)
    return mismatch


def make_nxgraph_seqlist( seq_list, max_mismatch = 3 ):
    '''
    Given list of sequences expected to form a Hamming 
    graph, create graph and return. 
    '''
    G = nx.Graph()
    for a,b in itertools.combinations(seq_list, 2):
        hd = calc_hamming(a,b)
        if hd <= max_mismatch :
            G.add_edge(a,b)
    logging.debug(f'made Graph: nodes={len(G.nodes)} edges={len(G.edges)} ')
    return G


def make_component_df_nx(components, parent_graph, outdir=None):
    '''
    Gather component info from native NX objects. Save to component info DF.
    Optionally write-out intermediate results since this can be a 
    very long-running operation. 
    
    '''
    COMP_INTERVAL = 10
    # write out intermediate results if outdir is not None
    # interval determined by total number of elements on all components. 
    TOTAL_ELEMENTS =  len(flatten_list(components))
    WRITE_INTERVAL = int( TOTAL_ELEMENTS / 10 )
    #   component info for DF:  cmp_idx, cmp_size, cmp_diam, cmp_max_deg 
    COMPINFO_COLUMNS = [ 'cmp_idx', 'size', 'diam', 'max_deg', 'max_deg_seq', 'max_count', 'max_count_seq']

    comps_handled = 0
    elements_handled = 0
    elements_floor = WRITE_INTERVAL
    intermediates_handled = 0
    comp_info_list = []
    intermediate_list = []        

    logging.debug(f'enumerating components...')
    if outdir is not None:
        logging.debug(f'will write out every {WRITE_INTERVAL} elements.')
    else:
        logging.debug('no intermediate output saves.')
    
    for i, comp in enumerate( components ):
        comp = list(comp)
        if comps_handled % COMP_INTERVAL == 0:
            logging.debug(f'[{i}] len={len(comp)} ...')
        sg = parent_graph.subgraph(comp)
        c_size = len(sg)
        c_diam = nx.diameter(sg)
        if comps_handled % COMP_INTERVAL == 0:
            logging.debug(f'[{i}] Calculating degree...')

        degree_sequence = sorted(( (d,n) for n, d in sg.degree()), reverse=True)
        (c_max_deg, c_max_deg_seq )  = degree_sequence[0]

        if comps_handled % COMP_INTERVAL == 0:
            logging.debug(f'[{i}] Calculating max node count ...')
        max_count_node = max( sg.nodes(data=True), key=lambda node: node[1]['count'])
        (c_max_count_seq, ndata) = max_count_node
        c_max_count = ndata['count']
        clist = [ i, c_size, c_diam, c_max_deg, c_max_deg_seq, c_max_count, c_max_count_seq ]
        comp_info_list.append(clist)
        intermediate_list.append(clist)
        comps_handled += 1
        if comps_handled % COMP_INTERVAL == 0:
            logging.debug(f'[{i}] Done. Handled {comps_handled} components...') 
        elements_handled += len(comp)
        if outdir is not None:
            if elements_handled > elements_floor:
                of = os.path.join(outdir, f'{intermediates_handled}.compinfo.partial.tsv')
                logging.debug(f'writing output to {of} ...')
                icidf = pd.DataFrame(intermediate_list, columns=COMPINFO_COLUMNS)
                icidf.sort_values(by='size', inplace=True, ascending=False)
                icidf.reset_index(inplace=True, drop=True)                
                icidf.to_csv(of, sep='\t')
                intermediate_list = []
                intermediates_handled += 1
                elements_floor = elements_handled + WRITE_INTERVAL
                
    if outdir is not None:
        intermediates_handled += 1
        of = os.path.join(outdir, f'{intermediates_handled}.compinfo.partial.tsv')
        logging.debug(f'writing output to {of} ...')
        icidf = pd.DataFrame(intermediate_list, columns=COMPINFO_COLUMNS)
        icidf.sort_values(by='size', inplace=True, ascending=False)
        icidf.reset_index(inplace=True, drop=True)                
        icidf.to_csv(of, sep='\t')        

    logging.info(f'handled {comps_handled} total components.')
    cidf = pd.DataFrame(comp_info_list, columns=COMPINFO_COLUMNS)
    cidf.sort_values(by='size', inplace=True, ascending=False)
    cidf.reset_index(inplace=True, drop=True)
    
    # write component assessment info. 
    #of = os.path.join( gdir , f'component_info.tsv')
    #cidf.to_csv( of, sep='\t')
    logging.debug(f'all components len={len(components)}')
    
    return cidf
    

def display_diff(a, b):
    '''
    Mark edits between two sequences, 
    
    CGACATCTACCCAGCGCCCTTTTTGGCCTT
    ||||||||-|||||-|||||||||||||-|
    CGACATCTCCCCAGTGCCCTTTTTGGCCCT
    '''
    dstr = ''
    for i, c in enumerate(a):
        if c != b[i]:
            dstr += '-'
        else:
            dstr += '|'
    out = f'{a}\n{dstr}\n{b}'
    return out


def parse_edges(edgefile):
    '''
    parse edges.txt file. 
    
    ['20163', '19273']
    ['20163', '14372']
    
    '''
    edge_lists = []
    
    with open(edgefile) as f:
        lines = f.readlines()
        logging.debug(f'got {len(lines)} edge lines')
        for line in lines:
            edgepair = eval(line.strip())
            edge_lists.append(edgepair)            
    logging.debug(f'got list of {len(edge_lists)} edges')
    return edge_lists

def make_edge_df(edgelist, unique_df, column='sequence'):
    '''
    make edge df containing 
    pairs of actual sequences 

    edgelist may have string index numbers

    '''
    edgelist_list = []
    
    for pairlist in edgelist:
        a = int( pairlist[0])
        b = int( pairlist[1])
        aseq = unique_df.iloc[a][column]
        bseq = unique_df.iloc[b][column]
        ham = calc_hamming(aseq, bseq)
        edgelist_list.append([aseq, bseq, ham])
    
    df = pd.DataFrame(edgelist_list, columns=['aseq', 'bseq', 'hamming'])
    return df 

 
def parse_components(compfile):
    '''
    parse components.txt file. 
    
    '''
    component_lists = []
    
    with open(compfile) as f:
        lines = f.readlines()
        logging.debug(f'got {len(lines)} component lines')
        for line in lines:
            complist = eval(line.strip())
            component_lists.append(complist)            
    logging.debug(f'got list of {len(component_lists)} components')
    return component_lists


def assess_topx_components( component_df, unique_df, component_lists, top_x = 10 ):
    '''
    @arg component_df        Product of make_component_df()
    @arg component_lists     List of Lists representation of components.txt 
    @arg unique_df           Pandas Dataframe with unique vbc_read and count columns. 
   
    
    '''
    for i in range(0,top_x):
        logging.debug(f'handling index {i}')
        comp_idx = comp_df.iloc[i]['cmp_idx'].astype(int)
        comp = component_lists[comp_idx]    
        slist = []
        for idx in comp:
            s = unique_df.iloc[idx].vbc_read
            slist.append(s)
        g, deg_df = assess_component(slist)


def assess_component( seq_list, max_mismatch=3 ):
    '''
    Make graph and node degree dataframe for one component. 
    
    '''
    g = make_nxgraph_seqlist(seq_list, max_mismatch )
    deg_df = make_degree_df(g)
    return g, deg_df


def check_components(uniques_file, 
                     components_file,
                     edges_file, 
                     column='vbc_read', cp=None):
    '''
    Called  via command line out of band of pipeline.    
    Handles filenames, params, and calls functions. 
    
    '''
    udf = load_mapseq_df( uniques_file, use_dask=False)
    logging.debug(f'loaded. len={len(udf)} dtypes =\n{udf.dtypes}')
    edges = parse_edges(edges_file)
    
    
     
    components = parse_components(components_file)   
    comp_df = make_component_df(udf, components, column=column )
    return comp_df

     
def make_component_df(uniques_df, 
                      component_lists, 
                      min_seq_count = None, 
                      top_x = None, 
                      column='vbc_read'):
    '''

    Take top_x / all components and calculate fraction that exceed 3 max_hamming. 
    
    @arg uniques_df          Pandas Dataframe with unique vbc_read and count columns. 
    @arg component_lists     List of Lists representation of components.txt 
    @arg min_seq_count       Only process components with more sequences.
    
    @return comp_df          Dataframe with component info. 
                                 cmp_idx      index within uniques DF
                                 n_seq        Number of sequences in component.
                                 max_ham      Maximum pairwise Hamming distance in component.
                                 n_pairs      Number of unique pairs
                                 n_exceed     Number that exceed 3
                                 good_prop    Proportion less than 3
    '''
    logging.debug(f'assessing uniques/components')
    start = dt.datetime.now()
   
    udf = uniques_df
    comps = component_lists
    logging.debug(f' {len(comps)} components ')
    logging.debug(f' {len(udf)} unique sequences ')
    
    sizemap = []
    for i, clist in enumerate(comps):
        csize = len(clist)
        sizemap.append( [ i, csize ])
    comp_df = pd.DataFrame(sizemap, columns=['cmp_idx','n_seq'])
    comp_df.sort_values('n_seq', ascending=False, inplace=True)
    comp_df.reset_index(inplace=True, drop=True)   
    logging.debug(f'stats: \n{comp_df.n_seq.describe()}')

    if min_seq_count is not None:
        comp_df = comp_df[comp_df['n_seq'] >= min_seq_count ]

    if top_x is None:
        top_x = len(comp_df)
    
    max_hamming_list = []
    n_pairs_list = []
    n_exceed_list = []
    good_prop_list = []
    
    for i in range(0, top_x):
        comp_idx = comp_df.iloc[i].cmp_idx
        comp = comps[comp_idx]    
        slist = []
        for idx in comp:
            s = udf.iloc[idx][column]
            slist.append(s)
        logging.debug(f'made list of len={len(slist)} component sequence list. checking...')
        ( hmax, n_pairs, n_exceed, max_ok)=  max_hamming(slist)
        logging.debug(f'got properties. hmax={hmax} n_pairs={n_pairs} n_exceed={n_exceed}')
        try:
            good_prop = 1.0 - (n_exceed / n_pairs)
        except ZeroDivisionError:
            good_prop = 0.0

        max_hamming_list.append(hmax)
        n_pairs_list.append(n_pairs)
        n_exceed_list.append(n_exceed)
        good_prop_list.append( good_prop)
    
    logging.debug(f'adding columns. max_hamming...')
    comp_df['max_ham'] = pd.Series(max_hamming_list)        
    logging.debug(f'adding columns. n_pairs...')
    comp_df['n_pairs'] = pd.Series(n_pairs_list, dtype='int')
    logging.debug(f'adding columns. n_exceed...')
    comp_df['n_exceed'] =  pd.Series(n_exceed_list, dtype='int')
    logging.debug(f'adding columns. good_prop...')
    comp_df['good_prop'] = pd.Series(good_prop_list, dtype='float')
    logging.debug(f'made component properties DF...')

    comp_df.fillna(0, inplace=True)
    for cn in ['max_ham','n_pairs','n_exceed']:
        comp_df[cn] = comp_df[cn].astype(int)   
    end = dt.datetime.now()
    td = end - start
    logging.info(f'handled {top_x} components in {td.total_seconds()} seconds.')
    return comp_df

#####################################################
#
#            NetworkX based processing
#
#####################################################


def make_degree_df( nxgraph):
    '''
    Make sorted list of node degree by sequence. 
    
    '''            
    dlist = []
    for n in nxgraph.nodes:
        dlist.append([n, nxgraph.degree[n]] )
    deg_df = pd.DataFrame(dlist, columns=['sequence','node_degree'])
    deg_df.sort_values('node_degree', ascending=False, inplace=True) 
    deg_df.reset_index(inplace=True, drop=True)
    
    return deg_df


def make_nxgraph_bt2df(btdf, 
                       seqdf, 
                       seq_col='vbc_read',
                       count_col = 'count', 
                       seq_len = 30, 
                       min_score = -18.6):
    '''
    Make networkx graph with node properties (read_count) and edge properties (positive 
    bowtie2 score) for more sophisticated algorithms with weighting. 
    
    seqdf has sequences corresponding to name_read and name_align values in btdf.  
    columns  vbc_read, count
    
    Bowtie 2 BTDF columns:
    'name_read', 'flagsum', 'name_align', 'offset', 'qual', 'cigar', 'mate',
       'mate_offset', 'fraglen', 'seq', 'quals', 'score', 'next', 'n_amb',
       'n_mismatch', 'n_gap', 'n_gapext', 'distance', 'md', 'yt']
    
    '''
    
    # Logging intervals...
    SEQ_INTERVAL = 10000
    EDGE_INTERVAL = 100000
    COMP_INTERVAL = 1000 
    
    # Edge info
    btdf['pscore'] = btdf['score'] + abs(min_score)
    
    # Node info
    seq_list = list( seqdf[seq_col] )
    count_list = list( seqdf[count_col])      
    
    # these are indices, so must be integers
    read_list = list( btdf['name_read'].astype(int))
    align_list = list( btdf['name_align'].astype(int))
    pscore_list = list( btdf['pscore'].astype(float))
    distance_list = list( btdf['distance'].astype(int))
    
    #G = nx.DiGraph()
    G = nx.Graph()
    edges_handled = 0
    for i in range(0, len(btdf)):
        rseq = seq_list[ read_list[i] ]
        aseq = seq_list[ align_list[i] ]
        #prop_dict = { "pscore": pscore_list[i] }
        G.add_edge(rseq, aseq, distance= distance_list[i] )
        #G.add_edge(aseq, rseq, distance= distance_list[i] )        
        edges_handled += 1
        if edges_handled % EDGE_INTERVAL == 0:
            logging.debug(f'handled {edges_handled} edges...')
    logging.info(f'handled all {edges_handled} edges.')
    
    # set nodes properties.
    nodes_handled = 0 
    n_singletons = 0
    for i in range(0, len(seq_list)):
        seq = seq_list[i]
        try:
            G.nodes[seq]['count'] = count_list[i]
            G.nodes[seq]['df_idx'] = i
            nodes_handled += 1
            if nodes_handled % SEQ_INTERVAL == 0:
                logging.debug(f'handled {nodes_handled} edges...')
        except KeyError:
            n_singletons += 1
    logging.debug(f'handled {nodes_handled} nodes. n_singletons={n_singletons}')
    logging.info(f'made Graph: nodes={len(G.nodes)} edges={len(G.edges)} ')
    return G


def plot_subgraph( subgraph, comp_idx = '1234'):
    '''
    matplotlib subgraph
    '''
    plt.rcParams["figure.figsize"] = [1.2 * 11.5, 1.2 * 8.5]
    
    max_count_node = max( subgraph.nodes(data=True), key=lambda node: node[1]['count'])

    degree_sequence = sorted(( (d,n) for n, d in subgraph.degree()), reverse=True)    
    (c_max_deg, c_max_deg_seq )  = degree_sequence[0]    
    max_degree_node = ( c_max_deg_seq, subgraph.nodes(data=True)[c_max_deg_seq] )
        
    for u,v, data in subgraph.edges(data=True):
        if 'distance' in data and data['distance'] != 0:
            # distance 3 = pull 1
            # distance 2 = pull 2
            # distance 1 = pull 3
            subgraph[u][v]['pull'] = 4 - data['distance']
    
    
    pos = nx.spring_layout(subgraph, 
                           seed=55, 
                           weight='pull')
    nx.draw_networkx_nodes(subgraph, pos, 
                           node_color='#1f78b4', 
                           node_size=2)
    distances = nx.get_edge_attributes(subgraph, 'distance')
    pulls = nx.get_edge_attributes(subgraph, 'pull')
    #nx.draw_networkx_edges(subgraph, pos, width=[w for w in weights.values()])
    nx.draw_networkx_edges(subgraph, pos, width=1)
    if len(subgraph) < 50:
        nx.draw_networkx_labels(subgraph,  pos, font_size=9, font_color="black")
    nx.draw_networkx_edge_labels(subgraph, pos, font_size=5, edge_labels=distances)
    plt.axis("off")
    plt.show()



def get_components_nx( graph):
    '''    
    complist is list of lists of node sequences
    idxlist  is parallel list o lists of node indexes in unique_df

    assumes "df_idx" data value in every graph node. 

    '''
    scc = nx.connected_components(graph)
    # scc is a generator, process...
    compsets = list( scc )
    # each element is a set, convert to lists
    complist = [  list(x) for x in compsets ]    
    idxlist = []    
    graphdict = graph.nodes(data=True)
    
    for comp in complist:
        nidxlist = []
        for node in comp:
            nidxlist.append( graphdict[node]['df_idx'] )
        idxlist.append(nidxlist)
        
    return ( complist, idxlist )


def get_subgraph( parentgraph, nodelist):
    '''
    
    '''    
    sg = G.subgraph( nodelist)
    return sg 


def check_components_nx( uniques_file, 
                         components_file,
                         edges_file, 
                         column='vbc_read', cp=None):
    '''
    Called  via command line out of band of pipeline.    
    Handles filenames, params, and calls functions. 
    
    '''
    udf = load_mapseq_df( uniques_file, use_dask=False)
    logging.debug(f'loaded. len={len(udf)} dtypes =\n{udf.dtypes}')
    edges = parse_edges(edges_file)
    
    
     
    components = parse_components(components_file)   
    comp_df = make_component_df(udf, components, column=column )
    return comp_df





def align_collapse_nx_grouped(df,
                              column='vbc_read',
                              pcolumn='read_count',
                              gcolumn='brain', 
                              max_distance = None,
                              max_recursion = None, 
                              outdir= None, 
                              datestr= None,
                              force= False,
                              min_reads = 1, 
                              cp=None):
    '''
    Groups alignment and collapse by gcolumn value. [brain]
    Uses sub-directories for standard intermediate output/scratch. 
    Assumes dataframe with sequence (vbc_read) and read_count columns
    
    column     which column to collapse on. 
    gcolumn    which column to group by 
    pcolumn    drives replacement policy. 
                count       number of appearances of unique vbc_read
                read_count  sequencing read_count behind unique 
                degree      node degree in graph of vbc_read.
    min_reads  threshold read_count by this value before collapse. 
    
    '''
    if cp is None:
        cp = get_default_config()
    
    if max_distance is None:
        max_distance = int(cp.get('collapse', 'max_distance'))
    else:
        max_distance = int(max_distance)

    if min_reads is None:
        min_reads = int(cp.get('collapse', 'min_reads', fallback=1))
    else:
        min_reads = int(min_reads)

    if max_recursion is not None:
        rlimit = int(max_recursion)
        logging.info(f'set new recursionlimit={rlimit}')
        sys.setrecursionlimit(rlimit)
    else:
        rlimit = int(cp.get('collapse','max_recursion'))
        logging.info(f'(from config) set new recursionlimit={rlimit}')
        sys.setrecursionlimit(rlimit)        

    aligner = cp.get('collapse','tool')

    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    
    if outdir is None:
        outdir = './'
    outdir = os.path.abspath(outdir)    
    os.makedirs(outdir, exist_ok=True)

    logging.debug(f'collapse: aligner={aligner} max_distance={max_distance} outdir={outdir}')    
    sh = get_default_stats()
    
    sh.add_value('/collapse','n_initial_sequences', len(df) )
    if min_reads > 1:
        logging.info(f'thresholding by min_reads={min_reads} before collapse. len={len(df)} ')
        df = df[ df['read_count'] >= min_reads  ]
        df.reset_index(inplace=True, drop=True)
        logging.info(f'after min_reads={min_reads} filter. len={len(df)} ')
    sh.add_value('/collapse','n_thresholded_sequences', len(df) )    
    sh.add_value('/collapse','n_full_sequences', len(df) )
    sh.add_value(f'/collapse','api_max_distance', str(max_distance) )

    gdflist = []
    gidlist = None
    nogroup_df = None

    if gcolumn is None:
        logging.info(f'gcolumn = None. Collapse all globally.')
        gidlist = [ None ]
    else:
        df[gcolumn] = df[gcolumn].astype('string')
        gidlist = list( df[gcolumn].dropna().unique() )
        gidlist = [ x for x in gidlist if len(x) > 0 ]
        gidlist.sort()
        logging.debug(f'gcolumn={gcolumn} group list: {gidlist}')

        logging.info(f'pulling out no group for merging at end...')
        nogroup_df = df[ df[gcolumn] == '' ]
        sh.add_value('/collapse','n_no_group', len(nogroup_df) )
    
    for gid in gidlist:
        if gid is None:
            gdf = df.copy()
            gdf.reset_index(inplace=True, drop=True)
            initial_len = len(gdf)
            gcolumn = 'all'
            gid = 'all'                          
        else:
            logging.info(f"collapsing '{column}' by {gcolumn} = '{gid}' ")

            gdf = df[df[gcolumn] == gid]
            gdf.reset_index(inplace=True, drop=True)
            initial_len = len(gdf)  
        
        gdir = os.path.join( outdir, f'{gcolumn}.{gid}' )
        os.makedirs(gdir, exist_ok=True)
        logging.info(f'[{gcolumn}:{gid}] initial len={len(gdf)} subdir={gdir}')        
        
        # get reduced dataframe of unique head sequences
        #logging.info('Getting unique DF...')    
        logging.debug('Getting unique DF with sum of counts...')
        vcdf = gdf[column].value_counts().reset_index()
        cdf = pd.merge( vcdf, gdf, on=column, how='left')
        udf = cdf.groupby(column).agg( {'count':'first', pcolumn:'sum'}).reset_index()
        udf.sort_values('count', ascending=False).reset_index(inplace=True, drop=True)
        logging.debug('Unique DF created...')
        
        #udf = gdf[column].value_counts().reset_index() 
        sh.add_value(f'/collapse/{gcolumn}_{gid}','n_unique_sequences', len(udf) )    
            
        # handle writing unique fasta    
        of = os.path.join( gdir , f'{column}.unique.fasta')      
        if os.path.exists(of) and not force:
            logging.debug(f'Output {of} exists and not force.')
            seqfasta = of
        else: 
            logging.info(f'Writing uniques as FASTA to {of}')
            seqfasta = write_fasta_from_df(udf, outfile=of, sequence=[column])        
     
        # run allXall bowtiex
        of = os.path.join( gdir , f'unique_sequences.bt2.sam')
        if os.path.exists(of) and not force:
            logging.debug(f'Output {of} exists and not force.')
            btfile = of
        else:
            logging.info(f'Running {aligner} on {seqfasta} file to {of}')
            btfile = run_bowtie(cp, seqfasta, of, tool=aligner) 
            
        logging.info(f'Bowtie done. Produced output {btfile}. Checking DF.')
        of = os.path.join( gdir , f'unique_sequences.btdf.tsv')
        if os.path.exists(of) and not force:
            logging.debug(f'Output {of} exists and not force.')
            btdf = load_bowtie_df(of)
        else:
            logging.info(f'Creating btdf dataframe from {btfile} max_distance={max_distance}')    
            btdf = make_bowtie_df(btfile, 
                                  max_distance=max_distance, 
                                  ignore_self=True)
            logging.debug(f'Writing output to {of}')
            btdf.to_csv(of, sep='\t')
        sh.add_value(f'/collapse/{gcolumn}_{gid}','n_bowtie_entries', len(btdf) )

        G = make_nxgraph_bt2df(btdf, udf, seq_col=column, count_col='count' )
        
        # Write edges
        of = os.path.join( gdir , f'edgelist.txt')
        edgelist = list(G.edges)
        writelist(of, edgelist)
        logging.debug(f'edgelist len={len(edgelist)}')
        sh.add_value(f'/collapse/{gcolumn}_{gid}','n_edges', len(edgelist) )

        # Set node degree in unique_df. 
        # fulldf[f'{column}_col'] = fulldf[column].map(smd, na_action='ignore')
        udf['degree'] = udf[column].map( dict( G.degree ) , na_action='ignore')
        # Missing values are not in graph, therefor singletons, degree=0
        udf['degree'] = udf['degree'].fillna(0)
        udf['degree'] = udf['degree'].astype(int)

        # Now that we have info, time to write unique df. 
        of = os.path.join( gdir , f'{column}.unique.tsv')
        logging.info(f'Writing unique DF to {of}')
        udf.to_csv(of, sep='\t')


        logging.info(f'finding components.')
        complist, idxlist = get_components_nx(G) 
        logging.info(f'done finding {len(complist)} components')
        
        # write components (sequences)
        of = os.path.join( gdir , f'components.txt')
        writelist(of, complist)
        logging.debug(f'all components len={len(complist)}')

        # write components (indexes)
        of = os.path.join( gdir , f'comp_indexes.txt')
        writelist(of, idxlist)
        logging.debug(f'all comp indexes len={len(idxlist)}')

        # Gather component information. 
        cidf = make_component_df_nx(complist, G, outdir=gdir)
        logging.debug(f'got component info DF len={len(cidf)}')

        # write component assessment info. 
        of = os.path.join( gdir , f'{gid}.component_info.tsv')
        cidf.to_csv( of, sep='\t')
        
        # collapse sequences to max_degree_seq
        newdf = collapse_by_components_nx(gdf, 
                                          udf, 
                                          idxlist, 
                                          cidf,
                                          column=column, 
                                          pcolumn=pcolumn, 
                                          outdir=gdir)
        gdflist.append(newdf)

    # merge all brains into one dataframe...
    logging.debug(f'sizes={[ len(x) for x in gdflist ]}' )
    if nogroup_df is not None:
        logging.debug(f'adding nogroup len={len(nogroup_df)}')
        gdflist.append(nogroup_df)
    outdf = pd.concat(gdflist, ignore_index = True)
    outdf.reset_index(inplace=True, drop=True)
    logging.info(f'All groups. Final DF len={len(outdf)}')
    return outdf


def collapse_by_components_nx(fulldf, 
                              uniqdf, 
                              comp_indexes,
                              component_info_df, 
                              column, 
                              pcolumn='read_count', 
                              outdir=None):
    '''
    comp info columns:
      cmp_idx        id in components.txt
      size           # sequences in component. 
      diam           # max pairwise edit distance.
      max_deg        # max node degree in component. 
      max_deg_seq    # seq with max node degree
      max_count      # max count in component. 
      max_count_seq  # seq with max count
      
    newdf = collapse_by_components_pd(gdf, udf, mcomponents, column=column, pcolumn=pcolumn, outdir=gdir)
    networkx version of applying collapse info. 

    '''
    logging.info('building sequence map dict...')
    smd = build_seqmapdict_nx(uniqdf, comp_indexes, column, pcolumn=pcolumn)
    logging.info('seqmapdict built. Applying.')
    if outdir is not None:
        outfile = f'{outdir}/seqmapdict.json'
        logging.debug(f'writing seqmapdict len={len(smd)} tp {outfile}')
        with open(outfile, 'w') as fp:
            json.dump(smd, fp, indent=4)
    else:
        logging.debug(f'no outdir given.')
    logging.info(f'applying seqmapdict...')
    # make deep copy of original sequence column
    fulldf.loc[:, f'{column}_col'] = fulldf.loc[:, column]    
    logging.info(f'mapping old {column} values to new {column}_col')
    fulldf[f'{column}_col'] = fulldf[column].map( smd, na_action='ignore')
    fulldf = fulldf.fillna( { 'vbc_read_col' : fulldf[column] })    
    logging.info(f'New collapsed df = \n{fulldf}')
    log_objectinfo(fulldf, 'fulldf')
    return fulldf


def build_seqmapdict_nx(udf, 
                        comp_idxlist, 
                        column='vbc_read', 
                        pcolumn='read_count'):
    '''
    May not be right approach for NX data. 
    Create mappings from all unique sequences to component sequence
    Can we do this faster? 
    dict should be oldsequence -> newsequence
    ???
    
    '''
    seqmapdict = {}
    comphandled = 0
    pcolumn='count'
    logging.debug(f'udf len={len(udf)} components len={len(comp_idxlist)} column={column} pcolumn={pcolumn} ')
    comphandled_interval = 1000
    comp_len = len(comp_idxlist)    
    for i, indexlist in enumerate( comp_idxlist):
        cdf = udf[[column, pcolumn]].iloc[indexlist]
        cdf.reset_index(inplace=True, drop=True)
        #logging.debug(f'component [{i}/{comp_len}]: len={len(cdf)}')
        maxid = int(cdf[pcolumn].idxmax())
        t = cdf[column].iloc[maxid]
        for compseq in list(cdf[column]): 
            #logging.debug(f'compseq={compseq} -> {t}')
            seqmapdict[compseq]= t    
        if comphandled % comphandled_interval == 0:
            logging.debug(f'[{i}/{comp_len}]: len={len(cdf)} seq = {t} ')
        comphandled += 1
    return seqmapdict


if __name__ == '__main__':
    
    graph = {
        0 : [1],
        1 : [2],
        2 : [1,3],
        3 : [3],
    }
    print( graph )



