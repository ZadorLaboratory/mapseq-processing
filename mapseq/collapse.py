import itertools
import logging
import os
import subprocess
import sys
import traceback

import datetime as dt

from configparser import ConfigParser
from collections import defaultdict

import pandas as pd
import numpy as np

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
    # newdf has sequence and newsequence columns, rename to orig_seq and sequence
    newcol = f'{column}_col'        
    newdf.rename( {'newsequence': newcol}, inplace=True, axis=1)    
    logging.info(f'Got collapsed DF. len={len(newdf)}')

    rdf = newdf[[ column, newcol, 'read_count' ]]
    of = os.path.join( outdir , f'read_collapsed.tsv')
    logging.info(f'Writing reduced mapping TSV to {of}')
    rdf.to_csv(of, sep='\t')
    
    newdf.drop(column, inplace=True, axis=1)
    newdf.rename( { newcol : column }, inplace=True, axis=1)       
    return newdf        


def align_collapse_pd_grouped(df,
                              column='vbc_read',
                              pcolumn='read_count',
                              gcolumn='brain', 
                              max_mismatch=None,
                              max_recursion=None, 
                              outdir=None, 
                              datestr=None,
                              force=False, 
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
            btdf = make_bowtie_df(btfile, max_mismatch=max_mismatch, ignore_self=True)
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
    
        rdf = newdf[[ column, newcol, 'read_count' ]]
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


def collapse_by_components_pd(fulldf, uniqdf, components, column, pcolumn, outdir=None):
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
    smd = build_seqmapdict_pd(uniqdf, components,  column, pcolumn)
      
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
    
    fulldf.fillna( { 'vbc_read_col' : fulldf['vbc_read'] }, inplace=True)    
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
    Calculate max hamming for set of sequences. 
    Keep track of how many exceed a given max. 
    
    '''
    hmax = 0
    n_pairs = 0
    n_exceed = 0
    for a,b in itertools.combinations(seq_list, 2):
        n_pairs += 1
        hd = calc_hamming(a,b)
        if hd > max_ok:
            n_exceed += 1
        if hd > hmax:
            hmax = hd
    logging.debug(f'considered {n_pairs} pairs, {n_exceed} were > {max_ok}')
    return ( hmax, n_pairs, n_exceed, max_ok)       
    
    
def calc_hamming( a, b):
    '''
    simple calc. 
    '''
    mismatch = 0
    for i,c in enumerate(a):
        if c != b[i]:
            mismatch += 1
    return mismatch


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






if __name__ == '__main__':
    
    graph = {
        0 : [1],
        1 : [2],
        2 : [1,3],
        3 : [3],
    }
    print( robust_topological_sort(graph) )



