

#from kneed import KneeLocator

def calc_kneed_idx(x, y , inflect, poly=2, sense=4):
   '''
  assumes convex, then concave, decreasing curve.
    inflect = 'knee'|'elbow'
    
   '''
   if inflect == 'knee':
       kl = KneeLocator(x=x, y=y, S=sense, curve='convex',direction='decreasing',interp_method='polynomial',polynomial_degree=poly)
       val = kl.knee
       logging.debug(f'got value {val} for knee from kneed...')
   elif inflect == 'elbow':
       # not validated!
       kl = KneeLocator(x=x, y=y, S=sense, curve='convex',direction='decreasing',interp_method='polynomial',polynomial_degree=poly)        
       val = kl.elbow
       logging.debug(f'got value {val} for knee from kneed...')
   return val

def calculate_threshold_kneed(config, cdf, site=None, inflect=None ):
    '''
    takes counts dataframe (with 'counts' column) 
    if 'label', use that. 
    and calculates 'knee' or 'elbow' threshold
    site = ['control','injection','target']   
        Will use relevant threshold. If None, will use default threshold
    
    target_threshold=100
    target_ctrl_threshold=1000
    inj_threshold=2
    inj_ctrl_threshold=2
    
    '''
    if inflect is None:
        inflect = config.get('ssifasta','threshold_heuristic')
    min_threshold = int(config.get('ssifasta','count_threshold_min'))
    label = 'BCXXX'
    
    try:
        label = cdf['label'].unique()[0]
    except:
        logging.warn(f'no SSI label in DF')
        
    # assess distribution.
    counts = cdf['counts']
    clength = len(counts)
    counts_max = counts.max()
    counts_min = counts.min()
    counts_mean = counts.mean()
    logging.info(f'handling {label} length={clength} max={counts_max} min={counts_min} ')
    
    val = calc_kneed_idx(cdf.index, cdf.counts, inflect='knee'  )
    if val < min_threshold:
        logging.warning(f'kneed calc threshold < min...')
    else:
        logging.debug(f'calculated count threshold={val} for SSI={label}')
    count_threshold=max(val, min_threshold)
        
    #if site is None:
    #    count_threshold = int(config.get('ssifasta', 'default_threshold'))
    #else:
    #    count_threshold = int(config.get('ssifasta', f'{site}_threshold'))
    #logging.debug(f'count threshold for {site} = {count_threshold}')
    return (count_threshold, label, clength, counts_max, counts_min)



def process_fastq_pairs_parallel(config, readfilelist, bclist, outdir, nthreads, force=False):
    '''
    
    nthreads:    use this number of CPUs. 0 means all. -1 means all but 1. 3 means 3. 
    
    '''
    ncpus, threads = calc_thread_count(nthreads)   
    logging.info(f'using {threads} of {ncpus} CPUs in parallel...')

    prog = os.path.expanduser('~/git/mapseq-processing/scripts/process_fastq_thread.py')
    readfilestr = ""
    for (a,b) in readfilelist:
        readfilestr += f" {a} {b} "
    
    logging.debug(f'readfilestr = {readfilestr}')
    
    # from cshlwork.utils import JobRunner, JobStack, JobSet
    
    jstack = JobStack()
    
    for bco in bclist:
        cmd = [ prog, 
               '-d',
               '-B', bco.label , 
               '-O' , outdir ,
               readfilestr 
               ]
        jstack.addjob(cmd)
    jset = JobSet(max_processes = threads, jobstack = jstack)
    jset.runjobs()
    
    
    filelist = []
    for bco in bclist:
        filelist.append(bco.filename)
    logging.info(f'Making counts df for {filelist} in {outdir}')
    make_counts_dfs(config, filelist, outdir)


def process_fastq_pairs_single(config, readfilelist, bclist, outdir, force=False):

    # if all the output files for bclist exist, don't recalc unless force=True. 
    if outdir is None:
        outdir = "."
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logging.debug(f'made outdir={outdir}')
    output_exists = check_output(bclist)
    logging.debug(f'output_exists={output_exists} force={force}')
    
    # list should have only one...
    bcho = bclist[0]
    
    if ( not output_exists ) or force:
        #outfile = os.path.abspath(f'{outdir}/unmatched.fasta')
        #pairedfile = os.path.abspath(f'{outdir}/paired.txt')
        #umf = open(outfile, 'w')
        #pf = open(pairedfile, 'w')
        r1s = int(config.get('fastq','r1start'))
        r1e = int(config.get('fastq','r1end'))
        r2s = int(config.get('fastq','r2start'))
        r2e = int(config.get('fastq','r2end'))
        
        seqhandled_interval = int(config.get('fastq','seqhandled_interval')) 
        matched_interval = int(config.get('fastq','matched_interval'))
        unmatched_interval = int(config.get('fastq','unmatched_interval'))

        seqshandled = 0
        pairshandled = 0
        unmatched = 0
        didmatch = 0
    
        #
        # handle pairs of readfiles from readfilelist
        #
        for (read1file, read2file) in readfilelist:
            pairshandled += 1
            logging.debug(f'handling file pair {pairshandled}')
            if read1file.endswith('.gz'):
                 read1file = gzip.open(read1file, "rt")
            if read2file.endswith('.gz'):
                 read2file = gzip.open(read2file, "rt")         
                
            recs1 = SeqIO.parse(read1file, "fastq")
            recs2 = SeqIO.parse(read2file, "fastq")
        
            while True:
                try:
                    r1 = next(recs1)
                    r2 = next(recs2)
                    sub1 = r1.seq[r1s:r1e]
                    sub2 = r2.seq[r2s:r2e]
                    fullread = sub1 + sub2
                    #pf.write(f'{fullread}\n')
                    
                    matched = False
                    r = bcho.do_match(seqshandled, fullread)
                    if r:
                        didmatch += 1
                        if didmatch % matched_interval == 0:
                            logging.debug(f'match {didmatch}: found SSI {bcho.label} in {fullread}!')
                        matched = True
                    else:
                        unmatched += 1
                        # when processing single, unmatched number not useful. 
                        #if unmatched % unmatched_interval == 0:
                        #    logging.debug(f'{unmatched} unmatched so far.')
                        #id = str(seqshandled)
                        #sr = SeqRecord( fullread, id=id, name=id, description=id)
                        #SeqIO.write([sr], umf, 'fasta')
                    
                    seqshandled += 1
                    if seqshandled % seqhandled_interval == 0: 
                        logging.debug(f'handled {seqshandled} reads from pair {pairshandled}. matched={didmatch} unmatched={unmatched}')
                
                except StopIteration as e:
                    logging.debug(f'iteration stopped?')
                    logging.warning(traceback.format_exc(None))
                    break
                
        #umf.close()
        #pf.close()
        #for bch in bclist:
        bcho.finalize()    
        # close possible gzip filehandles??
        #max_mismatch = bclist[0].max_mismatch
        logging.info(f'handled {seqshandled} sequences. {pairshandled} pairs. {didmatch} matched. {unmatched} unmatched')
    else:
        logging.warn('all output exists and force=False. Not recalculating.')




def collapse_by_components_faster(fulldf, uniqdf, components):
    #
    # faster than naive, but still slow
    logging.debug(f'all components len={len(components)}')
    components = remove_singletons(components)
    logging.debug(f'multi-element components len={len(components)}')
    logging.debug(f'fulldf length={len(fulldf)} uniqdf length={len(fulldf)} {len(components)} components.')
    
    groups = fulldf.groupby('sequence').groups
    # convert to python dict for speed
    gldict = {}   
    for k in groups.keys():
        gldict[k] = list(groups[k])
    gldict_len = len(gldict)
    logging.debug(f'grouplist by seq len={gldict_len}')
      
    # Make new full df:
    logging.debug('copying fulldf')
    newdf = fulldf.copy()
    
    #comphandled_interval = int(config.get('fasta','comphandled_interval')) 
    comphandled = 0
    comphandled_interval = 100
    for comp in components:
        ilist = []
        cslist = list(uniqdf['sequence'].iloc[comp])
        s = None
        for s in cslist:
            for i in gldict[s]:
                ilist.append(i)
        newdf['sequence'].iloc[ilist] = s
        
        comphandled += 1
        if comphandled % comphandled_interval == 0:
            logging.debug(f'setting all component seqs to {s}')
            logging.info(f'handled {comphandled}/{gldict_len} components.')
    logging.info(f'new collapsed df = \n{newdf}')
    return newdf

            
def collapse_by_components_naive(fulldf, uniqdf, components):
    #
    # initial version, slow iteration 
    #
    logging.debug(f'all components len={len(components)}')
    components = remove_singletons(components)
    logging.debug(f'multi-element components len={len(components)}')
    logging.debug(f'fulldf length={len(fulldf)} uniqdf length={len(fulldf)} {len(components)} components.')
    glist = fulldf.groupby('sequence').groups
    glist_len = len(glist)
    logging.debug(f'grouplist by seq len={glist_len}. e.g. {glist[list(glist.keys())[1]]}')
    
    # Make new full df:
    newdf = fulldf.copy()
    
    #comphandled_interval = int(config.get('fasta','comphandled_interval')) 
    comphandled = 0
    comphandled_interval = 1000
    
    for comp in components:
        compidx =  pd.Index([], dtype='int64')
        max_seq = None 
        n_max = 0
        for uniqidx in comp:
            seq = uniqdf.sequence.iloc[uniqidx]
            #logging.debug(f'uniqidx = {uniqidx} seq={seq}')
            fullidx = glist[seq]
            idxlen = len(fullidx)
            if idxlen > n_max:
                n_max = idxlen
                max_seq = seq
            #logging.debug(f'seq {seq} -> {fullidx}')
            compidx = compidx.union(fullidx)
            
        #logging.debug(f'fullidx len={len(compidx)}: {compidx}')
        #logging.debug(f'setting seq to {max_seq}')
        newdf.sequence.iloc[compidx] = max_seq
        
        comphandled += 1
        if comphandled % comphandled_interval == 0:
            logging.info(f'handled {comphandled}/{glist_len} components.')
    logging.info(f'new collapsed df = \n{newdf}')
    return newdf
