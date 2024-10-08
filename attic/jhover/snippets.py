

#
#  Biopython based code
#  Safe, but less efficient w/ memory. 
#


def trim_fasta(config, infile, outdir=None, length=44):
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    if outdir is not None:
        dirname = os.path.abspath(outdir)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)
    head = filename.split('.')[0]    
    logging.debug(f'handling {filepath}')
    
    ofpath = f'{dirname}/{head}.{length}.fasta'
    logging.debug(f'opening {ofpath}...')
    outfile = open(ofpath, 'w')    
    trimmed = []
    sfa = SeqIO.parse(filepath, "fasta")
    for sr in sfa:
        tseq = sr.seq[:length]
        tsr = SeqRecord( tseq, id=sr.id, name=sr.name, description=sr.description)
        trimmed.append(tsr)
    SeqIO.write(trimmed, outfile, 'fasta')
    logging.debug(f'wrote {len(trimmed)} to {ofpath}')
    return ofpath



def calc_min_target(config, braindf):
    '''
    how many molecules (unique UMIs) are in supposedly target-negative area?
    
    '''
    countlist = []
    min_target = 0
    braindf['umi_count'] = braindf['umi_count'].astype(int)
    tndf = braindf[ braindf['site'] == 'target-negative']
    tndf = tndf[ tndf['type'] == 'real']
    lablist = list(tndf['label'].dropna().unique())
    for label in lablist:
        ldf = tndf[tndf['label'] == label]
        if len(ldf) > 0:
            countlist.append( ldf['umi_count'].sum())
    if len(countlist) > 0:
        min_target = max(countlist)
        logging.debug(f'calculated min_target={min_target}')
    return min_target


def max_hamming(sequence, sequencelist):
    '''
    calculates maximum mismatch between sequence and all sequences in sequencelist. 
    assumes all sequences are same length
    no indels, just substitutions. 
    '''
    #logging.debug(f'seq={sequence}')
    #logging.debug(f'sequencelist={sequencelist}')
    max_dist = 0
    for s in sequencelist:
        dist = 0
        for i in range(0,len(s)):
            if sequence[i] != s[i]:
                dist += 1
        if dist > max_dist:
            max_dist = dist
    return max_dist









#
# previously used to calculate read thresholds. 
#

def calc_final_thresholds(config, threshdf):
    '''
    take threshold df for all sites, and derive final thresholds df for
    
    threshdf columns used:   site  count_threshold   
    
    target_threshold = 100
    target-control_threshold = 1000
    target-negative_threshold = 100
    target-lone_threshold = 100
    injection_threshold = 2
    injection-control_threshold=2
    
    
        'site'  'threshold'
    
    
    '''
    tdf = pd.concat( [threshdf[threshdf['site'] == 'target-negative'], 
                      threshdf[threshdf['site'] == 'target']] )
    idf = threshdf[threshdf['site'] == 'injection']
    
    target_thresh = int(tdf['count_threshold'].min())
    inj_thresh = int(idf['count_threshold'].min())

    finaldf = pd.DataFrame( data=[ ['target', target_thresh ],['injection', inj_thresh] ], 
                            columns= ['site','threshold'] )
                  
    return finaldf
    

def calc_thresholds_all(config, sampdf, filelist, fraction=None ):
    '''
    reads in all counts.df (assumes counts column).
     
    calculates thresholds for 'target' and 'injection'
    
    returns 2 dfs. one general info, one with final thresholds
    '''
    if fraction is not None:
        config.set('ssifasta','count_threshold_fraction', fraction)
    
    outlist = []
    
    for filename in filelist:
       logging.debug(f'handling {filename}') 
       (rtprimer, site, brain, region) = guess_site(filename, sampdf)
       cdf = pd.read_csv(filename ,sep='\t', index_col=0)
       (count_threshold, label, clength, counts_max, counts_min)   = calculate_threshold(config, cdf, site )
       outlist.append( [rtprimer, site, count_threshold, label, clength, counts_max, counts_min    ])
    threshdf = pd.DataFrame(data=outlist, columns=['rtprimer', 'site', 'count_threshold', 'label', 'counts_length', 'counts_max', 'counts_min'  ])
    finaldf = calc_final_thresholds(config, threshdf)   
    
    return (finaldf, threshdf)
     
    

def calculate_threshold(config, cdf, site=None):
    '''
    takes counts dataframe (with 'counts' column) 
    if 'label', use that. 
    and calculates 'shoulder' threshold
    site = ['control','injection','target']   
        Will use relevant threshold. If None, will use default threshold
    
    target_threshold=100
    target_ctrl_threshold=1000
    inj_threshold=2
    inj_ctrl_threshold=2
    
    '''
    count_pct = float(config.get('ssifasta','count_threshold_fraction'))
    min_threshold = int(config.get('ssifasta','count_threshold_min'))
    label = 'BCXXX'
    
    try:
        label = cdf['label'].unique()[0]
    except:
        logging.warn(f'no SSI label in DF')
        
    # assess distribution.
    counts = cdf['read_count']
    clength = len(counts)
    counts_max = counts.max()
    counts_min = counts.min()
    counts_mean = counts.mean()
    logging.info(f'handling {label} length={clength} max={counts_max} min={counts_min} ')
    
    val =  cumulative_fract_idx(counts, count_pct)
    if val < min_threshold:
        logging.warning(f'calc threshold < min...')
    else:
        logging.debug(f'calculated count threshold={val} for SSI={label}')
    count_threshold=max(val, min_threshold)
    
    
    #if site is None:
    #    count_threshold = int(config.get('ssifasta', 'default_threshold'))
    #else:
    #    count_threshold = int(config.get('ssifasta', f'{site}_threshold'))
    #logging.debug(f'count threshold for {site} = {count_threshold}')
    return (count_threshold, label, clength, counts_max, counts_min)



def cumulative_fract_idx_naive(ser, fract):
    '''
    value at index of row that marks cumulative fraction of total. 
    assumes series sorted in descending order. 
    starts with largest value. 
    '''
    sum_total = ser.sum()
    fraction_int = int(fract * sum_total)
    cum_total = 0
    val = 0
    idx = 0
    for idx in range(0, len(ser)):
        val = ser[idx]
        cum_total = cum_total + val
        fraction =  cum_total / sum_total
        if cum_total > fraction_int:
            break
        else:
            logging.debug(f'at idx={idx} fraction is {fraction}')
    logging.debug(f'val={val} idx={idx} cum_total={cum_total} ')
    return val    
        

def cumulative_fract_idx(ser, fract):
    '''
    value at index of row that marks cumulative fraction of total. 
    assumes series sorted in descending order. 
    starts with largest value. 
    '''
    sum_total = ser.sum()
    fraction_int = int(fract * sum_total)
    cumsum = ser.cumsum()
    ltser = cumsum[cumsum < fraction_int]
    if len(ltser) < 1:
        idx = 0
        val = cumsum
    else:
        idx =  ltser.index[-1]
        val = ser[idx]
    logging.debug(f'val={val} idx={idx} ')
    return val    






def filter_min_target(df, min_target=1):
    '''
    
    '''
    logging.debug(f'before threshold inj df len={len(ridf)}')
    ridf = ridf[ridf.counts >= min_injection]
    ridf.reset_index(inplace=True, drop=True)
    logging.debug(f'before threshold inj df len={len(ridf)}')   
    
    mdf = pd.merge(rtdf, ridf, how='inner', left_on='sequence', right_on='sequence')
    incol = mdf.columns
    outcol = []
    selcol =[]
    for c in incol:
        if not c.endswith('_y'):
            selcol.append(c)
            outcol.append(c.replace('_x',''))
    mdf = mdf[selcol]
    mdf.columns = outcol
    logging.debug(f'created merged/joined DF w/ common sequence items.  df=\n{mdf}')
    return mdf    



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


def counts_freq(matlabfile, logscale = 'log10', logcol = 'counts' ):
    '''
    '''
    df = pd.read_csv('barcodematrix.tsv',sep='\t', header=None)
    rowsum = df.sum(axis=1)
    rowsort = rowsum.sort_values(ascending=False)
    df = pd.DataFrame(data=rowsort)
    df.columns = ['counts']
    df  = df.reset_index(drop=True)
       
    df[f'log_{logcol}'] = np.log10(df[f'{logcol}'])

    
    #if logscale == 'log2':
    #    df = np.log2(df)
    #elif logscale == 'log10':
    #    df = np.log10(df)
    
    ax = sns.lineplot(data=df, x=df.index, y=df[ f'log_{logcol}' ])
    ax.set_title('HZ120Vamp2 counts frequency')
    ax.set(xlabel='Sequence Rank', ylabel='log10(BC molecule count)')
    
    #   OR 
    # plt.xlabel('x-axis label')
    # plt.ylabel('y-axis label')
    
    plt.savefig('Z120Vamp2.log10.countsfreq.png')
    #plt.show()

def process_ssifasta_withcollapse(config, infile, outdir=None, site=None, datestr=None):
    '''
    by default, outdir will be same dir as infile
    assumes infile fasta has already been trimmed to remove SSI
    
    site = ['target-control','injection-control','target','target-negative',target-lone']   
    Will use relevant threshold. If None, will use default threshold
    
    '''
    aligner = config.get('ssifasta','tool')
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    
    if outdir is None:
        outdir = "."
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logging.debug(f'made outdir={outdir}')

    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    sh = StatsHandler(config, outdir=outdir, datestr=datestr)
    #sh = get_default_stats()
    
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    logging.debug(f'handling {filepath} base={base}')
    
    # make raw fasta TSV of barcode-splitter output for one barcode. 
    # trim to 44 nt since we know last 8 are SSI  
    logging.debug('calc counts...')
    seqdf = make_fasta_df(config, infile)
    of = os.path.join(dirname , f'{base}.read.seq.tsv')
    seqdf.to_csv(of, sep='\t')
    
    # to calculate threshold we need counts calculated. 
    cdf = make_read_counts_df(config, seqdf, label=base)  
    logging.debug(f'initial counts df {len(cdf)} all reads.')
    
    # these are ***READ*** counts
    of = os.path.join(dirname , f'{base}.read.counts.tsv')
    cdf.to_csv(of, sep='\t') 
        
    threshold = get_read_count_threshold(config, cdf, site)
    logging.debug(f'got threshold={threshold} for site {site}')
    tdf = threshold_read_counts(config, cdf, threshold=threshold)
    logging.info(f'at threshold={threshold} {len(tdf)} unique molecules.')
    
    # thresholded raw counts. duplicates are all one UMI, so set counts to 1. 
    # each row (should be) a distinct UMI, so trim to 32. 
    #tdf['counts'] = 1          
    tdf['sequence'] = tdf['sequence'].str[:32]
    tdf['umi_count'] = 1
    
    # this contains duplicate VBCs with *different* UMIs
    of = os.path.join(dirname , f'{base}.umi.seq.tsv')
    tdf.to_csv(of, sep='\t') 
    
    # now we have actual viral barcode df with *unique molecule counts.*
    vbcdf = make_umi_counts_df(config, tdf)
    of = os.path.join(dirname , f'{base}.umi.counts.tsv')
    vbcdf.to_csv(of, sep='\t')     
    
    # split out spike, real, lone, otherwise same as 32.counts.tsv    
    #spikedf, realdf, lonedf = split_spike_real_lone_barcodes(config, vbcdf)
    spikedf, realdf, lonedf, unmatched = split_spike_real_lone_barcodes(config, vbcdf)
    
    # write out this step...
    realdf.to_csv(os.path.join(outdir , f'{base}.real.counts.tsv'), sep='\t')
    lonedf.to_csv(os.path.join(outdir , f'{base}.lone.counts.tsv'), sep='\t')
    spikedf.to_csv(os.path.join(outdir , f'{base}.spike.counts.tsv'), sep='\t')
    unmatched.to_csv(os.path.join(outdir , f'{base}.unmatched.counts.tsv', sep='\t' ))

    # remove homopolymers in real sequences.
    max_homopolymer_run=int(config.get('ssifasta', 'max_homopolymer_run')) 
    realdf = remove_base_repeats(realdf, col='sequence', n=max_homopolymer_run)
     
    # align and collapse all.         
    acrealdf = align_and_collapse(config, realdf, outdir, base, 'real')
    acspikedf = align_and_collapse(config, spikedf, outdir, base, 'spike')
    aclonedf = align_and_collapse(config, lonedf, outdir, base, 'lone')
    acrealdf.to_csv(os.path.join(outdir , f'{base}.real.tsv'), sep='\t')
    acspikedf.to_csv(os.path.join(outdir , f'{base}.spike.tsv'), sep='\t')
    aclonedf.to_csv(os.path.join(outdir , f'{base}.lone.tsv'), sep='\t')

    # add labels for merging...
    acrealdf['type'] = 'real'
    acspikedf['type'] = 'spike'
    aclonedf['type'] = 'lone'

    outdf = merge_dfs([ acrealdf, acspikedf, aclonedf ])
    outdf['label'] = base
    outdf.sort_values(by = ['type', 'umi_count'], ascending = [True, False], inplace=True)
    outdf.reset_index(drop=True, inplace=True)
    return outdf
    
    
def align_and_collapse(config, countsdf, outdir, base, label):
    '''
    countsdf  'sequence' 'read_count' 'umi_count' columns
    outdir    working dir or temp dir. 
    base      leading file name, e.g. barcode label, e.g. 'SSI4'
    label     type of sequence, e.g. real, spike, L1 (lone)
    
    '''
    sh = get_default_stats()
    newdf = None
    logging.debug(f'handling {base} {label}s...')
    aligner = config.get('ssifasta','tool')
    num_sequences = len(countsdf)
    num_reads = countsdf.read_count.sum()
    logging.info(f'{base} {label} {num_sequences} sequences, representing {num_reads} reads.')      
    # sh.add_value(f'/fastq/pair{pairshandled}','reads_total', num_handled)
    sh.add_value(f'/ssifasta/{base}/{label}','num_sequences',len(countsdf))
    sh.add_value(f'/ssifasta/{base}/{label}','num_reads',len(countsdf))
    
    of = os.path.join( outdir , f'{base}.{label}.seq.fasta')
    logging.debug(f'make fasta for {aligner} = {of}') 
    seqfasta = write_fasta_from_df(countsdf, outfile=of)
    of = os.path.join(outdir , f'{base}.{label}.{aligner}')
    logging.debug(f'running {aligner}...')
    try:
        afile = run_bowtie(config, seqfasta, of, tool=aligner )  
        logging.debug(f'handle {aligner} align file: {afile}')
        btdf = make_bowtie_df(afile)
        of = os.path.join(outdir , f'{base}.{label}.btdf.tsv')
        btdf.to_csv(of, sep='\t') 
        edgelist = edges_from_btdf(btdf)
        components = get_components(edgelist)
        logging.debug(f'countdf columns are {countsdf.columns}')
        newdf = collapse_umi_counts_df(countsdf, components)
        logging.debug(f'orig len={num_sequences}, {len(components)} components, collapsed len={len(newdf)}')
        sh.add_value(f'/ssifasta/{base}/{label}','num_collapsed', len(newdf) )

    except NonZeroReturnException:
        logging.warning(f'NonZeroReturn Exception. Probably no {label}s found. ')
        newdf = pd.DataFrame(columns = ['sequence','counts'])
    return newdf



def collapse_umi_counts_df(countsdf, components):
    '''
    takes components consisting of indices
    determines sequence with largest count columns: 'sequence', 'umi_counts'
    collapses all other member components to the sequence of the largest.
    adds their counts to that of that sequence.
    retain columns and values for highest counts row. 
    
    *** need to *sum*  read_counts column as well... 
    
     
    '''
    # list of lists to collect values..
    logging.debug(f'collapsing countsdf len={len(countsdf)} w/ {len(components)} components.')
    lol = []
    colnames = list(countsdf.columns)
    for component in components:
        logging.debug(f'component={component}')
        # make new df of only component sequence rows
        cdf = countsdf.iloc[component].reset_index(drop=True)
        logging.debug(f'cdf=\n{cdf}')
        # which sequence has highest count?
        maxid = cdf['umi_count'].idxmax()
        # extract sequence and count as python list
        row = list(cdf.iloc[maxid])
        # set counts as sum of all collapse sequences. 
        row[1] = cdf['read_count'].sum()
        row[2] = cdf['umi_count'].sum()
        lol.append(row)
    
        if logging.getLogger().level <= logging.DEBUG:
            slist = list(cdf['sequence'])
            logging.debug(f'slist={slist}')
            if len(slist) > 1:
                s = row[0]
                maxdiff = max_hamming(s, slist)
                logging.debug(f'max_hamming = {maxdiff} n_seqs={len(slist)}')
            else:
                 logging.debug(f'skip distance calc, one sequence in component.')
        
    newdf = pd.DataFrame(data=lol, columns=colnames)
    logging.debug(f'original len={len(countsdf)} collapsed len={len(newdf)}')
    newdf.sort_values('umi_count',ascending=False, inplace=True)
    newdf.reset_index(drop=True, inplace=True)
    return newdf

#


def process_fasta_old(config, sampdf, infile, bclist, outdir, force=False, countsplots=False, readtsvs=False, datestr=None):
    '''
    Take paired FASTA file as input rather than raw FASTQ
    Use combined barcode sorter structure to minimize character comparisons. 

    '''
    if outdir is None:
        outdir = "."
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logging.debug(f'made outdir={outdir}')
    
    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    
    output_exists = check_output(bclist)
    cfilename = f'{outdir}/process_fasta.config.txt'
    bc_length=len(bclist[0].barcode)
    unmatched = os.path.abspath(f'{outdir}/unmatched.fasta')
    
    seqhandled_interval = int(config.get('fastq','seqhandled_interval')) 
    matched_interval = int(config.get('fastq','matched_interval'))
    unmatched_interval = int(config.get('fastq','unmatched_interval'))
    
    logging.info(f'performing split. outdir={outdir} output_exists={output_exists} force={force} countsplots={countsplots} readtsvs={readtsvs}')
    
    
    if ( not output_exists ) or force:
        write_config(config, cfilename, timestamp=True, datestring=datestr)
        sh = StatsHandler(config, outdir=outdir, datestr=datestr)        
        # create structure to search for all barcode simultaneously
        labeldict = {}
        for bco in bclist:
            labeldict[bco.barcode] = bco.label
        matchdict, seqdict, unmatched_file = build_bcmatcher(bclist) 

        pairshandled = 0
        num_handled_total = 0
        num_matched_total = 0
        num_unmatched_total = 0

        # handle all sequences in input 
        with open(infile) as f:
            num_handled = 0
            num_matched = 0
            num_unmatched = 0
            
            logging.debug(f'handling file {infile}')
            while True:
                try:
                    line = f.readline()
                    if line.startswith('>'):
                        pass
                    else:
                        if len(line) == 0:
                            raise StopIteration
                        
                        fullread = line.strip()
                        # handle sequence
                        matched = False
                        seq = fullread[-bc_length:]
                        # id, seq, matchdict, fullseq, unmatched
                        matched = do_match(num_handled, seq, matchdict, fullread, unmatched_file)
                        num_handled += 1
                        num_handled_total += 1
    
                        # report progress...                    
                        if not matched:
                            num_unmatched += 1
                            num_unmatched_total += 1                        
                            if num_unmatched % unmatched_interval == 0:
                                logging.debug(f'{num_unmatched} unmatched so far.')
                        else:
                            num_matched += 1
                            num_matched_total += 1
                            if num_matched % matched_interval == 0:
                                logging.debug(f'match {num_matched}: found SSI in {fullread}!')                        
    
                        if num_handled % seqhandled_interval == 0: 
                            logging.info(f'handled {num_handled} matched={num_matched} unmatched={num_unmatched}')
                    
                except StopIteration as e:
                    logging.debug('iteration stopped')    
                    break
            
        logging.debug(f'finished with {infile}')
        sh.add_value('/fasta','reads_handled', num_handled_total )
        sh.add_value('/fasta','reads_unmatched', num_unmatched_total )
        f.close()
        
        matchrate = 0.0
        if num_matched_total > 0: 
            unmatchrate = num_unmatched_total / num_matched_total              
            matchrate = 1.0 - unmatchrate
        logging.info(f'handled {num_handled_total} sequences. {num_matched_total} matched. {num_unmatched_total} unmatched matchrate={matchrate}')
    else:
        logging.warn('All FASTA output exists and force=False. Not recalculating.')
    
    
    if readtsvs:
        filelist = []
        for bch in bclist:
            if os.path.exists(bch.filename):
                filelist.append(bch.filename)
            else:
                logging.warning(f'missing FASTA file!: {of}')
        logging.info(f'Making counts dfs/TSVs for {filelist} in {outdir}')
        make_read_counts_dfs(config, filelist, outdir)

    if countsplots:
        logging.info('Making combined counts plots PDF...')
        countsfilelist = []
        for bch in bclist:
            base = get_mainbase(bch.filename)   
            of = os.path.join(outdir , f'{base}.read.counts.tsv')
            if os.path.exists(bch.filename):
                countsfilelist.append(of)
            else:
                logging.warning(f'missing countsfile needed for plot!:  {of}')
        make_read_countsplot_combined_sns(config, sampdf, countsfilelist, outdir=outdir, expid=None )




#
# Old fastq merging making a lot of MAPseq-specific assumptions. 
# Removed in favor of more generic approach. No SSIs, sampldf...

def process_fastq_pairs(config, sampdf, readfilelist, bclist, outdir, force=False, countsplots=False, readtsvs=False, datestr=None):
    '''
    Do not use BioConda data structures. 
    Use combined barcode sorter structure to minimize character comparisons. 

    '''
    if outdir is None:
        outdir = "."
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logging.debug(f'made outdir={outdir}')
    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")

    if countsplots:
        readtsvs = True
    
    output_exists = check_output(bclist)
    logging.debug(f'output_exists={output_exists} force={force}')

    cfilename = f'{outdir}/process_fastq.config.txt'
    bc_length=len(bclist[0].barcode)
    outfile = os.path.abspath(f'{outdir}/unmatched.fasta')
    pairedfile = os.path.abspath(f'{outdir}/paired.txt')    
    r1s = int(config.get('fastq','r1start'))
    r1e = int(config.get('fastq','r1end'))
    r2s = int(config.get('fastq','r2start'))
    r2e = int(config.get('fastq','r2end'))    

    seqhandled_interval = int(config.get('fastq','seqhandled_interval')) 
    matched_interval = int(config.get('fastq','matched_interval'))
    unmatched_interval = int(config.get('fastq','unmatched_interval'))
    #seqhandled_interval = 1000000 
    #matched_interval = 1000000
    #unmatched_interval = 1000000
    
    if ( not output_exists ) or force:
        write_config(config, cfilename, timestamp=True, datestring=datestr)
        sh = StatsHandler(config, outdir=outdir, datestr=datestr)        
        # create structure to search for all barcode simultaneously
        labeldict = {}
        for bco in bclist:
            labeldict[bco.barcode] = bco.label
        matchdict, seqdict, unmatched_file = build_bcmatcher(bclist) 

        pf = open(pairedfile, 'w')
        pairshandled = 0
        num_handled_total = 0
        num_matched_total = 0
        num_unmatched_total = 0

        # handle all pairs of readfiles from readfilelist
        for (read1file, read2file) in readfilelist:
            pairshandled += 1
            num_handled = 0
            num_matched = 0
            num_unmatched = 0
            
            logging.debug(f'handling file pair {pairshandled}')
            if read1file.endswith('.gz'):
                r1f = gzip.open(read1file, "rt")
            else:
                r1f = open(read1file)
            
            if read2file.endswith('.gz'):
                r2f = gzip.open(read2file, "rt")         
            else:
                r2f = open(read2file)
        
            while True:
                try:
                    meta1 = r1f.readline()
                    if len(meta1) == 0:
                        raise StopIteration
                    seq1 = r1f.readline().strip()
                    sep1 = r1f.readline()
                    qual1 = r1f.readline().strip()

                    meta2 = r2f.readline()
                    if len(meta2) == 0:
                        break
                    seq2 = r2f.readline().strip()
                    sep2 = r2f.readline()
                    qual2 = r2f.readline().strip()                    

                    sub1 = seq1[r1s:r1e]
                    sub2 = seq2[r2s:r2e]
                    fullread = sub1 + sub2
                    pf.write(f'{fullread}\n')
                    
                    matched = False
                    seq = fullread[-bc_length:]
                    # id, seq, matchdict, fullseq, unmatched
                    matched = do_match(num_handled, seq, matchdict, fullread, unmatched_file)
                    
                    num_handled += 1
                    num_handled_total += 1

                    # report progress...                    
                    if not matched:
                        num_unmatched += 1
                        num_unmatched_total += 1                        
                        if num_unmatched % unmatched_interval == 0:
                            logging.debug(f'{num_unmatched} unmatched so far.')
                    else:
                        num_matched += 1
                        num_matched_total += 1
                        if num_matched % matched_interval == 0:
                            logging.debug(f'match {num_matched}: found SSI in {fullread}!')                        

                    if num_handled % seqhandled_interval == 0: 
                        logging.info(f'handled {num_handled} reads from pair {pairshandled}. matched={num_matched} unmatched={num_unmatched}')
                
                except StopIteration as e:
                    logging.debug('iteration stopped')    
                    break
            
            logging.debug(f'finished with pair {pairshandled} saving pair info.')
            sh.add_value(f'/fastq/pair{pairshandled}','reads_total', num_handled)
            sh.add_value(f'/fastq/pair{pairshandled}','reads_matched', num_matched) 
            sh.add_value(f'/fastq/pair{pairshandled}','reads_unmatched', num_unmatched)
            num_unmatched_total += num_unmatched
   
        sh.add_value('/fastq','reads_handled', num_handled_total )
        sh.add_value('/fastq','reads_unmatched', num_unmatched_total )
               
        pf.close()              
        logging.info(f'handled {num_handled_total} sequences. {pairshandled} pairs. {num_matched_total} matched. {num_unmatched_total} unmatched')
    else:
        logging.warn('All FASTA output exists and force=False. Not recalculating.')
    
    if readtsvs:
        filelist = []
        for bch in bclist:
            filelist.append(bch.filename)
        logging.info(f'Making counts dfs/TSVs for {filelist} in {outdir}')
        make_read_counts_dfs(config, filelist, outdir)

    
    if countsplots:
        logging.info('Making combined counts plots PDF...')
        countsfilelist = []
        for bch in bclist:
            base = get_mainbase(bch.filename)   
            of = os.path.join(dirname , f'{base}.read.counts.tsv')
            countsfilelist.append(of)
        make_read_countsplot_combined_sns(config, sampdf, countsfilelist, outdir=outdir, expid=None )


def make_read_countsplots(config, filelist ): 
    '''
    makes individual read counts plots from reads.counts.tsv files. 
    
    '''   
    for bcfile in filelist:
        logging.debug(f'handling {bcfile}')
        filepath = os.path.abspath(bcfile)    
        dirname = os.path.dirname(filepath)   
        filename = os.path.basename(filepath)
        (base, ext) = os.path.splitext(filename)
        base = base.split('.')[0] 
        
        bcdata = pd.read_csv(bcfile, sep='\t')
        plt.figure()
        plt.plot(np.log10(bcdata['Unnamed: 0']), np.log10(bcdata['read_count']))
        plt.title(base)
        plt.xlabel("log10(BC index)")
        plt.ylabel("log10(BC counts)")
        plt.savefig(bcfile.replace('tsv', 'pdf'))    




def counts_axis_plot_sns(ax, bcdata, labels):
    '''
    Creates individual axes for single plot within figure. 
    
    '''
    bcdata['log10index'] = np.log10(bcdata.index)
    bcdata['log10counts'] = np.log10(bcdata['read_count'])
    sns.lineplot(ax=ax, x=bcdata['log10index'], y=bcdata['log10counts'] )
    s = bcdata['read_count'].sum()
    n = len(bcdata)
    t = bcdata['read_count'].max()

    title = f"{bcdata['label'][0]}"      
    ax.set_title(title, fontsize=10)
    ax.text(0.15, 0.2, f"site={labels['site']}\nn={n}\ntop={t}\nsum={s}\nthreshold={labels['threshold']}", fontsize=9) #add text
    

def make_read_countsplot_combined_sns(config, sampdf, filelist, outdir, expid=None ):    
    '''
     makes combined figure with all plots. 
     assumes column 'label' for title. 
     
    '''
    min_ssi_count = int(config.get('analysis','min_ssi_count')) 
    from matplotlib.backends.backend_pdf import PdfPages as pdfpages
    
    datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    
    outfile = f'countsplots.{datestr}.pdf'
    if expid is not None:
        outfile = f'{expid}.{outfile}'
    outfile = os.path.join(outdir, outfile)
    
    # do nine per figure...
    page_dims = (11.7, 8.27)
    with pdfpages(outfile) as pdfpages:
        #fig_n = math.ceil( math.sqrt(len(filelist)) )
        #fig, axes = plt.subplots(nrows=fig_n, ncols=fig_n, figsize=a4_dims,  layout='constrained')
        plots_per_page = 9
        num_figs = float(len(filelist)) / float(plots_per_page)
        if num_figs % 9 == 0:
            num_figs = int(num_figs)
        else:
            num_figs = int(num_figs) + 1
        logging.debug(f'with {plots_per_page} plots/page, need {num_figs} for {len(filelist)} file plots.')
        
        figlist = []
        axlist = []
        for i in range(0,num_figs):
            fig,axes = plt.subplots(nrows=3, ncols=3, figsize=page_dims,  layout='constrained') 
            if expid is not None:
                fig.suptitle(f'{expid} read counts frequency plots.')
            else:
                fig.suptitle(f'Read counts frequency plots')
            figlist.append(fig)
            # numpy.flatirator doesn't handle indexing
            for a in axes.flat:
                axlist.append(a)
        logging.debug(f'created {len(figlist)} figures to go on {num_figs} pages. ')
                  
        #fig.set_xlabel("log10(BC index)")
        #fig.set_ylabel("log10(BC counts)")
        filelist = natsorted(filelist)
        logging.debug(f'handling {len(filelist)} files...')            
        for i, bcfile in enumerate(filelist):
            logging.debug(f'handling {bcfile}')
            bcdata = pd.read_csv(bcfile, sep='\t')
            if len(bcdata) > min_ssi_count:
                (rtprimer_num, site, brain, region ) = guess_site(bcfile, sampdf )           
                count_threshold, label, clength, counts_max, counts_min = calculate_threshold(config, bcdata)
                labels = {'rtprimer':rtprimer_num,
                          'site':site,
                          'brain':brain,
                          'region': region,
                          'threshold' : count_threshold
                          }
                
                ax = axlist[i]
                counts_axis_plot_sns(ax, bcdata, labels=labels)
            else:
                ax = axlist[i]
                # make empty axis?
        for f in figlist:
            pdfpages.savefig(f)
    logging.info(f'saved plot PDF to {outfile}')
    

#  handling seqmandict merging...
#fulldf['newsequence'] = fulldf.apply(apply_setcompseq, axis=1, seqmapdict=seqmapdict)
# make deep copy of original sequence column
fulldf.loc[:, f'{column}_col'] = fulldf.loc[:, column]
#fulldf[f'{column}_col'] = fulldf[column]

# map old to new
fulldf[f'{column}_col'] = fulldf[column].map(smd, na_action='ignore')
# fill in NaN with original values. 
df['vbc_read_col'].fillna(df['vbc_read'], inplace=True)



# replace is *VERY* slow ?
#fulldf[f'{column}_col'] = fulldf[column].replace(smd)

def aggtest():
    '''
    
    '''
    columns = [ 'vbc', 'vbc_read_col', 'umi','label','read_count']
    lol = [     [ 'AAA', 'AAA', 'M',  'BC1', 2],
                [ 'AAB', 'AAA', 'M',  'BC1', 15],        
                [ 'AAC', 'AAA', 'N',  'BC1', 4],
                [ 'BBB', 'BBB', 'O',  'BC1', 7],
                [ 'BBA', 'BBB', 'P',  'BC1', 3],        
                [ 'BBC', 'BBB', 'Q',  'BC1', 2],        
                [ 'BBC', 'BBB', 'R',  'BC2', 8],
                [ 'BBC', 'BBB', 'R',  'BC2', 9],
                [ 'BBC', 'BBB', 'S',  'BC2', 7],
                [ 'CCC', 'CCC', 'O',  'BC1', 14],
                [ 'CCA', 'CCC', 'O',  'BC1', 13],
                [ 'CCC', 'CCC', 'O',  'BC1', 41],
                [ 'CCC', 'CCC', 'P',  'BC1', 4],        
                [ 'CCA', 'CCC', 'P',  'BC1', 6],        
                [ 'CCB', 'CCC', 'Q',  'BC1', 5], 
                [ 'CCD', 'CCC', 'R',  'BC2', 3],
            ]
                # does it matter if there is another column that always corresponds to an existing? 
    columns2 = [ 'vbc', 'vbc_read_col', 'umi','label','site','read_count']
    lol2 = [     [ 'AAA', 'AAA', 'M',  'BC1', 'target', 2],
                [ 'AAB', 'AAA', 'M',  'BC1', 'target', 15],        
                [ 'AAC', 'AAA', 'N',  'BC1', 'target', 4],
                [ 'BBB', 'BBB', 'O',  'BC1', 'target', 7],
                [ 'BBA', 'BBB', 'P',  'BC1', 'target', 3],        
                [ 'BBC', 'BBB', 'Q',  'BC1', 'target', 2],        
                [ 'BBC', 'BBB', 'R',  'BC2', 'injection', 8],
                [ 'BBC', 'BBB', 'R',  'BC2', 'injection', 9],
                [ 'BBC', 'BBB', 'S',  'BC2', 'injection', 7],
                [ 'CCC', 'CCC', 'O',  'BC1', 'target', 14],
                [ 'CCA', 'CCC', 'O',  'BC1', 'target', 13],
                [ 'CCC', 'CCC', 'O',  'BC1', 'target', 41],
                [ 'CCC', 'CCC', 'P',  'BC1', 'target', 4],        
                [ 'CCA', 'CCC', 'P',  'BC1', 'target', 6],        
                [ 'CCB', 'CCC', 'Q',  'BC1', 'target', 5], 
                [ 'CCD', 'CCC', 'R',  'BC2', 'injection', 3],
            ]    
    df = pd.DataFrame(lol, columns=columns)
    #df[['vbc_col','umi','label']].groupby(['label','vbc_col']).nunique().reset_index()
    df.groupby(['label','vbc_read_col']).agg( {'umi' : 'nunique','read_count':'sum'}).reset_index()
    
    df2 =  pd.DataFrame(lol2, columns=columns2)
    df2.groupby(['label','vbc_read_col']).agg( {'umi' : 'nunique','read_count':'sum', 'site':'first'}).reset_index()

    # same output. 

