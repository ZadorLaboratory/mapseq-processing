

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




def process_fastq_pairs_pd(infilelist, 
                            outdir,                         
                            force=False, 
                            cp = None):
    '''
    only parse out read lines to pandas, then join with pandas. 
    CPU times: user 9min 14s,
    Output:
        sequence    vbc_read    umi    ssi
    
    '''
    r1s = int(cp.get('fastq','r1start'))
    r1e = int(cp.get('fastq','r1end'))
    r2s = int(cp.get('fastq','r2start'))
    r2e = int(cp.get('fastq','r2end'))
    logging.debug(f'read1[{r1s}:{r1e}] + read2[{r2s}:{r2e}]')
    df = None
    sh = get_default_stats()
    for (read1file, read2file) in infilelist:
        logging.info(f'handling {read1file}, {read2file} ...')
        if df is None:
            logging.debug(f'making new read DF...')
            df = pd.DataFrame(columns=['read1_seq', 'read2_seq'])
            logging.debug(f'handling {read1file}')
            df['read1_seq'] = read_fastq_sequence_pd(read1file, r1s, r1e )
            logging.debug(f'handling {read2file}')
            df['read2_seq'] = read_fastq_sequence_pd(read2file, r2s, r2e )
        else:
            logging.debug(f'making additional read DF...')
            ndf = pd.DataFrame(columns=['read1_seq', 'read2_seq'], dtype="string[pyarrow]")
            logging.debug(f'handling {read1file}')
            ndf['read1_seq'] = read_fastq_sequence_pd(read1file, r1s, r1e )
            logging.debug(f'handling {read2file}')
            ndf['read2_seq'] = read_fastq_sequence_pd(read2file, r2s, r2e )
            logging.debug(f'appending dataframes...')
            df = pd.concat([df, ndf], copy=False, ignore_index=True)
    
    of = f'{outdir}/read1read2.tsv'
    logging.debug(f'writing read1/2 TSV {of} for QC.')
    df.to_csv(of, sep='\t')
  
    df['sequence'] = df['read1_seq'] + df['read2_seq']
    df.drop(['read1_seq','read2_seq'], inplace=True, axis=1)
    
    of = f'{outdir}/fullread.tsv'
    logging.debug(f'writing fullread TSV {of} for QC.')
    df.to_csv(of, sep='\t')
    
    logging.info(f'pulling out MAPseq fields...')
    df['vbc_read'] = df['sequence'].str.slice(0,30)
    df['spikeseq'] = df['sequence'].str.slice(24,32)
    df['libtag'] = df['sequence'].str.slice(30,32)    
    df['umi'] = df['sequence'].str.slice(32,44)
    df['ssi'] = df['sequence'].str.slice(44,52)
    logging.info(f'df done. len={len(df)} returning...')
    sh.add_value('/fastq','reads_handled', len(df) )
    return df


def process_fastq_pairs_pd_agg(infilelist, 
                                    outdir,                         
                                    force=False, 
                                    cp = None):
    '''
    only parse out read lines to pandas, then join with pandas. 
    also AGGREGATE by identical read to minimize data. read_count column is counts of full sequence. 

    Output:
        sequence    read_count
    
    '''
    r1s = int(cp.get('fastq','r1start'))
    r1e = int(cp.get('fastq','r1end'))
    r2s = int(cp.get('fastq','r2start'))
    r2e = int(cp.get('fastq','r2end'))
    logging.debug(f'read1[{r1s}:{r1e}] + read2[{r2s}:{r2e}]')
    df = None
    sh = get_default_stats()
    for (read1file, read2file) in infilelist:
        logging.info(f'handling {read1file}, {read2file} ...')
        if df is None:
            logging.debug(f'making new read DF...')
            df = pd.DataFrame(columns=['read1_seq', 'read2_seq'])
            logging.debug(f'handling {read1file}')
            df['read1_seq'] = read_fastq_sequence_pd(read1file, r1s, r1e )
            logging.debug(f'handling {read2file}')
            df['read2_seq'] = read_fastq_sequence_pd(read2file, r2s, r2e )
        else:
            logging.debug(f'making additional read DF...')
            ndf = pd.DataFrame(columns=['read1_seq', 'read2_seq'], dtype="string[pyarrow]")
            logging.debug(f'handling {read1file}')
            ndf['read1_seq'] = read_fastq_sequence_pd(read1file, r1s, r1e )
            logging.debug(f'handling {read2file}')
            ndf['read2_seq'] = read_fastq_sequence_pd(read2file, r2s, r2e )
            logging.debug(f'appending dataframes...')
            df = pd.concat([df, ndf], copy=False, ignore_index=True)
    
    #of = f'{outdir}/read1read2.tsv'
    #logging.debug(f'writing read1/2 TSV {of} for QC.')
    #df.to_csv(of, sep='\t')
 
    df['sequence'] = df['read1_seq'] + df['read2_seq']
    df.drop(['read1_seq','read2_seq'], inplace=True, axis=1)
    sh.add_value('/fastq','reads_handled', len(df) )  
    
    df = aggregate_reads_pd(df, pcolumn='sequence')  
    of = f'{outdir}/fullreadcounts.tsv'
    logging.debug(f'writing fullreadcounts TSV {of} for QC.')
    df.to_csv(of, sep='\t')
    return df    



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


def split_spike_real_lone_barcodes(config, df):
    '''
    df has  sequence  counts
    should be length 32 ( 30 + YY ) Y= C or T

    # T or C = YY
    # last 2nt of spike-ins AND reals
    
    # A or G = RR
    # last 2nt of L1 controls 
    
    spikeinregex= CGTCAGTC$
    realregex = [TC][TC]$
    loneregex = [AG][AG]$
          
    '''
    #  df[df["col"].str.contains("this string")==False]
    sire = config.get('ssifasta', 'spikeinregex')
    realre = config.get('ssifasta','realregex')
    lonere = config.get('ssifasta', 'loneregex')
    
    logging.debug(f'before filtering: {len(df)}')   
    logging.debug(f"spike-in regex = '{sire}' ")
    simap = df['sequence'].str.contains(sire, regex=True) == True
        
    spikedf = df[simap].copy()
    spikedf.reset_index(inplace=True, drop=True)
    
    remaindf = df[~simap]
    logging.debug(f'spikeins={len(spikedf)} remaindf={len(remaindf)}')
    
    # split real from L1s, and track any that fail both. 
    logging.debug(f"realre = '{realre}' lonere = '{lonere}' ")
    realmap = remaindf['sequence'].str.contains(realre, regex=True) == True
    realdf = remaindf[realmap].copy()
    realdf.reset_index(inplace=True, drop=True)
    
    # remove reals from input. 
    remaindf = remaindf[~realmap]
    logging.debug(f'realdf={len(realdf)} remaindf={len(remaindf)}')
        
    lonemap = remaindf['sequence'].str.contains(lonere, regex=True) == True 
    lonedf = remaindf[lonemap].copy()
    lonedf.reset_index(inplace=True, drop=True)
    logging.debug(f'lonedf={len(lonedf)} remaindf={len(remaindf)}')
    
    unmatcheddf = remaindf[~lonemap].copy()
    unmatcheddf.reset_index(inplace=True, drop=True)
        
    logging.info(f'initial={len(df)} spikeins={len(spikedf)} real={len(realdf)} lone={len(lonedf)} unmatched={len(unmatcheddf)}')    
    return (spikedf, realdf, lonedf, unmatcheddf)


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
def load_readstsv(infile):
    '''
    handle reads output of process_fastq_pairs/aggregate/filter.
     
    '''
    logging.debug(f'loading reads TSV from {infile}')
    STR_COLUMNS = [ 'vbc_read','spikeseq', 'ssi','umi','libtag']
    INT_COLUMNS = ['read_count']
    df = pd.read_csv(infile, sep='\t', index_col=0)
    logging.debug(f'dtypes={df.dtypes}')
    for col in STR_COLUMNS:
        logging.debug(f'converting {col} to string[pyarrow]')
        try:
            df[col] = df[col].astype('string[pyarrow]')
        except KeyError:
            logging.warning(f'no {col} column. continue...')
        
    for col in INT_COLUMNS:
        logging.debug(f'converting {col} to category')
        try:
            df[col] = df[col].astype('uint32')    
        except KeyError:
            logging.warning(f'no {col} column. continue...')
            
    log_objectinfo(df, 'reads-df')
    return df


def load_readtable(infile):
    '''
    very large CSV/TSV files cause some OS to kill. 
    'vbc_read_col', 'libtag','umi','ssi', 'read_count', 'label','rtprimer','type', 'brain', 'region', 'site']
    
    '''
    logging.debug(f'loading readtable TSV from {infile}')
    STR_COLUMNS = [ 'ssi','umi','libtag','label','brain','region','site','vbc_read_col']
    CAT_COLUMNS = ['label','site','type','brain','region','rtprimer']
    df = pd.read_csv(infile, sep='\t', index_col=0)
    logging.debug(f'dtypes={df.dtypes}')
    for col in STR_COLUMNS:
        logging.debug(f'converting {col} to string[pyarrow]')
        df[col] = df[col].astype('string[pyarrow]')

    for col in CAT_COLUMNS:
        logging.debug(f'converting {col} to category')
        df[col] = df[col].astype('category')    
    log_objectinfo(df, 'readtable-df')
    return df

def load_vbctable(infile):
    '''
    very large CSV/TSV files cause some OS to kill. 
    vbc_read_col    label    type    umi_count    read_count    brain    region    site   
    '''
    logging.debug(f'loading vbctable TSV from {infile}')
    STR_COLUMNS = [ 'vbc_read_col']
    CAT_COLUMNS = ['label','site','type','brain','region']
    INT_COLUMNS = ['read_count','umi_count']
    df = pd.read_csv(infile, sep='\t', index_col=0)
    logging.debug(f'dtypes={df.dtypes}')
    for col in STR_COLUMNS:
        logging.debug(f'converting {col} to string[pyarrow]')
        df[col] = df[col].astype('string[pyarrow]')

    for col in CAT_COLUMNS:
        logging.debug(f'converting {col} to category')
        df[col] = df[col].astype('category')

    for col in INT_COLUMNS:
        logging.debug(f'converting {col} to category')
        df[col] = df[col].astype('uint32')
              
    log_objectinfo(df, 'vbctable-df')
    return df


def load_collapse(infile):
    '''
    very large CSV/TSV files cause some OS to kill. 
    'vbc_read_col', 'libtag','umi','ssi', 'read_count', 'label','rtprimer','type', 'brain', 'region', 'site']
    
    For 385610984 rows:   115GB -> 40GB
    
    '''
    logging.debug(f'loading collapse table TSV from {infile}')
    STR_COLUMNS = [ 'spikeseq', 'libtag', 'umi',  'ssi', 'vbc_read_col']
    INT_COLUMNS = ['read_count']
    df = pd.read_csv(infile, sep='\t', index_col=0)
    logging.debug(f'dtypes={df.dtypes}')
    for col in STR_COLUMNS:
        logging.debug(f'converting {col} to string[pyarrow]')
        df[col] = df[col].astype('string[pyarrow]')

    for col in INT_COLUMNS:
        logging.debug(f'converting {col} to integer')
        df[col] = df[col].astype('uint32')    
    log_objectinfo(df, 'collapse-df')
    return df

def oldsnippet():      
    for (read1file, read2file) in infilelist:
        logging.info(f'handling {read1file}, {read2file} ...')
        if df is None:
            logging.debug(f'making new read DF...')
            df = pd.DataFrame(columns=['read1_seq', 'read2_seq'])
            logging.debug(f'handling {read1file}')
            df['read1_seq'] = read_fastq_sequence_pd(read1file, r1s, r1e )
            logging.debug(f'handling {read2file}')
            df['read2_seq'] = read_fastq_sequence_pd(read2file, r2s, r2e )
        else:
            logging.debug(f'making additional read DF...')
            ndf = pd.DataFrame(columns=['read1_seq', 'read2_seq'], dtype="string[pyarrow]")
            logging.debug(f'handling {read1file}')
            ndf['read1_seq'] = read_fastq_sequence_pd(read1file, r1s, r1e )
            logging.debug(f'handling {read2file}')
            ndf['read2_seq'] = read_fastq_sequence_pd(read2file, r2s, r2e )
            logging.debug(f'appending dataframes...')
            df = pd.concat([df, ndf], copy=False, ignore_index=True)
    
    #of = f'{outdir}/read1read2.tsv'
    #logging.debug(f'writing read1/2 TSV {of} for QC.')
    #df.to_csv(of, sep='\t')
 
    df['sequence'] = df['read1_seq'] + df['read2_seq']
    df.drop(['read1_seq','read2_seq'], inplace=True, axis=1)
    sh.add_value('/fastq','reads_handled', len(df) )  
    
    df = aggregate_reads_pd(df, pcolumn='sequence')  
    of = f'{outdir}/fullreadcounts.tsv'
    logging.debug(f'writing fullreadcounts TSV {of} for QC.')
    df.to_csv(of, sep='\t')
    return df    

def load_seqtsv_dd(infile):
    '''
    handle reads output of process_fastq_pairs.
    52nt string sequence 
    USE DASK as dataframe type. 
    
    '''
    logging.debug(f'loading reads TSV from {infile}')
    STR_COLUMNS = [ 'vbc_read','spikeseq', 'ssi','umi','libtag']
    INT_COLUMNS = ['read_count']
    df = pd.read_csv(infile, sep='\t', index_col=0)
    logging.debug(f'dtypes={df.dtypes}')
    for col in STR_COLUMNS:
        logging.debug(f'converting {col} to string[pyarrow]')
        try:
            df[col] = df[col].astype('string[pyarrow]')
        except KeyError:
            logging.warning(f'no {col} column. continue...')
        
    for col in INT_COLUMNS:
        logging.debug(f'converting {col} to category')
        try:
            df[col] = df[col].astype('uint32')    
        except KeyError:
            logging.warning(f'no {col} column. continue...')
            
    log_objectinfo(df, 'reads-df')
    return df


def read_fastq_sequence(infile, start=0, end=-1 ):
    '''
    pull out sequence line, returns series.  
    dtype="string[pyarrow]"
    '''
    logging.info(f'handling FASTQ {infile}')
    slist = []
    interval = 40000000
    if infile.endswith('.gz'):
        fh = gzip.open(infile, "rt")
    else:
        fh = open(infile, 'rt')
    try:
        i = 0
        for line in fh:
            if i % 4 == 1:
                slist.append(line[start:end])  # strip linefeed. 
            if i % interval == 0:
                logging.info(f'sequence {int(i/4)}')
            i += 1
    except:
        pass
    finally:
        fh.close()
    log_objectinfo(slist, 'slist')
    ser = pd.Series(slist, dtype="string[pyarrow]")
    logging.debug(f'series dtype={ser.dtype}')
    log_objectinfo(ser, 'series')
    logging.info(f'done. {len(ser)} sequences extracted.')    
    return ser

def read_fastq_sequence_pd(infile, start=0, end=-1 ):
    '''
    pull out sequence line, returns series.  
    dtype="string[pyarrow]"
    '''
    logging.info(f'handling FASTQ {infile}')
    ser = None
    if infile.endswith('.gz'):
        fh = gzip.open(infile, "rt")
        logging.debug('handling gzipped file...')
    else:
        fh = open(infile, 'rt')
    
    try:
        logging.info(f'reading sequence lines of FASTQ file')
        # investigate engine=pyarrow  (once it supports skiprows lambda)
        # investigate low_memory = False for high-memory nodes
        # investigate memory_map effects on performance.         
        df = pd.read_csv(fh, header=None, skiprows = lambda x: x % 4 != 1, dtype="string[pyarrow]")
        logging.info(f'got sequence df len={len(df)} slicing...')
        df[0] = df[0].str.slice(start, end)
        logging.debug(f'done. pulling series...') 
        ser = df[0]
        logging.debug(f'series defined.')
    except Exception as ex:
        logging.warning(f'exception thrown: {ex} ')
        logging.info(traceback.format_exc(None))

    finally:
        fh.close()
        
    logging.debug(f'series dtype={ser.dtype}')
    log_objectinfo(ser, 'series')
    logging.info(f'done. {len(ser)} sequences extracted.')    
    return ser

def align_collapse_split_pd(df,
                      column='sequence',
                      n_column='vbc_read',
                      p_column='read_count', 
                      max_mismatch=None,
                      n_bases=None, 
                      outdir=None, 
                      datestr=None, 
                      cp=None):
    '''
    Assumes dataframe with unique sequence and read_count columns. 
    Use read_count to choose parent sequence.
    Split out initial <vbc_bases> for alignment/collapse. Save as 'ncolumn'
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
    
    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    
    if n_bases is None:
        n_bases = int(cp.get('collapse', 'n_bases'))
    else:
        n_bases = int(n_bases)
    
    if outdir is None:
        outdir = './'
    outdir = os.path.abspath(outdir)    
    os.makedirs(outdir, exist_ok=True)
    logging.debug(f'collapse: aligner={aligner} n_bases={n_bases} max_mismatch={max_mismatch} outdir={outdir}')    
    
    sh = get_default_stats()      
    sh.add_value('/collapse','n_full_sequences', len(df) )
    sh.add_value('/collapse','n_bases', n_bases )

    # split out vbc_read 
    df[n_column] = df[column].str.slice(0,n_bases)

    # get reduced dataframe of unique head sequences
    logging.info('Getting unique DF...')    
    #udf = pd.DataFrame(df['sequence'].unique(), columns=['sequence'])
    udf = df[n_column].value_counts().reset_index() 
    sh.add_value('/collapse','n_unique_sequences', len(udf) )    



    of = os.path.join( outdir , f'{n_column}.unique.tsv')
    logging.info(f'Writing unique DF to {of}')
    udf.to_csv(of, sep='\t') 
    
    of = os.path.join( outdir , f'{n_column}.unique.fasta')
    logging.info(f'Writing uniques as FASTA to {of}')
    seqfasta = write_fasta_from_df(udf, outfile=of, sequence=[n_column])    

    of = os.path.join(outdir, f'{n_column}.fulldf.tsv')
    logging.info(f'Writing slimmed full DF to {of}')    
    df.to_csv(of, sep='\t', columns=[n_column, p_column])

    # run allXall bowtiex
    of = os.path.join( outdir , f'unique_sequences.bt2.sam')
    logging.info(f'Running {aligner} on {seqfasta} file to {of}')
    btfile = run_bowtie(cp, seqfasta, of, tool=aligner)
    logging.info(f'Bowtie done. Produced output {btfile}. Creating btdf dataframe...')
    btdf = make_bowtie_df(btfile, max_mismatch=max_mismatch, ignore_self=True)
    of = os.path.join( outdir , f'unique_sequences.btdf.tsv')
    btdf.to_csv(of, sep='\t') 
    sh.add_value('/collapse','n_bowtie_entries', len(btdf) )

    # perform collapse...      
    logging.info('Calculating Hamming components...')
    edgelist = edges_from_btdf(btdf)
    btdf = None  # help memory usage
    sh.add_value('/collapse','n_edges', len(edgelist) )
    
    of = os.path.join( outdir , f'edgelist.txt')
    writelist(of, edgelist)
    logging.debug(f'edgelist len={len(edgelist)}')
    
    components = get_components(edgelist)
    logging.debug(f'all components len={len(components)}')
    sh.add_value('/collapse','n_components', len(components) )
    of = os.path.join( outdir , f'components.txt')
    writelist(of, components)
    edgelist = None  # help memory usage    

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
    newdf = collapse_by_components_pd(df, udf, mcomponents, column=n_column, pcolumn=p_column, outdir=outdir)
    # newdf has sequence and newsequence columns, rename to orig_seq and sequence
    newcol = f'{n_column}_col'        
    newdf.rename( {'newsequence': newcol}, inplace=True, axis=1)    
    logging.info(f'Got collapsed DF. len={len(newdf)}')

    rdf = newdf[[ n_column, newcol, 'read_count' ]]
    of = os.path.join( outdir , f'read_collapsed.tsv')
    logging.info(f'Writing reduced mapping TSV to {of}')
    rdf.to_csv(of, sep='\t')
    
    # newdf.drop(column, inplace=True, axis=1)
    #joindf = pd.DataFrame( newdf['new_seq'] + newdf['tail'], columns=['sequence'])
    #of = os.path.join( outdir , f'collapsed.fasta')
    #logging.info(f'Writing collapsed fasta to {of}')
    #write_fasta_from_df(joindf, of)        
    return newdf   


def filter_non_inj_umi_merge(rtdf, ridf, inj_min_umi=1, write_out=False):
    '''
    rtdf and ridf should already be filtered by brain, type, and anything else that might complicate matters.
    remove rows from rtdf that do not have at least <min_injection> value in the row 
    of ridf with the same index (VBC sequence)
    Does an inner join() on the dataframes, keyed on sequence. 
    Keeps values and columns from first argument (rtdf)
    
    '''
    logging.info(f'filtering non-injection. inj_min_umi={inj_min_umi}')
    logging.debug(f'before threshold inj df len={len(ridf)}')
    ridf = ridf[ridf.umi_count >= inj_min_umi]
    ridf.reset_index(inplace=True, drop=True)
    logging.debug(f'before threshold inj df len={len(ridf)}')   
    
    if write_out:
        ridf.to_csv('./ridf.tsv', sep='\t')
        rtdf.to_csv('./rtdf.tsv', sep='\t')
    
    mdf = pd.merge(rtdf, ridf, how='inner', left_on='vbc_read_col', right_on='vbc_read_col')
    incol = mdf.columns
    outcol = []
    selcol =[]
    for c in incol:
        if not c.endswith('_y'):
            selcol.append(c)
            outcol.append(c.replace('_x',''))
    mdf = mdf[selcol]
    mdf.columns = outcol
    logging.debug(f'before drop_duplicates. len={len(mdf)}')
    mdf.drop_duplicates(inplace=True)
    mdf.reset_index(drop=True, inplace=True)
    logging.debug(f'after drop_duplicates. len={len(mdf)} columns={mdf.columns}')
    logging.debug(f'created merged/joined DF w/ common sequence items.  df=\n{mdf}')
    return mdf

def findrepeats(s, n=2):
    '''
    Finds strings of repeated characters longer than <num>
    Typically used to find homopolymers in sequences for exclusion from results
    as repeats interfere with amplification/replication.
    '''
    rcount = 1
    last = None
    res = False
    for c in s:
        if last is None:
            last = c
        else:
            if c == last:
                rcount += 1
            else:
                last = c
                rcount = 1
        if rcount >= n:
            # short-circuit if found. 
            return True
        else:
            pass
    return res


def remove_base_repeats(df, col='sequence', n=7):
    '''
    removes rows where column string repeats C,G,A, or T <n> times or more.    
    '''
    startlen=len(df)
    logging.debug(f'removing n>{n} repeats from df len={startlen}')
    bases = ['C','G','A','T']
    for b in bases:
        pat = b*n
        logging.debug(f'searching for {pat}')
        df = df[ ~df[col].str.contains(pat)]
    
    endlen=(len(df))
    logging.debug(f'df without {n} repeats len={endlen}')
    df.reset_index(drop=True, inplace=True)
    return df


def has_base_repeats(seqstring, n=7):
    '''
    if string repeats C,G,A, or T <n> times or more.    
    '''
    bases = ['C','G','A','T']
    for b in bases:
        pat = b*n
        if pat in seqstring:
            return True
    return False

def has_n_bases(seqstring, n=0):
    '''
    if string contains 'N' <n> times or more.    
    '''
    ncount = seqstring.count('N')
    if ncount > n:
        return True
    return False
 
 
def gzip_decompress(filename):
    """
    default for copyfileobj is 16384
    https://blogs.blumetech.com/blumetechs-tech-blog/2011/05/faster-python-file-copy.html

    """
    log = logging.getLogger('utils')
    if filename.endswith('.gz'):
        targetname = filename[:-3]
        bufferlength = 10 * 1024 * 1024  # 10 MB
        with gzip.open(filename, 'rb') as f_in:
            with open(targetname, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out, length=bufferlength)
    else:
        log.warn(
            f'tried to gunzip file without .gz extension {filename}. doing nothing.')

def load_df(filepath):
    """
    Convenience method to load DF consistently across modules. 
    """
    logging.debug(f'loading {filepath}')
    filepath = os.path.expanduser(filepath)
    df = pd.read_csv(filepath, sep='\t', index_col=0, keep_default_na=False, dtype="string[pyarrow]", comment="#")
    #df.fillna(value='', inplace=True)
    logging.debug(f'initial load done. converting types...')
    df = df.convert_dtypes(convert_integer=False)
    for col in df.columns:
        logging.debug(f'trying column {col}')
        try:
            df[col] = df[col].astype('uint32')
        except ValueError:
            logging.debug(f'column {col} not int')
    logging.debug(f'{df.dtypes}')
    return df


def merge_fastq_pairs(config, readfilelist, outdir):
    logging.debug(f'processing {readfilelist}')
    if outdir is None:
        outdir = "."
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            logging.debug(f'made outdir={outdir}')
    pairshandled = 0
    for (read1file, read2file) in readfilelist:
        logging.debug(f'read1file = {read1file}')
        logging.debug(f'read2file = {read2file}')
        pairshandled += 1 

def filter_split_pd(df, 
                    max_repeats=None,
                    max_n_bases=None,
                    column='sequence',
                    drop=True,
                    cp=None
                    ):
    '''
    filter sequence by repeats, Ns. 
    split into MAPseq fields
    optionally drop sequence column
    
    '''
    if cp is None:
        cp = get_default_config()

    logging.info(f'Filtering by read quality. Repeats. Ns.')
    df = filter_reads_pd(df, 
                       max_repeats=max_repeats,
                       max_n_bases=max_n_bases, 
                       column=column )
    logging.info(f'Splitting into mapseq fields. ')    
    df = split_mapseq_fields(df, 
                             drop=True)
    return df

def make_read_countplot(config, df, outfile, title=None ): 
    '''
    makes individual read count plot from sequence read_count DF 
    assumes 'sequence' and 'read_count' columns. 
    
    '''   
    plt.set_loglevel (level = 'warning')
    
    logging.debug(f'handling sequence df len={len(df)}')
    outfile = os.path.abspath(outfile)    

    if title is None:
        title='Read count frequence plot.'

    df.sort_values(by='read_count', ascending=False, inplace=True)
    df.reset_index(inplace=True)
    df['index'] = df['index'].astype('int64')

    plt.figure()
    plt.plot(np.log10(df['index']), np.log10(df['read_count']))
    plt.title(title)
    plt.xlabel("log10(index)")
    plt.ylabel("log10(read_count)")
    logging.info(f'Saving count plot to {outfile}')
    plt.savefig(outfile)


def make_read_countplot_sns(cp, df, outfile='count-frequency-plot.pdf' ):
    '''
    makes two figures, one log-normalized, one straight for same counts df. 
    '''
    from matplotlib.backends.backend_pdf import PdfPages as pdfpages
    df.sort_values(by='read_count', ascending=False, inplace=True)
    df.reset_index(inplace=True)
    df['index'] = df.index.astype('int64')

    #page_dims = 
    #  A4 landscape   (11.69,8.27) 
    #  A3 landscape  (16.53,11.69)
    #  A4 portrait  (8.27, 11.69)  
    page_dims = (8.27, 11.69)
    with pdfpages(outfile) as pdfpages:
        axlist = []
        fig,axes = plt.subplots(nrows=2, ncols=1, figsize=page_dims, layout='constrained') 
        fig.suptitle(f'Read counts frequency plots.')
        for a in axes.flat:
            axlist.append(a)
        ax = axlist[0]
        counts_axis_plot_sns(ax, df, scale=None)
        ax = axlist[1]
        counts_axis_plot_sns(ax, df, scale='log10')
        
        pdfpages.savefig(fig)




def make_shoulder_plot_sns(df, 
                           title='frequency plot', 
                           site='target', 
                           outfile='frequency-plot.pdf', 
                           column='read_count' ):
    '''
    makes two figures, one log-normalized, one straight for same counts df. 
    '''
    from matplotlib.backends.backend_pdf import PdfPages as pdfpages
    
    cdf = pd.DataFrame( df[df['site'] == site ][column].copy(), columns=[column])
    cdf.sort_values(by=column, ascending=False, inplace=True)
    cdf.reset_index(inplace=True,  drop=True)
    cdf['index'] = cdf.index.astype('int64')

    #page_dims = 
    #  A4 landscape   (11.69,8.27) 
    #  A3 landscape  (16.53,11.69)
    #  A4 portrait  (8.27, 11.69)  
    page_dims = (8.27, 11.69)
    with pdfpages(outfile) as pdfpages:
        axlist = []
        fig,axes = plt.subplots(nrows=2, ncols=1, figsize=page_dims, layout='constrained') 
        fig.suptitle(f'{column} freq: site={site}')
        for a in axes.flat:
            axlist.append(a)
        ax = axlist[0]
        counts_axis_plot_sns(ax, cdf, scale=None, column=column)
        ax = axlist[1]
        counts_axis_plot_sns(ax, cdf, scale='log10', column=column)        
        pdfpages.savefig(fig)

def make_shoulder_plots(df, outdir=None, cp=None):
    # make shoulder plots. injection, target
    logging.info('making shoulder plots...')
    if outdir is None:
        outdir = os.path.abspath('./')
    logging.getLogger('matplotlib.font_manager').disabled = True
    if len(df[df['site'] == 'injection'] ) > 1:
        make_shoulder_plot_sns(df, site='injection', outfile=f'{outdir}/inj-counts.pdf')
    else:
        logging.info(f'no injection sites, so no plot.')
    if len(df[df['site'] == 'target'] ) > 1:    
        make_shoulder_plot_sns(df, site='target', outfile=f'{outdir}/target-counts.pdf')   
    else:
        logging.info(f'no target sites, so no plot.') 
        
        
def process_make_matrices_pd(df,
                          outdir=None,
                          exp_id = 'M001',  
                          inj_min_umi = None,
                          target_min_umi = None,
                          target_min_umi_absolute = None, 
                          label_column='label',
                          cp = None):
    '''
    -- per brain, pivot real VBCs (value=umi counts)
    -- create real, real normalized by spikein, and  
    -- use label_column to pivot on, making it the y-axis, x-axis is vbc sequence.  
    -- optionally require VBCs to be in injection to be included
    -- optionally include injections in matrices.
    
    theshold logic. 
    inj_min_umi                VBC UMI must exceed to be kept.
    target_min_umi             if ANY target area exceeds, keep all of that VBC targets. 
    target_min_umi_absolute    hard threshold cutoff
 

    '''
    if cp is None:
        cp = get_default_config()
    if outdir is None:
        outdir = os.path.abspath('./')
    require_injection = cp.getboolean('matrices','require_injection')
    include_injection = cp.getboolean('matrices','include_injection')

    max_negative = 1
    max_water_control = 1

    if inj_min_umi is None:
        inj_min_umi = int(cp.get('matrices','inj_min_umi'))
    if target_min_umi is None:
        target_min_umi = int(cp.get('matrices','target_min_umi'))   
    if target_min_umi_absolute is None:
        target_min_umi_absolute = int(cp.get('matrices','target_min_umi_absolute'))

    use_target_negative=cp.getboolean('matrices','use_target_negative')
    use_target_water_control=cp.getboolean('matrices','use_target_water_control')    
    clustermap_scale = cp.get('matrices','clustermap_scale')
    logging.debug(f'running exp_id={exp_id} inj_min_umi={inj_min_umi} target_min_umi={target_min_umi} use_target_negative={use_target_negative} use_target_water_control={use_target_water_control}')

    df['brain'] = df['brain'].astype('string')
    bidlist = list(df['brain'].dropna().unique())
    bidlist = [ x for x in bidlist if len(x) > 0 ]
    bidlist.sort()
    logging.debug(f'handling brain list: {bidlist}')
    
    sh = get_default_stats()

    norm_dict = {}

    for brain_id in bidlist:
        valid = True
        logging.debug(f'handling brain_id={brain_id}')
        bdf = df[df['brain'] == brain_id]

        # extract and threshold injection areas
        idf = bdf[bdf['site'].str.startswith('injection')]
        ridf = idf[idf['type'] == 'real'] 
        if inj_min_umi > 1:
            before = len(ridf)
            ridf = filter_all_lt(ridf, 'vbc_read_col', 'umi_count', inj_min_umi)            
            #
            # INCORRECT: if require_injection is false, there may not be any injection but thats OK.
            # Let the filter_non_inj_umi() function do its thing...
            #
            #if not len(ridf) > 0:
            #    valid = False
            #    logging.warning(f'No injection VBCs passed min_inj_umi filtering! Skip brain.')
            logging.debug(f'filtering by inj_min_umi={inj_min_umi} before={before} after={len(ridf)}')
        else:
            logging.debug(f'inj_min_umi={inj_min_umi} no filtering.')
                   
        # extract and do absolute thesholding of target areas
        tdf = bdf[bdf['site'].str.startswith('target')]
        rtdf = tdf[tdf['type'] == 'real'] 

        if target_min_umi_absolute > 1:
            before = len(rtdf)
            rtdf = rtdf[rtdf['umi_count'] >= target_min_umi_absolute ]
            rtdf.reset_index(drop=True, inplace=True)
            after = len(rtdf)
            logging.debug(f'filtering by target_min_umi_absolute={target_min_umi_absolute} before={before} after={after}')

        # threshold by min_target or threshold by target-negative
        # if use_target_negative is true, but no target negative site 
        # defined, use min_target and throw warning. 
        if use_target_negative:
            logging.info(f'use_target_negative is {use_target_negative}')
            max_negative = calc_min_umi_threshold(bdf, 'target-negative', cp)
            logging.debug(f'target-negative UMI count = {max_negative}')

        if use_target_water_control:
            logging.info(f'use_target_water_control is {use_target_water_control}')
            max_water_control = calc_min_umi_threshold(bdf, 'target-water-control',cp)
            logging.debug(f'target_water_control UMI count = {max_water_control}')

        target_min_umi = max([target_min_umi, max_negative, max_water_control ])
        logging.debug(f'min_target UMI count after all constraints = {target_min_umi}')   

        # min_target is now either calculated from target-negative, or from config. 
        if target_min_umi > 1:
            before = len(rtdf)
            frtdf = filter_all_lt(rtdf, 'vbc_read_col', 'umi_count', target_min_umi)            
            if not len(frtdf) > 0:
                valid = False
                logging.warning(f'No VBCs passed min_target filtering! Skip brain.')
            logging.debug(f'filtering by target_min_umi={target_min_umi} before={before} after={len(rtdf)}')
        else:
            logging.debug(f'target_min_umi={target_min_umi} no filtering.')
            frtdf = rtdf
  
        # get injection-filtered real target table, and target-filtered real injection table
        # in case either is needed. 
        (ifrtdf, tfridf) = filter_non_inj_umi(frtdf, ridf, inj_min_umi=inj_min_umi)            
        logging.debug(f'{len(ifrtdf)} real target VBCs after injection filtering.')
        logging.debug(f'{len(tfridf)} real injection VBCs after target filtering.')

        if require_injection:
            logging.debug(f'require_injection={require_injection} inj_min_umi={inj_min_umi}')
            if len(ridf) == 0:
                logging.warning('require_injection=True but no real VBCs from any injection site.')
                valid = False
            elif not len(ifrtdf) > 0:
                logging.warning(f'No VBCs passed injection filtering! Skip brain.')
                valid = False
            else:
                # set filtered real to injection-filtered real targets. 
                frtdf = ifrtdf    
        else:
            # frtdf remains unfiltered...
            logging.debug(f'require_injection={require_injection} proceeding...')

        vbcdf = frtdf
        
        if include_injection:
            logging.debug(f'include_injection={include_injection} including mutually present injection VBCs') 
            vbcdf = pd.concat( [vbcdf, tfridf], ignore_index=True ) 

        # make matrices if per-brain data is valid... 
        if valid:                               
            vbcdf.to_csv(f'{outdir}/{brain_id}.vbctable.tsv', sep='\t')
            vbcdf.to_parquet(f'{outdir}/{brain_id}.vbctable.parquet')
            
            target_columns = list( vbcdf[vbcdf['site'].str.startswith('target')]['label'].unique())
            injection_columns = list( vbcdf[vbcdf['site'].str.startswith('injection')]['label'].unique())
            
            rbcmdf = vbcdf.pivot(index='vbc_read_col', columns=label_column, values='umi_count')
            scol = natsorted(list(rbcmdf.columns))
            rbcmdf = rbcmdf[scol]
            rbcmdf.fillna(value=0, inplace=True)
            logging.debug(f'brain={brain_id} raw real barcode matrix len={len(rbcmdf)}')            
            
            # filtered reals
            fbcmdf = vbcdf.pivot(index='vbc_read_col', columns=label_column, values='umi_count')
            scol = natsorted(list(fbcmdf.columns))
            fbcmdf = fbcmdf[scol]
            fbcmdf.fillna(value=0, inplace=True)
            logging.debug(f'brain={brain_id} real barcode matrix len={len(fbcmdf)}')

            # spikes
            sdf = tdf[tdf['type'] == 'spike']
            if include_injection: 
                logging.debug(f'include_injection={include_injection} need injection spikes...')
                isdf = idf[idf['type']== 'spike']
                sdf = pd.concat( [sdf, isdf], ignore_index=True  )
                
            sbcmdf = sdf.pivot(index='vbc_read_col', columns=label_column, values='umi_count')
            spcol = natsorted(list(sbcmdf.columns))
            sbcmdf = sbcmdf[spcol]
            sbcmdf.fillna(value=0, inplace=True)    
            logging.debug(f'brain={brain_id} spike barcode matrix len={len(sbcmdf)}')
    
            (fbcmdf, sbcmdf) = sync_columns(fbcmdf, sbcmdf)
            
            if include_injection:
                logging.debug(f'include_injection={include_injection} need grouped normalization.')
                nbcmdf = normalize_weight_grouped(fbcmdf, sbcmdf, columns = [target_columns, injection_columns])
            else:
                logging.debug(f'include_injection={include_injection} ungrouped normalization.')
                nbcmdf = normalize_weight_grouped(fbcmdf, sbcmdf)
            
            # put columns in natural sort order...
            nbcmdf = nbcmdf[ natsorted( list(nbcmdf.columns) ) ]            
            logging.debug(f'nbcmdf.describe()=\n{nbcmdf.describe()}')
            
            norm_dict[brain_id] = nbcmdf
            
            # make scaled normalized matrix
            scbcmdf = normalize_scale(nbcmdf, logscale=clustermap_scale)
            scbcmdf.fillna(value=0, inplace=True)
            logging.debug(f'scbcmdf.describe()=\n{scbcmdf.describe()}')
            
            sh.add_value(f'/matrices/brain_{brain_id}','valid', 'True' )            
            sh.add_value(f'/matrices/brain_{brain_id}','n_vbcs_real', len(rbcmdf) )
            sh.add_value(f'/matrices/brain_{brain_id}','n_vbcs_real_filtered', len(fbcmdf) )
            sh.add_value(f'/matrices/brain_{brain_id}','n_vbcs_spike', len(sbcmdf) )
            
            # raw real
            rbcmdf.to_csv(f'{outdir}/{brain_id}.rbcm.tsv', sep='\t')
            # filtered real
            fbcmdf.to_csv(f'{outdir}/{brain_id}.fbcmdf.tsv', sep='\t')
            # spike-in matrix
            sbcmdf.to_csv(f'{outdir}/{brain_id}.sbcm.tsv', sep='\t')
            # filtered normalized by spike-ins.    
            nbcmdf.to_csv(f'{outdir}/{brain_id}.nbcm.tsv', sep='\t')
            # log-scaled normalized 
            scbcmdf.to_csv(f'{outdir}/{brain_id}.scbcm.tsv', sep='\t')
                       
        else:
            logging.info(f'brain {brain_id} data not valid.')
            sh.add_value(f'/matrices/brain_{brain_id}','valid', 'False' )
            
        logging.info(f'done with brain={brain_id}')    
    logging.info(f'got dict of {len(norm_dict)} normalized barcode matrices. returning.')
    return norm_dict


def load_sample_info_old(config, file_name, sheet_name='Sample information'):
    #
    # Parses Excel spreadsheet to get orderly sample metadata, saves as sampleinfo.tsv.     
    # OR Reads in sampleinfo.tsv
    # Assumes various properties of spreadsheet that need to stay static. 
    #
    #   ['Tube # by user', 'Our Tube #', 'Sample names provided by user',
    #   'Site information', 'RT primers for MAPseq', 'Brain ', 'Column#']
    #
    # If brain is not given, or is empty, all are set to 'brain1'. 
    # If region is not given, or is empty, all are set to <rtprimer>
    # 
    
    # Mappings for excel columns. 
    sheet_to_sample = {
            'Tube # by user'                  : 'usertube', 
            'Our Tube #'                      : 'ourtube', 
            'Sample names provided by user'   : 'samplename', 
            'Site information'                : 'siteinfo',
            'RT primers for MAPseq'           : 'rtprimer',
            'Brain'                           : 'brain',
            'Region'                          : 'region',
            'Matrix Column'                   : 'matrixcolumn',
        }
    
    sample_columns = ['usertube', 'ourtube', 'samplename', 'siteinfo', 'rtprimer', 'brain', 'region', 'matrixcolumn'] 
    int_sample_col = ['usertube', 'ourtube', 'rtprimer', 'region', 'matrixcolumn']     # brain is sometimes not a number. 
    str_sample_col = ['usertube', 'ourtube', 'samplename', 'siteinfo', 'rtprimer', 'brain' ,'region']

    if file_name.endswith('.xlsx'):
        edf = pd.read_excel(file_name, sheet_name=sheet_name, header=1)        
        sdf = pd.DataFrame()
        
        for ecol in edf.columns:
            ecol_stp = ecol.strip()    
            try:
                # map using stripped column name, retrieve using actual excel column name
                # which may have trailing spaces...
                scol = sheet_to_sample[ecol_stp]
                logging.debug(f'found mapping {ecol} -> {scol}')
                cser = edf[ecol]
                #logging.debug(f'column for {scol}:\n{cser}')
                sdf[scol] = cser
            
            except KeyError:
                logging.debug(f'no mapping for {ecol} continuing...')
            
            except Exception as ex:
                logging.error(f'error while handling {ecol} ')
                logging.info(traceback.format_exc(None))
                
        #logging.debug(f"rtprimer = {sdf['rtprimer']} dtype={sdf['rtprimer'].dtype}")
        
        # Only keep rows with rtprimer info. 
        sdf = sdf[sdf['rtprimer'].isna() == False]

        # Fix brain column
        sdf.loc[ sdf.brain.isna(), 'brain'] = 0
        try:
            sdf.brain = sdf.brain.astype('int')
        except ValueError:
            pass
        sdf.brain = sdf.brain.astype('string')
        sdf.loc[ sdf.brain == '0', 'brain'] = ''

        # fix rtprimer column, if it was read as float string (e.g '127.0' )
        #sdf['rtprimer'] =  sdf['rtprimer'].astype(float).astype(int).astype(str)
        
        for scol in sample_columns:
            try:
                ser = sdf[scol]
            except KeyError as ke:
                logging.info(f'no column {scol}, required. Creating...')
                if scol == 'samplename':
                    sdf[scol] = sdf['ourtube']
                elif scol == 'region':
                    sdf[scol] = sdf['rtprimer']
        
        # fix empty rows. 
        sdf.brain = sdf.brain.astype(str)
        #sdf.brain[sdf.brain == 'brain-nan'] = ''
        sdf.loc[sdf.brain == 'brain-nan','brain'] = ''
                    
        sdf.replace(r'^s*$', float('NaN'), regex = True, inplace=True)
        sdf.dropna(how='all', axis=0, inplace=True)        
        sdf = fix_columns_int(sdf, columns=int_sample_col)
        sdf = fix_columns_str(sdf, columns=str_sample_col)

    elif file_name.endswith('.tsv'):
        sdf = pd.read_csv(file_name, sep='\t', index_col=0, keep_default_na=False, dtype =str, comment="#")
        #df.fillna(value='', inplace=True)
        sdf = sdf.astype('str', copy=False)    
        sdf = fix_columns_int(sdf, columns=int_sample_col)
    else:
        logging.error(f'file {file_name} neither .xlsx or .tsv')
        sdf = None
        
    logging.debug(f'created reduced sample info df:\n{sdf}')
    return sdf


def read_threshold_all(cp, alldf):
    '''
    we did not threshold BC by BC when processing fastqs, so we can threshold by proxy
    using a ratio of UMIs to read_count behind them. 
    '''
    target_ratio = float( cp.get('ssifasta','target_read_ratio'))
    injection_ratio = float( cp.get('ssifasta','injection_read_ratio'))
    
    alldf['read_ratio'] = alldf['read_count'] / alldf['umi_count']
        
    injdf = alldf[ alldf['site'].str.startswith('injection') ]
    tardf = alldf[ alldf['site'].str.startswith('target') ]
    injdf = injdf[ injdf['read_ratio'] > injection_ratio ]
    tardf = tardf[ tardf['read_ratio'] > target_ratio ]
    df_merged = pd.concat([injdf, tardf], ignore_index = True, sort=False)
    logging.info(f'target_ratio = {target_ratio} injection_ratio = {injection_ratio} before len={len(alldf)} after len={len(df_merged)}')
    return df_merged
       
def make_read_counts_df(config, seqdf, label=None):
    '''
    input dataframe with 'sequence' column
    make counts column for identical sequences.  
    optionally assign a label to set in new column
    
    '''
    logging.debug(f'seqdf=\n{seqdf}')
    ser = seqdf['sequence'].value_counts()
    df = pd.DataFrame(columns=['sequence','read_count'])
    df['sequence'] = ser.index
    df['read_count'] = ser.values
    logging.debug(f'counts df = \n{df}')
    if label is not None:
        df['label'] = label
    return df

def make_umi_counts_df(config, seqdf, label=None):
    '''
    input dataframe with 'sequence' column
    make counts column for identical sequences.  
    optionally assign a label to set in new column
    combine values from collapsed read_counts 
    sequence   read_count  label     
    AAA        23           A
    AAA        5            A
    BBB        7            A
    BBB        6            A
    CCC        5            A
    
    ->
    AAA    read_count  umi_count label
    BBB    13           2          A
    CCC    28           2          A
    CCC    5            1          A
    
    https://stackoverflow.com/questions/73874908/value-counts-then-sum-of-a-different-column
    https://stackoverflow.com/questions/46431243/how-to-get-groupby-sum-of-multiple-columns
     
    '''
    logging.debug(f'seqdf=\n{seqdf}')
    # keep original label
    label = seqdf['label'].unique()[0]
    vdf = seqdf.drop('label',axis=1).groupby('sequence').agg({'read_count':'sum','umi_count':'sum'}).reset_index()
    #vdf = seqdf.groupby('sequence')['read_counts'].agg(count='count',sum='sum').reset_index()
    #vdf.rename(columns={'count': 'umi_counts','sum': 'read_counts'}, inplace=True)
    vdf.sort_values(by=['umi_count'], ascending=False, inplace=True)
    vdf.reset_index(inplace=True, drop=True)
    vdf['label'] = label
    logging.debug(f'counts df = \n{vdf}')
    return vdf




def make_read_counts_dfs(config, filelist, outdir):
    '''
    
    '''
    dflist = []
    for filepath in filelist:
        logging.info(f'calculating counts for file {filepath} ...')    
        dirname = os.path.dirname(filepath)
        
        if outdir is not None:
            dirname = outdir
        
        filename = os.path.basename(filepath)
        (base, ext) = os.path.splitext(filename)   
        logging.debug(f'handling {filepath} base={base}')
        
        # make raw fasta TSV of barcode-splitter output for one barcode. 
        # trim to 44 unique w/ counts. 
        seqdf = make_fasta_df(config, filepath)
        of = os.path.join(dirname , f'{base}.read.seq.tsv')
        seqdf.to_csv(of, sep='\t')
        
        # to calculate threshold we need counts calculated. 
        cdf = make_read_counts_df(config, seqdf, label=base)  
        logging.debug(f'made counts df: {len(cdf)} reads.')
        of = os.path.join(dirname , f'{base}.read.counts.tsv')
        cdf.to_csv(of, sep='\t')
        dflist.append(cdf)
    logging.debug(f'returning list of {len(dflist)} counts DFs...')
    return dflist 


def filter_counts_df(cp, countsdf, min_count):
    '''
    Assumes read_count column
    '''      
    countsdf = countsdf[countsdf['read_count'] >= min_count]
    countsdf.reset_index(drop=True, inplace=True)
    logging.info(f'filtering by count: before={initlen} after={len(seqdf)} min_count={min_count} ')
    return cdf    
            
def make_fasta_df(config, infile, ignore_n=True):
    '''
    input fasta 
    ignore 'N' sequences.
    '''   
    slist = []
    rcs = SeqIO.parse(infile, "fasta")
    handled = 0
    for sr in rcs:
        s = sr.seq
        if ('N' in sr) and ignore_n :
            pass
        else:
            slist.append(str(s))
        handled += 1    
    logging.debug(f"kept {len(slist)} sequences out of {handled}")    
    df = pd.DataFrame(slist, columns=['sequence'] )
    return df

def normalize_scale(df, columns = None, logscale='log2', min=0.0, max=1.0):
    '''
    Log scale whole matrix.   log10 or log2 ???
    Set -inf to 0
    Set NaN to 0
    
    
    '''
    #logging.debug(f'making rows sum to one...')
    if logscale == 'log2':
        ldf = np.log2(df)
    elif logscale == 'log10':
        ldf = np.log10(df)
    
    for c in ldf.columns:
        ldf[c] = np.nan_to_num(ldf[c], neginf=0.0)
    return ldf

def set_siteinfo(df, sampdf, column='sequence', cp=None):
    '''
    This is a purely utility function to determine site type for
    shoulder plots. 
    
    '''
    if cp is None:
        cp=get_default_config()
    ssi_st = int(cp.get('mapseq','ssi_st') )
    ssi_end = int(cp.get('mapseq','ssi_end') )
    bcfile = os.path.expanduser( cp.get('barcodes','ssifile'))

    df['ssi'] = df[column].str.slice(ssi_st,ssi_end).astype('string[pyarrow]')
    df.drop( column, axis=1, inplace=True)

    # set site
    # map SSIs, set unknown to unmatched.
    logging.debug(f'getting rt labels...')
    labels = get_rtlist(sampdf)
    logging.debug(f'rtlabels={labels}')
    #bcdict = get_barcode_dict(bcfile, labels)
    rtdict = get_rtprimer_dict(bcfile, labels)
   
    #logging.info('filling in labels by SSI sequences...')
    #df['label'] = df['ssi'].map(bcdict)
    #logging.info('labelling unmatched...')
    #df.fillna({'label': 'nomatch'}, inplace=True)

    logging.info('filling in rtprimer number by SSI sequences...')    
    df['rtprimer'] = df['ssi'].map(rtdict)
    logging.info('labelling unmatched...')
    df.fillna({'rtprimer': 'nomatch'}, inplace=True)
    
    sdf = sampdf[['rtprimer','siteinfo']]
    sdf = sdf[sdf['siteinfo'] != '']
    smap = dict(zip(sdf['rtprimer'],sdf['siteinfo']))
    df['site'] = df['rtprimer'].map(smap)
    df.fillna({'site': 'nomatch'}, inplace=True)
    sdf = None    
    
    return df

def set_counts_df(seqdf, column='', cp=None):
    '''
    assumed sequence column 
    calculates duplicates and sets read_count column
 
    '''
    initlen = len(seqdf)
    logging.debug(f'setting read counts for sequence DF len={len(seqdf)}')
    cdf = make_read_counts_df(cp, seqdf, label=None)
    fmap = cdf.set_index('sequence')['read_count'].to_dict()
    seqdf['read_count'] = seqdf['sequence'].map(fmap)
    # change made to inbound df, but return anyway
    return seqdf
    

def process_mapseq_all_native(infiles, sampleinfo, outdir=None, force=False, cp=None):    
    '''
    DO NOT USE. NEEDS UPDATING FOR LARGE DATA
    
    performs end-to-end default processing. fastq pairs to matrices
    
    '''
    if cp is None:
        cp = get_default_config()
    bcfile = cp.get('barcodes','ssifile')
    expid = cp.get('project','project_id')
    
    if outdir is None:
        outdir = os.path.abspath('./')
        logging.debug(f'outdir not provided, set to {outdir}')
    else:
        outdir = os.path.abspath(outdir)
        logging.debug(f'outdir = {outdir} ')   
   
    logging.debug(f'ensuring outdir: {outdir} ')
    os.makedirs(outdir, exist_ok=True)
    logging.info(f'outdir={outdir}')

   
    # process_fastq_pairs
    if (len(infiles) < 2)  or (len(infiles) % 2 != 0 ):
        parser.print_help()
        print('error: the following arguments are required: 2 or multiple of 2 infiles')
        sys.exit(1)
    infiles = package_pairs(infiles) 
    
    outfile = f'{outdir}/{expid}.reads.tsv'
    df = process_fastq_pairs_pd(infiles, 
                                 outdir,  
                                 cp=cp)    

    logging.debug(f'filtering by read quality. repeats. Ns.')
    df = filter_reads_pd(df, 
                           column='sequence' )
    logging.debug(f'calc/set read counts on original reads.')    
    df = set_counts_df(df, column='sequence')
    logging.info(f'dropping sequence column to slim.')
    df.drop('sequence', axis=1, inplace=True)    
    logging.info(f'Got dataframe len={len(df)} Writing to {outfile}')
    logging.debug(f'dataframe dtypes:\n{df.dtypes}\n')
    df.to_csv(outfile, sep='\t')
    dir, base, ext = split_path(outfile)
    outfile = os.path.join( dir , f'{base}.parquet')
    logging.info(f'df len={len(df)} as parquet to {outfile}...')
    df.to_parquet(outfile)
    logging.info(f'process_fastq_pairs done.')
    
    # align_collapse    
    outfile = f'{outdir}/{expid}.collapse.tsv'
    df = align_collapse_pd(df, 
                           outdir=outdir, 
                           cp=cp)
    logging.info(f'Saving len={len(df)} as TSV to {outfile}...')
    df.to_csv(outfile, sep='\t')
    dir, base, ext = split_path(outfile)
    outfile = os.path.join(dir, f'{base}.parquet')
    logging.info(f'df len={len(df)} as parquet to {outfile}...')
    df.to_parquet(outfile)    
    logging.info(f'align_collapse done.')
    
    # make_readtable
    outfile = f'{outdir}/{expid}.readtable.tsv'
    logging.debug(f'loading sample DF...')
    sampdf = load_sample_info(cp, sampleinfo)
    logging.debug(f'\n{sampdf}')
    sampdf.to_csv(f'{outdir}/sampleinfo.tsv', sep='\t')
    
    df = process_make_readtable_pd(df,
                                   sampdf,
                                   bcfile=bcfile, 
                                   outdir=outdir, 
                                   cp=cp)

    logging.info(f'Got dataframe len={len(df)} Writing to {outfile}')
    df.to_csv(outfile, sep='\t')
    dir, base, ext = split_path(outfile)
    outfile = os.path.join(dir, f'{base}.parquet')
    logging.info(f'df len={len(df)} as parquet to {outfile}...')
    df.to_parquet(outfile)
    logging.info(f'make_readtable done.')    

    # make_vbctable
    outfile = f'{outdir}/{expid}.vbctable.tsv'
    df = process_make_vbctable_pd(df,
                               outdir=outdir,
                               cp=cp)

    logging.debug(f'inbound df len={len(df)} columns={df.columns}')
    logging.info(f'Got dataframe len={len(df)} Writing to {outfile}')
    df.to_csv(outfile, sep='\t')
    
    dir, base, ext = split_path(outfile)
    outfile = os.path.join( dir , f'{base}.parquet')
    logging.info(f'df len={len(df)} as parquet to {outfile}...')
    df.to_parquet(outfile)    
    
    # make_matrices
    process_make_matrices_pd(df,
                               exp_id = expid,  
                               outdir=outdir, 
                               cp=cp)
    

def filter_counts_fasta(cp, infile, 
                        outfile, 
                        datestr=None, 
                        min_count = 2 ):
    '''
    we want to go FASTA to FASTA, but filter by read counts.  
    '''
    sdf = read_fasta_to_df(infile)
    logging.debug(f'got {len(sdf)} sequences.')
    fdf = filter_counts_df(cp, sdf, min_count)
    logging.info(f'got {len(fdf)} filtered sequences. Writing to {outfile} FASTA.')
    write_fasta_from_df(fdf, outfile)
    logging.debug(f'Wrote {outfile}')
    
def calc_min_threshold(config, df, site='target-negative'):
    '''
    retrieve maximum umi_count value for provided site type. 
    returns max max value...

    '''
    logging.debug(f'getting UMI threshold type={site}')
    countlist = []
    min_threshold = 0
    df['umi_count'] = df['umi_count'].astype(int)
    tndf = df[ df['site'] == site]
    lablist = list(tndf['label'].dropna().unique())
    for label in lablist:
        ldf = tndf[tndf['label'] == label]
        if len(ldf) > 0:
            maxumi = ldf['umi_count'].max()
            logging.debug(f'max umi_count for label {label} = {maxumi}')
            countlist.append( maxumi )
    if len(countlist) > 0:
        min_threshold = max(countlist)
        logging.debug(f'calculated UMI threshold={min_threshold} type={site} ')    
    return min_threshold


def threshold_read_counts(config, df, threshold=1):
    '''
    
    '''
    logging.debug(f'threshold counts threshold={threshold}')
    threshold= int(threshold)   
    df = df[df['read_count'] >= threshold].copy()
    return df


def get_read_count_threshold(config, cdf, site=None):
    '''
    site = ['target-X','injection-X'] where X is '', lone, control, negative  
        Will use relevant threshold. If None, will use default threshold
    
    target_threshold=100
    target_ctrl_threshold=1000
    inj_threshold=2
    inj_ctrl_threshold=2
    
    '''
    logging.debug(f'getting info for site={site}')
    count_threshold=1
    if site is None:
        count_threshold = int(config.get('ssifasta', 'default_threshold'))
    else:
        count_threshold = int(config.get('ssifasta', f'{site}_threshold'))
    logging.debug(f'count threshold for {site} = {count_threshold}')
    return count_threshold

                 
def add_rowlist_column(rowlist, colval):
    """
    For use during dataframe construction. Adds col to list of rows with specified.
       
    """
    for row in rowlist:
        row.append(colval)
    return rowlist



