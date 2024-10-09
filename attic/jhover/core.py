def process_ssifasta(config, infile, outdir=None, site=None, datestr=None):
    '''
    by default, outdir will be same dir as infile
    assumes infile fasta has already been trimmed to remove SSI
    
    site = ['target-control','injection-control','target','target-negative',target-lone']   
    Will use relevant threshold. If None, will use default threshold
    
    '''
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
    unmatched.to_csv(os.path.join(outdir , f'{base}.unmatched.counts.tsv'), sep='\t')

    # remove homopolymers in real sequences.
    #
    #  base repeats already assumed removed in newer processing. 
    #
    # max_homopolymer_run=int(config.get('ssifasta', 'max_homopolymer_run')) 
    # realdf = remove_base_repeats(realdf, col='sequence', n=max_homopolymer_run)
     
    # align and collapse all.         

    # add labels for merging...
    realdf['type'] = 'real'
    spikedf['type'] = 'spike'
    lonedf['type'] = 'lone'
    unmatched['type'] = 'unmatched'

    outdf = merge_dfs([ realdf, spikedf, lonedf, unmatched ])
    outdf['label'] = base
    outdf.sort_values(by = ['type', 'umi_count'], ascending = [True, False], inplace=True)
    outdf.reset_index(drop=True, inplace=True)
    return outdf



def collapse_by_components(fulldf, uniqdf, components, outdir=None):
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
    seqmapdict = build_seqmapdict(uniqdf, components)
      
    # Make new full df:
    logging.info('seqmapdict built. Applying.')
    if outdir is not None:
        outfile = f'{outdir}/seqmapdict.json'
        logging.debug(f'writing seqmapdict len={len(seqmapdict)} tp {outfile}')
        with open(outfile, 'w') as fp:
            json.dump(seqmapdict, fp, indent=4)
    else:
        logging.debug(f'no outdir given.')
    logging.info(f'applying seqmapdict...')
    fulldf['newsequence'] = fulldf.apply(apply_setcompseq, axis=1, seqmapdict=seqmapdict)
    logging.info(f'New collapsed df = \n{fulldf}')
    log_objectinfo(fulldf, 'fulldf')
    return fulldf


def align_collapse_df(df, seq_length=None, max_mismatch=None, outdir=None, datestr=None, cp=None):
    '''
    Assumes dataframe with sequence and read_count columns
    
    '''
    # housekeeping...
    aligner = cp.get('collapse','tool')   
    if seq_length is None:
        seq_length = int(cp.get('collapse', 'seq_length'))    
    if max_mismatch is None:
        max_mismatch = int(cp.get('collapse', 'max_mismatch'))
    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")
    if outdir is None:
        outdir = './'
    outdir = os.path.abspath(outdir)    
    os.makedirs(outdir, exist_ok=True)
    logging.debug(f'collapse: aligner={aligner} seq_length={seq_length} max_mismatch={max_mismatch} outdir={outdir}')
    
    sh = StatsHandler(cp, outdir=outdir, datestr=datestr)      

    sh.add_value('/collapse','n_full_sequences', len(df) )


    # rename column, and split by sequence length
    df.rename( {'sequence': 'full_read'}, inplace=True, axis=1)
    df['sequence'] = df['full_read'].str.slice(0,seq_length)
    df['tail'] = df['full_read'].str.slice(start=seq_length)

    # get reduced dataframe of unique head sequences
    logging.info('Getting unique DF...')    
    udf = pd.DataFrame(df['sequence'].unique(), columns=['sequence'])
    #udf = df['sequence'].value_counts() 
    sh.add_value('/collapse','n_unique_sequences', len(udf) )    

    of = os.path.join( outdir , 'unique_sequences.tsv')
    logging.info(f'Writing unique DF to {of}')
    udf.to_csv(of, sep='\t') 
    
    of = os.path.join( outdir , 'unique_sequences.fasta')
    logging.info(f'Writing uniques as FASTA to {of}')
    seqfasta = write_fasta_from_df(udf, outfile=of, sequence=['sequence'])    

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

    components = remove_singletons(components)
    logging.debug(f'multi-element components len={len(components)}')
    sh.add_value('/collapse','n_multi_components', len(components) )
    of = os.path.join( outdir , f'multi_components.txt')
    writelist(of, components)        

    logging.info(f'Collapsing {len(components)} components...')
    newdf = collapse_by_components(df, udf, components, outdir=outdir)
    # newdf has sequence and newsequence columns, rename to orig_seq and sequence
    newdf.rename( {'sequence':'orig_seq', 'newsequence':'new_seq'}, inplace=True, axis=1)    
    logging.info(f'Got collapsed DF. len={len(newdf)}')
        
    joindf = pd.DataFrame( newdf['new_seq'] + newdf['tail'], columns=['sequence'])
    of = os.path.join( outdir , f'collapsed.fasta')
    logging.info(f'Writing collapsed fasta to {of}')
    write_fasta_from_df(joindf, of)        
    return newdf



def process_fasta(config, infile, bclist, outdir, force=False, datestr=None):
    '''
    Take paired FASTA file as input rather than raw FASTQ
    Use combined barcode sorter structure to minimize character comparisons. 

    [splitfasta]
    ssifile = ~/git/mapseq-processing/etc/barcode_v2.txt
    seqhandled_interval = 1000000
    matched_interval = 1000000
    unmatched_interval = 1000000

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
    
    seqhandled_interval = int(config.get('splitfasta','seqhandled_interval')) 
    matched_interval = int(config.get('splitfasta','matched_interval'))
    unmatched_interval = int(config.get('splitfasta','unmatched_interval'))
    
    logging.info(f'performing split. outdir={outdir} output_exists={output_exists} force={force}')
    
    
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
    
   

def process_ssifasta_files(config, sampleinfo, infilelist, numthreads=1, outdir=None):
    '''
    Process each infile in separate process.
    Produce {outdir}/<BASE>.all.tsv 
    These are then combined to product outfile.   
    
    '''
    logging.debug(f'called: sampleinfo={sampleinfo} outdir={outdir}')
    if outdir is None:
        afile = infilelist[0]
        filepath = os.path.abspath(afile)    
        dirname = os.path.dirname(filepath)
        outdir = dirname
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
        logging.debug(f'made outdir={outdir}')

    ncpus, threads = calc_thread_count(numthreads)   
    logging.info(f'using {threads} of {ncpus} CPUs in parallel...')
    proglog = '-v'
    if logging.getLogger().level == logging.INFO:
        proglog = '-v'
    elif logging.getLogger().level == logging.DEBUG:
        proglog = '-d'
    
   
    cfilename =  f'{outdir}/process_ssifasta.config.txt'
    configfile = write_config(config, cfilename, timestamp=True)
    
    prog = os.path.expanduser('~/git/mapseq-processing/scripts/process_ssifasta_single.py')

    jstack = JobStack()
    outfilelist = []
    for infile in infilelist:
        base = get_mainbase(infile)
        logfile = f'{outdir}/{base}.log'
        outfile = f'{outdir}/{base}.all.tsv'
        outfilelist.append(outfile)
        cmd = [ prog, 
               proglog,
               '-c', configfile , 
               '-L', logfile ,
               '-s', sampleinfo ,  
               '-O' , outdir ,
               '-o' , outfile,
               infile
               ]
        jstack.addjob(cmd)
    jset = JobSet(max_processes = threads, jobstack = jstack)
    jset.runjobs()

    logging.info(f'finished jobs. merging output files:{outfilelist} ')    
    outdf = merge_tsvs(outfilelist)
    logging.info(f'made final DF len={len(outdf)}')
    return outdf





def process_merged(config, infile, outdir=None, expid=None, recursion=200000, label_column='region' ):
    '''
     takes in combined 'all' TSVs. columns=(sequence, counts, type, label, brain, site) 
     outputs brain-specific SSI x target matrix dataframes, with counts normalized to spikeins by target.  
     writes all output to outdir (or current dir). 
     
    '''
    #sh = get_default_stats()
    logging.debug(f'infile={infile}')
    alldf = load_df(infile)
    logging.debug(f'alldf=\n{alldf}')
    if outdir is None:
        outdir = './'

    require_injection = config.getboolean('merged','require_injection')
    min_injection = int(config.get('merged','min_injection'))
    min_target = int(config.get('merged','min_target'))   
    use_target_negative=config.getboolean('merged','use_target_negative')
    use_target_water_control=config.getboolean('merged','use_target_water_control')
    
    clustermap_scale = config.get('analysis','clustermap_scale')
      
    if expid is None:
        expid = 'M000'
    
    logging.debug(f'running exp={expid} min_injection={min_injection} min_target={min_target} use_target_negative={use_target_negative} use_target_water_control={use_target_water_control}')

    alldf['brain'] = alldf['brain'].astype('string')
    bidlist = list(alldf['brain'].dropna().unique())
    bidlist = [ x for x in bidlist if len(x) > 0 ]
    bidlist.sort()
    logging.debug(f'handling brain list: {bidlist}')
    
    for brain_id in bidlist:
        valid = True
        logging.debug(f'handling brain_id={brain_id}')
        bdf = alldf[alldf['brain'] == brain_id]
                    
        # handle target areas...
        tdf = bdf[bdf['site'].str.startswith('target')]
        rtdf = tdf[tdf['type'] == 'real'] 

        # threshold by min_target or threshold by target-negative
        # if use_target_negative is true, but no target negative site 
        # defined, use min_target and throw warning. 
        if use_target_negative:
            logging.info(f'use_target_negative is {use_target_negative}')
            max_negative = calc_min_threshold(config, bdf, 'target-negative')
            logging.debug(f'target-negative UMI count = {max_negative}')

        if use_target_water_control:
            logging.info(f'use_target_water_control is {use_target_water_control}')
            max_water_control = calc_min_threshold(config, bdf, 'target-water-control')
            logging.debug(f'target_water_control UMI count = {max_water_control}')
        
        min_target = max([ min_target, max_negative, max_water_control ])
        logging.debug(f'min_target UMI count after all constraints = {min_target}')   

        # min_target is now either calculated from target-negative, or from config. 
        if min_target > 1:
            before = len(rtdf)
            frtdf = filter_all_lt(rtdf, 'sequence', 'umi_count', min_target)            
            if not len(frtdf) > 0:
                valid = False
                logging.warning(f'No VBCs passed min_target filtering! Skip brain.')
            logging.debug(f'filtering by min_target={min_target} before={before} after={len(rtdf)}')
        else:
            logging.debug(f'min_target={min_target} no filtering.')
            frtdf = rtdf

        if require_injection:
            # extract and filter injection areas.
            logging.debug(f'require_injection={require_injection} min_injection={min_injection}') 
            idf = bdf[bdf['site'].str.startswith('injection')]
            ridf = idf[idf['type'] == 'real']  
            if len(ridf) == 0:
                logging.warning('require_injection=True but no real VBCs from any injection site.')
            logging.debug(f'{len(frtdf)} real target VBCs before filtering.')      
            frtdf = filter_non_injection(frtdf, ridf, min_injection=min_injection)
            logging.debug(f'{len(frtdf)} real target VBCs after injection filtering.')
            if not len(frtdf) > 0:
                logging.warning(f'No VBCs passed injection filtering! Skip brain.')
                valid = False
        else:
            logging.debug(f'require_injection={require_injection} proceeding...')
               
        # make matrices if brain data is valid... 
        if valid:       
            
            # raw reals
            rbcmdf = rtdf.pivot(index='sequence', columns=label_column, values='umi_count')
            scol = natsorted(list(rbcmdf.columns))
            rbcmdf = rbcmdf[scol]
            rbcmdf.fillna(value=0, inplace=True)
            logging.debug(f'brain={brain_id} raw real barcode matrix len={len(rbcmdf)}')            
            
            # filtered reals
            fbcmdf = frtdf.pivot(index='sequence', columns=label_column, values='umi_count')
            scol = natsorted(list(fbcmdf.columns))
            fbcmdf = fbcmdf[scol]
            fbcmdf.fillna(value=0, inplace=True)
            logging.debug(f'brain={brain_id} real barcode matrix len={len(fbcmdf)}')

            # spikes
            sdf = tdf[tdf['type'] == 'spike']
            sbcmdf = sdf.pivot(index='sequence', columns=label_column, values='umi_count')
            spcol = natsorted(list(sbcmdf.columns))
            sbcmdf = sbcmdf[spcol]
            sbcmdf.fillna(value=0, inplace=True)    
            logging.debug(f'brain={brain_id} spike barcode matrix len={len(sbcmdf)}')
    
            (fbcmdf, sbcmdf) = sync_columns(fbcmdf, sbcmdf)
            
            
            nbcmdf = normalize_weight(fbcmdf, sbcmdf)
            logging.debug(f'nbcmdf.describe()=\n{nbcmdf.describe()}')
            
            scbcmdf = normalize_scale(nbcmdf, logscale=clustermap_scale)
            scbcmdf.fillna(value=0, inplace=True)
            logging.debug(f'scbcmdf.describe()=\n{scbcmdf.describe()}')
                        
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
            
        logging.info(f'done with brain={brain_id}')
        
        


def align_collapse_fasta(config, infile, seq_length=None, max_mismatch=None, outdir=None, datestr=None):
    '''
    Algorithm:
    
    take input FASTA to FULL dataframe. ( already removed Ns and homopolymer runs ) 
    split first seq_length part (VBC) and remainder to 2 columns. 
    
    get UNIQUE DF dataframe on (VBC) sequence      
    save VBCs to fasta
    align all to all bowtie
    make bowtie DF
    drop mismatch > max_mismatch
        for each component, create map from all variants back to a single virtual parent sequence
        (this is not necessarily the true biological parent sequence, but it doesn't matter. they just need to all
        be the same).
    in FULL DF, set all variant sequences to virtual parent sequence
    re-assemble VBC + remainder sequence
    output to adjusted FULL FASTA file ready for input to SSI splitting.  
   
    can use iloc and range of indices to set value for all elements of components:
    
    testdf = pd.DataFrame([[0, 2, 3], [0, 4, 1], [10, 20, 30]],columns=['A', 'B', 'C'])
    testdf
            A   B   C
        0   0   2   3
        1   0   4   1
        2  10  20  30
    testdf.iloc[[0,2],[0]] = 15
            A   B   C
        0  15   2   3
        1   0   4   1
        2  15  20  30

https://stackoverflow.com/questions/46204521/pandas-get-unique-values-from-column-along-with-lists-of-row-indices-where-the
    
    '''
    aligner = config.get('collapse','tool')
    
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    head = base.split('.')[0]
    
    if seq_length is None:
        seq_length = int(config.get('collapse', 'seq_length'))
    
    if max_mismatch is None:
        max_mismatch = int(config.get('collapse', 'max_mismatch'))
    
    if outdir is None:
        outdir = dirname    
    os.makedirs(outdir, exist_ok=True)

    if datestr is None:
        datestr = dt.datetime.now().strftime("%Y%m%d%H%M")

    sh = StatsHandler(config, outdir=outdir, datestr=datestr)
    #  sh.add_value('/collapse','n_components', len(components) )
    
    # need to strip to seq_length
    logging.info(f'Reading {infile} to df...')
    fdf = read_fasta_to_df(infile, seq_length)
    logging.debug(f'fdf=\n{fdf}')
    of = os.path.join( outdir , f'{base}.fulldf.tsv')
    logging.info(f'Writing full DF to {of}')
    fdf.to_csv(of, sep='\t')
    sh.add_value('/collapse','n_full_sequences', len(fdf) )

    # get reduced dataframe of unique sequences
    logging.info('Getting unique DF...')    
    udf = pd.DataFrame(fdf['sequence'].unique(), columns=['sequence'])
    
    sh.add_value('/collapse','n_unique_sequences', len(udf) )
    
    logging.debug(f'udf = \n{udf}')
    of = os.path.join( outdir , f'{base}.udf.tsv')
    logging.info(f'Writing unique DF to {of}')
    udf.to_csv(of, sep='\t') 
    
    of = os.path.join( outdir , f'{base}.udf.fasta')
    logging.info(f'Writing uniques as FASTA to {of}')
    seqfasta = write_fasta_from_df(udf, outfile=of)
    
    # run allXall bowtiex
    of = os.path.join( outdir , f'{base}.bt2.sam')
    logging.info(f'Running {aligner} on {seqfasta} file to {of}')
    # switch to generic bowtie later... JRH
    #afile = run_bowtie(config, seqfasta, seqfasta, of, tool=aligner)
    btfile = run_bowtie(config, seqfasta, of, tool=aligner)
    logging.info(f'Bowtie done. Produced output {btfile}. Creating btdf dataframe...')
    btdf = make_bowtie_df(btfile, max_mismatch=max_mismatch, ignore_self=True)
    of = os.path.join( outdir , f'{base}.btdf.tsv')
    btdf.to_csv(of, sep='\t') 
    sh.add_value('/collapse','n_bowtie_entries', len(btdf) )
    
    # perform collapse...      
    logging.info('Calculating Hamming components...')
    edgelist = edges_from_btdf(btdf)
    btdf = None  # help memory usage
    sh.add_value('/collapse','n_edges', len(edgelist) )
    
    logging.debug(f'edgelist len={len(edgelist)}')
    components = get_components(edgelist)
    logging.debug(f'all components len={len(components)}')
    sh.add_value('/collapse','n_components', len(components) )
    edgelist = None  # help memory usage

    # components is list of lists.
    data = [ len(c) for c in components]
    data.sort(reverse=True)
    ccount = pd.Series(data)
    of = os.path.join( outdir , f'{base}.comp_count.tsv')
    ccount.to_csv(of, sep='\t') 
 
    components = remove_singletons(components)
    logging.debug(f'multi-element components len={len(components)}')
    sh.add_value('/collapse','n_multi_components', len(components) )
    of = os.path.join( outdir , f'{base}.multi.components.txt')
    writelist(of, components)
    
    logging.info(f'Collapsing {len(components)} components...')
    newdf = collapse_by_components(fdf, udf, components, outdir=outdir)
    # newdf has sequence and newsequence columns, rename to orig_seq and sequence
    newdf.rename( {'sequence':'orig_seq', 'newsequence':'sequence'}, inplace=True, axis=1)
    
    del fdf 
    del udf
    del components
        
    of = os.path.join( outdir , f'{base}.collapsed.tsv')
    logging.info(f'Got collapsed DF. Writing to {of}')
    newdf.to_csv(of, sep='\t')     
    logging.info('Done. Calculating fasta.')
    of = os.path.join( outdir , f'{base}.collapsed.fasta')
    cdf = pd.DataFrame( newdf['sequence'] + newdf['tail'], columns=['sequence'])
    logging.info(f'Writing fasta to {of}')
    write_fasta_from_df(cdf, of)
    logging.info(f'Wrote re-joined sequences to {of}')



def apply_setcompseq(row, seqmapdict ):
    ''' 
      Set Component Sequence
      Applies mapping from old to new sequence from dict. 
      
      seqmapdict = { oldseq : newseq }
      
      Usage example:
            df['sequence'] = df.apply(apply_setcompseq, axis=1, seqmapdict=seqmapdict)
      
    '''
    try:
        #logging.debug(f"getting mapping for {row['sequence']}")
        a = seqmapdict[row['sequence']]
    except KeyError:
        # A KeyError means we don't want to remap, just use original value
        a = row['sequence']
    return a

