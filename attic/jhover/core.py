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

