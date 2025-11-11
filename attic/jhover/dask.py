def aggregate_reads_dd_clientXXX(df, 
                              column='sequence', 
                              outdir=None, 
                              min_reads=1, 
                              chunksize=50000000, 
                              dask_temp='./temp'):
    '''
    
    WARNING dask 
    
    ASSUMES INPUT IS DASK DATAFRAME
    retain other columns and keep first value

    Suggestions for partitions. 
    https://stackoverflow.com/questions/44657631/strategy-for-partitioning-dask-dataframes-efficiently

    Limiting resource usage (CPU, memory)
    https://stackoverflow.com/questions/59866151/limit-dask-cpu-and-memory-usage-single-node
    
    '''
    if dask_temp is not None:
        logging.info(f'setting dask temp to {os.path.expanduser(dask_temp)} ')
        dask.config.set(temporary_directory=os.path.abspath( os.path.expanduser(dask_temp)))
    else:
        logging.info(f'no dask_temp specified. letting dask use its default')


    client = Client( memory_limit='16GB', 
                     processes=True,
                     n_workers=8, 
                     threads_per_worker=2)
       
    before_len = len(df)
    logging.debug(f'collapsing with read counts for col={column} len={before_len}')
    sent = client.submit(sequence_value_counts, df)
    ndf = sent.result().compute()
    client.close()
    
    ndf.reset_index(inplace=True, drop=True)
    ndf.rename({'count':'read_count'}, inplace=True, axis=1)
    logging.info(f'computed counts. new DF len={len(ndf)}')
    
    # get back all other columns from original set. e.g. 'source'
    logging.info(f'merging to recover other columns from original DF')
    ndf = dd.from_pandas(ndf)
    #ndf = pd.merge(ndf, seqdf.drop_duplicates(subset=column,keep='first'),on=column, how='left')  
    result = ndf.merge(seqdf.drop_duplicates(subset=column,keep='first'), on=column, how='left')  
    ndf = result.compute()
    logging.info(f'got merged DF=\n{ndf}')
  
    if min_reads > 1:
        logging.info(f'Dropping reads with less than {min_reads} read_count.')
        logging.debug(f'Length before read_count threshold={len(ndf)}')
        ndf = ndf[ndf['read_count'] >= min_reads]
        ndf.reset_index(inplace=True, drop=True)
        logging.info(f'Length after read_count threshold={len(ndf)}')    
    else:
        logging.info(f'min_reads = {min_reads} skipping initial read count thresholding.')  
    logging.info(f'final output DF len={len(ndf)}')    
    return ndf

