





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