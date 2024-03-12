







def process_mapseq_dir(exp_id, loglevel, force):
    '''
    E.g. 
    
    process_fastq.py -v -b barcode_v2.txt -s M205_sampleinfo.xlsx -O fastq.20231204.1213.out fastq/M205_HZ_S1_R1_002.fastq.gz fastq/M205_HZ_S1_R2_002.fastq.gz
    process_ssifasta.py -v -a bowtie2 -s M205_sampleinfo.xlsx -O ssifasta.20231204.1213.out -o M205.all.20231204.1213.tsv fastq.20231204.1213.out/BC*.fasta
    process_merged.py -v -s M205_sampleinfo.xlsx -e M205_20231204.1213 -l label -O merged.20231204.1213.out M205.all.20231204.1213.tsv 
    make_heatmaps.py -v -s M205_sampleinfo.xlsx -e M205_20231204.1213 -O merged.20231204.1213.out   merged.20231204.1213.out/*.nbcm.tsv  
    
    '''   
    d = os.path.abspath(exp_id)
    if not os.path.exists(d):
        sys.exit(f'Experiment directory {d} does not exist.')

    logging.info(f'processing experiment dir: {d}')
    
    config = get_default_config()
    expconfig = f'{d}/mapseq.conf'

    if os.path.exists(expconfig):
         config.read(expconfig)
         logging.debug(f'read {expconfig}')
    
         
    try:
       samplefile = f'{d}/{exp_id}_sampleinfo.jrh.xlsx'
       sampdf = load_sample_info(config, samplefile)
       rtlist = get_rtlist(sampdf)
       
       outdir = f'{d}/fastq.out'
       bcfile = f'{d}/barcode_v2.{exp_id}.txt'       
       bclist = load_barcodes(config, bcfile, labels=rtlist, outdir=outdir)
       readfilelist = package_pairfiles( glob.glob(f'{d}/fastq/*.fastq*'))

       logging.info(f'running process_fastq_pairs. readfilelist={readfilelist} outdir={outdir}')
       process_fastq_pairs(config, readfilelist, bclist, outdir, force=False)
       
         
       #process_ssifasta(config, infile, outdir=None, site=None)
       #process_merged(config, filelist, outdir=None, expid=None, recursion=100000, combined_pdf=True)
       #process_qc(config, exp_dir)

    except Exception as ex:
        logging.error(f'error while handling {d} ')
        logging.warning(traceback.format_exc(None))
        
        
