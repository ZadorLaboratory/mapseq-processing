#  
#  Wrapper for fast demultiplexing utility ultraplex. 
#
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8287537/
#
#    conda install ultraplex 
#
#    ultraplex 
#           --barcodes           barcodes.csv
#           --inputfastq         M205_HZ_S1_R1_001.fastq.gz
#           --input_2            M205_HZ_S1_R2_001.fastq.gz
#           --adapter            GTAC
#           --adapter2           TGTACAGCTAGCGGTG 
#        
#
#  barcodes file needs 5'
#
#
#
# read1
#   30ntRAND NN GTAC
# read2
#   12ntUMI  8ntSSI 