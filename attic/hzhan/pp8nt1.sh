#!/bin/bash 
# run job in the current working directory where qsub is executed from
#$ -cwd
# specify that the job requires 4GB of memory
#$ -l m_mem_free=4G
#example preprocessing of seqeuncing libary ZL145

zcat Mseq204_YC_inj2_S1_R1_001.fastq.gz | awk "NR%4==2" - | cut -b 1-32 > Mseq204_YC_inj2_1_stripped.txt;
zcat Mseq204_YC_inj2_S1_R2_001.fastq.gz | awk "NR%4==2" - | cut -b 1-20 | paste -d '' Mseq204_YC_inj2_1_stripped.txt - > Mseq204_YC_inj2_PE.txt;


#split dataset according to inline indexes using fastx toolkit; this by default allows up to 1 missmatch. we could go higher if we want, though maybe not neccessary
mkdir barcodesplitter
cd barcodesplitter

nl ../Mseq204_YC_inj2_PE.txt | awk '{print ">" $1 "\n" $2}'| fastx_barcode_splitter.pl --bcfile ../barcode_v2.txt --prefix ../barcodesplitter/ --eol --exact

#from here on, do everything for every sample individually

BCidx=($(seq 240 1 270)) #the first number should be n-1
for i in {1..30}; do #this number should be exactly the same as the total RT primer number used in BCidx
	#filter out reads with Ns, cut off indexes and unique datafiles
	awk "NR%2==0" BC${BCidx[$i]} | grep -v N | cut -b 1-44 | sort | uniq -c | sort -nr > Mseq204_YC_inj2processedBC${BCidx[$i]}.txt
	#split output files into two files per index, one that is containing the read counts of each unique sequnce, the other the unique sequences themselves.
	awk '{print $1}' Mseq204_YC_inj2processedBC${BCidx[$i]}.txt > Mseq204_YC_inj2BC${BCidx[$i]}_counts.txt
	awk '{print $2}' Mseq204_YC_inj2processedBC${BCidx[$i]}.txt > Mseq204_YC_inj2_BC${BCidx[$i]}seq.txt
done


mkdir thresholds
cd thresholds



