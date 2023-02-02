#!/bin/bash -x
# run job in the current working directory where qsub is executed from
#$ -cwd
# specify that the job requires 4GB of memory
#$ -l m_mem_free=4G
#example preprocessing of seqeuncing libary ZL145

#unzip original datafiles                                                                                
gunzip -k *.gz

FQ1="M205_HZ_S1_R1_001"
FQ2="M205_HZ_S1_R2_001"
SAMPLE="M205_HZ_S1"
BCFILE="barcode_v2.txt"
  
#strip fastq files and clip sequences

awk "NR%4==2" ${FQ1}.fastq | cut -b 1-32 > ${FQ1}.stripped.txt # barcode + YY
rm ${FQ1}.fastq
awk "NR%4==2" ${FQ2}.fastq | cut -b 1-20 > ${FQ2}.stripped.txt #12nt tag + 8nt index
rm ${FQ2}.fastq


#make a new file that contains only one sequence per sequenced cluster
paste -d '\0' ${FQ1}.stripped.txt ${FQ2}.stripped.txt > ${SAMPLE}.paired.txt


#split dataset according to inline indexes using fastx toolkit; this by default allows up to 1 missmatch. we could go higher if we want, though maybe not neccessary
mkdir barcodesplitter
cd barcodesplitter

nl ../${SAMPLE}.paired.txt |awk '{print ">" $1 "\n" $2}'| fastx_barcode_splitter.pl --bcfile ../${BCFILE} --prefix ../barcodesplitter/ --eol --exact

#from here on, do everything for every sample individually

#BCidx=($(seq 0 1 81; seq 97 1 177; seq 193 1 233; seq 249 1 279)) #the first number should be n-1
BCidx=($(seq 0 1 26)) 
for i in {1..26}; do #this number should be exactly the same as the total RT primer number used in BCidx
	#filter out reads with Ns, cut off indexes and unique datafiles
	awk "NR%2==0" BC${BCidx[$i]} | grep -v N | cut -b 1-44 | sort | uniq -c | sort -nr > ${SAMPLE}.processed.BC${BCidx[$i]}.txt
	#split output files into two files per index, one that is containing the read counts of each unique sequence, the other the unique sequences themselves.
	awk '{print $1}' ${SAMPLE}.processed.BC${BCidx[$i]}.txt > ${SAMPLE}.BC${BCidx[$i]}.counts.txt
	awk '{print $2}' ${SAMPLE}.processed.BC${BCidx[$i]}.txt > ${SAMPLE}.BC${BCidx[$i]}.seq.txt
done


mkdir thresholds
cd thresholds
