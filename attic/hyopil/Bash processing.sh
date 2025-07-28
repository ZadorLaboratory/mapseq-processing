#!/bin/bash 
#raw files are already demultiplexed by the core according to index i5 i7 pairs
#Give name for your experiment
Prefix='dat8'
#prepare CSI txt file
#Check whether the fastq.gz file name fits the below script and check CSI (or SSI)for i in {1..96};do    zcat raw/${Prefix}_${i}_*_L00[0-9]_R1_001.fastq.gz | awk "NR%4==2" | cut -b 1-32 > ${Prefix}_${i}_R1_stripped.txt    zcat raw/${Prefix}_${i}_*_L00[0-9]_R2_001.fastq.gz | awk "NR%4==2" | cut -b 1-20 | paste -d '' ${Prefix}_${i}_R1_stripped.txt - > ${Prefix}_BC${i}_PE.txtdone
# if you want check CSI, run this.#If everything goes well, most of samples should have correct CSI
for i in {1..96};do    cut -b 45-52 ${Prefix}_BC${i}_PE.txt | sort | uniq -c | sort -nr > ${Prefix}_BC${i}_CSIcheck.txt    echo ${i}
    head -5 ${Prefix}_BC${i}_CSIcheck.txt
done
#filtering reads only have correct CSI
for i in {1..96};do	awk "NR==${i}" CSI.txt > CSI	C=$(cut -b 1-8 CSI)	grep "$C"$ ${Prefix}_BC${i}_PE.txt | grep -v N | cut -b 1-44 | sort | uniq -c | sort -nr > ${Prefix}_BC${i}_BU.txt	awk '{print $1}' ${Prefix}_BC${i}_BU.txt > ${Prefix}_BC${i}_BU_counts.txt	awk '{print $2}' ${Prefix}_BC${i}_BU.txt > ${Prefix}_BC${i}_BU_seq.txt
done
#Determine threshold using rankplot in matlabmkdir thresholdscd thresholds
#Set the first position as 0 and give your thresholds from the second position, for example, threshold=(0 3 3 3 3 3 3 3 3)#do a quick collpase around the UMI tags
for i in {1..96}; do
    j=${threshold[$i]} 
    echo ${j}
    if [ "$j" != "1" ];then
	awk '$1 < '$j' {print NR}' ../${Prefix}_BC${i}_BU_counts.txt | head -1 >t
	thresh=$(cat t)
	echo $thresh
	head ../${Prefix}_BC${i}_BU_seq.txt -n $thresh | cut -b 1-32 | sort | uniq -c | sort -nr > ${Prefix}_BC${i}_quickout.txt
   else	cut -b 1-32 ../${Prefix}_BC${i}_BU_seq.txt | sort | uniq -c | sort -nr > ${Prefix}_BC${i}_quickout.txt
   fi
done
#ok, thats a lot of preprocessing done. Now we have to error correct these sequences. The all against all mapping of barcodes in using bowtie is done in the command line. after this we move into MATLAB.
# this script suppose that you used L4 spikein, if you didn't you should change the sequence 'CGTC' in the below script as your spikein 

mkdir indexes
for i in {1..96};do	echo ${i}	in=${Prefix}_BC${i}_quickout.txt	grep -v 'CGTCAGTC$' ${in} > ${Prefix}_BC${i}_quickprocessed.txt	awk '{print $1}' ${Prefix}_BC${i}_quickprocessed.txt > ${Prefix}_${i}_counts.txt	awk '{print $2}' ${Prefix}_BC${i}_quickprocessed.txt > ${Prefix}_${i}_seq.txt	nl ${Prefix}_${i}_seq.txt | awk '{print ">" $1 "\n" $2}' > ${Prefix}_${i}fasta2u.txt	bowtie-build -q ${Prefix}_${i}fasta2u.txt indexes/${i}fasta2u	bowtie -v 3 -p 10 -f --best -a indexes/${i}fasta2u ${Prefix}_${i}fasta2u.txt bowtiealignment${i}_2u.txt	awk '{print $1}' bowtiealignment${i}_2u.txt > bowtie${i}_2u_1.txt
	awk '{print $3}' bowtiealignment${i}_2u.txt > bowtie${i}_2u_3.txt
done
#  you want to deal with spike insfor i in {1..96};do	echo ${i}	in=${Prefix}_BC${i}_quickout.txt	grep 'CGTCAGTC$' ${in} > ${Prefix}_spikesBC${i}_quickprocessed.txt	awk '{print $1}' ${Prefix}_spikesBC${i}_quickprocessed.txt > ${Prefix}_spikes${i}_counts.txt	awk '{print $2}' ${Prefix}_spikesBC${i}_quickprocessed.txt > ${Prefix}_spikes${i}_seq.txt	nl ${Prefix}_spikes${i}_seq.txt | awk '{print ">" $1 "\n" $2}' > ${Prefix}_spikes${i}fasta2u.txt	bowtie-build -q ${Prefix}_spikes${i}fasta2u.txt indexes/spikes${i}fasta2u	bowtie -v 3 -p 10 -f --best -a indexes/spikes${i}fasta2u ${Prefix}_spikes${i}fasta2u.txt bowtiealignmentspikes${i}_2u.txt	awk '{print $1}' bowtiealignmentspikes${i}_2u.txt > bowtiespikes${i}_2u_1.txt	awk '{print $3}' bowtiealignmentspikes${i}_2u.txt > bowtiespikes${i}_2u_3.txt
done