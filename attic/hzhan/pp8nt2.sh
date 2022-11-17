#!/bin/bash 
# run job in the current working directory where qsub is executed from
#$ -cwd
# specify that the job requires 4GB of memory
#$ -l m_mem_free=4G
#example preprocessing of seqeuncing libary ZL145


FQ1="M205_HZ_S1_R1_001"
FQ2="M205_HZ_S1_R2_001"
SAMPLE="M205_HZ_S1"

cd barcodesplitter
cd thresholds
  
#pick thresholds from matlab plots, based on a steep drop of the 32+12 read counts. avoids too many PCR and sequnecning errors in the data themselves. note, bash seems to start counting at 0, so put a random number at the first position of this array, to get the order right.

BCidx=($(seq 240 1 270))
threshold=(0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2) #you can put more than the actual RT number

#just do a quick and dirty collpase around the UMI tags
for i in {1..30}; do
    j=${threshold[$i]} 
    echo $j
    if [ "$j" != "1" ];then
		awk '$1< '$j' {print NR}' ../${SAMPLE}2BC${BCidx[$i]}_counts.txt | head -n 1 >t
		thresh=$(cat t)
		echo $thresh
		head ../${SAMPLE}2_BC${BCidx[$i]}seq.txt -n $thresh | cut -b 1-32 | sort | uniq -c | sort -nr > ${SAMPLE}2${i}quickout.txt

   	else
		grep -nr ^${threshold[$i]}$  ../${SAMPLE}2BC${BCidx[$i]}_counts.txt -m 1 | cut -f1 -d":" > t
		thresh=$(cat t)
		echo $thresh

		head ../${SAMPLE}2_BC${BCidx[$i]}seq.txt -n $thresh | cut -b 1-32 | sort | uniq -c | sort -nr > ${SAMPLE}2${i}quickout.txt
   fi
done


#ok, thats a lot of preprocessing done. Now we have to error correct these sequneces. The all against all mapping of barcodes in using bowtie is done in the command line. after this we move into MATLAB.

mkdir indexes

for i in {1..30}
do
	echo $i
	in=${SAMPLE}2${i}quickout.txt
	#split off real barcodes from spike-ins

	grep -v 'CGTCAGTC$' $in | grep '[TC][TC]$' > ${SAMPLE}2BC${i}_quickprocessed.txt
	awk '{print $1}' ${SAMPLE}2BC${i}_quickprocessed.txt > ${SAMPLE}2${i}_counts.txt
	awk '{print $2}' ${SAMPLE}2BC${i}_quickprocessed.txt > ${SAMPLE}2${i}_seq.txt


	nl ${SAMPLE}2${i}_seq.txt | awk '{print ">" $1 "\n" $2}' > ${SAMPLE}2_BC${i}fasta2u.txt; 
	bowtie-build -q ${SAMPLE}2_BC${i}fasta2u.txt indexes/BC${i}fasta2u; 
	bowtie -v 3 -p 10 -f --best -a indexes/BC${i}fasta2u ${SAMPLE}2_BC${i}fasta2u.txt bowtiealignment${i}_2u.txt
	awk '{print $1}' bowtiealignment${i}_2u.txt > bowtie${i}_2u_1.txt;awk '{print $3}' bowtiealignment${i}_2u.txt > bowtie${i}_2u_3.txt
done

# now deal with spike ins

for i in {1..30}; do 
	echo $i; in=${SAMPLE}2${i}quickout.txt 
	grep 'CGTCAGTC$' $in > ${SAMPLE}2spikes${i}_quickprocessed.txt; 
	awk '{print $1}' ${SAMPLE}2spikes${i}_quickprocessed.txt > ${SAMPLE}2spikes${i}_counts.txt; 
	awk '{print $2}' ${SAMPLE}2spikes${i}_quickprocessed.txt > ${SAMPLE}2spikes${i}_seq.txt;  
	nl ${SAMPLE}2spikes${i}_seq.txt | awk '{print ">" $1 "\n" $2}' > ${SAMPLE}2_spikes${i}fasta2u.txt; 
	bowtie-build -q ${SAMPLE}2_spikes${i}fasta2u.txt indexes/spikes${i}fasta2u; bowtie -v 3 -p 10 -f --best -a indexes/spikes${i}fasta2u ${SAMPLE}2_spikes${i}fasta2u.txt bowtiealignmentspikes${i}_2u.txt; 
	awk '{print $1}' bowtiealignmentspikes${i}_2u.txt > bowtiespikes${i}_2u_1.txt;
	awk '{print $3}' bowtiealignmentspikes${i}_2u.txt > bowtiespikes${i}_2u_3.txt; 
done

# L1 barcodes
for i in {1..30}
do
	echo $i
	in=${SAMPLE}2${i}quickout.txt
	#split off real barcodes from spike-ins

	grep -v 'CGTCAGTC$' $in | grep '[AG][AG]$' > ${SAMPLE}2BC${i}_quickprocessedL1.txt
	awk '{print $1}' ${SAMPLE}2BC${i}_quickprocessedL1.txt > ${SAMPLE}2${i}_countsL1.txt
	awk '{print $2}' ${SAMPLE}2BC${i}_quickprocessedL1.txt > ${SAMPLE}2${i}_seqL1.txt


	nl ${SAMPLE}2${i}_seqL1.txt | awk '{print ">" $1 "\n" $2}' > ${SAMPLE}2_BC${i}fasta2uL1.txt;
	bowtie-build -q ${SAMPLE}2_BC${i}fasta2uL1.txt indexes/BC${i}fasta2uL1;
	bowtie -v 3 -p 10 -f --best -a indexes/BC${i}fasta2uL1 ${SAMPLE}2_BC${i}fasta2uL1.txt bowtiealignment${i}_2uL1.txt
	awk '{print $1}' bowtiealignment${i}_2uL1.txt > bowtie${i}_2u_1L1.txt;awk '{print $3}' bowtiealignment${i}_2uL1.txt > bowtie${i}_2u_3L1.txt
done

