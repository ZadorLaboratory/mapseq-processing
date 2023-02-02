#!/bin/bash -x
#example preprocessing of sequencing libary M205

FQ1="M205_HZ_S1_R1_001"
FQ2="M205_HZ_S1_R2_001"
SAMPLE="M205_HZ_S1"

cd barcodesplitter/thresholds
 
# pick thresholds from matlab plots, based on a steep drop of the 32+12 read counts. 
# avoids too many PCR and sequencing errors in the data themselves. Note, bash seems 
# to start counting at 0, so put a random number at the first position of this array, 
# to get the order right.

BCidx=($(seq 0 1 26))
threshold=(2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ) #you can put more than the actual RT number

#just do a quick and dirty collapse around the UMI tags
for i in {1..28}; do
    j=${threshold[$i]} 
    #echo " threshval = $j"
    if [ "$j" != "1" ]; then
    	echo "awk -v j=$j '$1<j {print NR}' ../${SAMPLE}.BC${BCidx[$i]}.counts.txt | head -n 1 > t$i"
		awk -v j=$j '$1<j {print NR}' ../${SAMPLE}.BC${BCidx[$i]}.counts.txt | head -n 1 > t$i
		thresh=$(cat t$i)
		echo "thresh = ${thresh}"
		echo "head  -n $thresh ../${SAMPLE}.BC${BCidx[$i]}.seq.txt | cut -b 1-32 | sort | uniq -c | sort -nr > ${SAMPLE}.${i}.quickout.txt"
		head  -n $thresh ../${SAMPLE}.BC${BCidx[$i]}.seq.txt | cut -b 1-32 | sort | uniq -c | sort -nr > ${SAMPLE}.${i}.quickout.txt
    else
		grep -nr ^${threshold[$i]}$  ../${SAMPLE}.BC${BCidx[$i]}.counts.txt -m 1 | cut -f1 -d":" > t$i
		thresh=$(cat t$i)
		echo "thresh = ${thresh}"
		head ../${SAMPLE}.BC${BCidx[$i]}.seq.txt -n $thresh | cut -b 1-32 | sort | uniq -c | sort -nr > ${SAMPLE}.${i}.quickout.txt
    fi
done

#ok, thats a lot of preprocessing done. Now we have to error correct these sequences. 
# The all against all mapping of barcodes using bowtie is done on the command line. 
# after this we move into MATLAB.

mkdir indexes

for i in {1..26}
do
	echo $i
	in=${SAMPLE}.${i}.quickout.txt
	echo "handling real ${in} ...}"
	#split off real barcodes from spike-ins

	grep -v 'CGTCAGTC$' $in | grep '[TC][TC]$' > ${SAMPLE}.BC${i}.quickprocessed.txt
	awk '{print $1}' ${SAMPLE}.BC${i}.quickprocessed.txt > ${SAMPLE}.${i}.counts.txt
	awk '{print $2}' ${SAMPLE}.BC${i}.quickprocessed.txt > ${SAMPLE}.${i}.seq.txt

	nl ${SAMPLE}.${i}.seq.txt | awk '{print ">" $1 "\n" $2}' > ${SAMPLE}.BC${i}.fasta2u.txt 
	bowtie-build -q ${SAMPLE}.BC${i}.fasta2u.txt indexes/BC${i}fasta2u 
	bowtie -v 3 -p 10 -f --best -a indexes/BC${i}fasta2u ${SAMPLE}.BC${i}.fasta2u.txt bowtiealign${i}.2u.txt
	awk '{print $1}' bowtiealign${i}.2u.txt > bowtie${i}.2u.1.txt
	awk '{print $3}' bowtiealign${i}.2u.txt > bowtie${i}.2u.3.txt
done


# now deal with spike ins
for i in {1..26}; do 
	echo $i
	in=${SAMPLE}.${i}.quickout.txt 
	echo "handling spikeins ${in} ...}"	
	
	grep 'CGTCAGTC$' $in > ${SAMPLE}.spikes.${i}.quickprocessed.txt 
	awk '{print $1}' ${SAMPLE}.spikes.${i}.quickprocessed.txt > ${SAMPLE}.spikes.${i}.counts.txt 
	awk '{print $2}' ${SAMPLE}.spikes.${i}.quickprocessed.txt > ${SAMPLE}.spikes.${i}.seq.txt  
	nl ${SAMPLE}.spikes.${i}.seq.txt | awk '{print ">" $1 "\n" $2}' > ${SAMPLE}.spikes.${i}.fasta2u.txt 
	bowtie-build -q ${SAMPLE}.spikes.${i}.fasta2u.txt indexes/spikes${i}fasta2u 
	bowtie -v 3 -p 10 -f --best -a indexes/spikes${i}fasta2u ${SAMPLE}.spikes.${i}.fasta2u.txt bowtiealignspikes.${i}.2u.txt 
	awk '{print $1}' bowtiealign.spikes.${i}.2u.txt > bowtiespikes.${i}.2u.1.txt
	awk '{print $3}' bowtiealign.spikes.${i}.2u.txt > bowtiespikes.${i}.2u.3.txt 
done

# L1 barcodes
for i in {1..26}; do
	echo $i
	in=${SAMPLE}.${i}.quickout.txt
	echo "handling L1 ${in} ...}"	
	
	#split off real barcodes from spike-ins
	grep -v 'CGTCAGTC$' $in | grep '[AG][AG]$' > ${SAMPLE}.BC${i}.quickprocessed.L1.txt
	awk '{print $1}' ${SAMPLE}.BC${i}.quickprocessed.L1.txt > ${SAMPLE}.${i}.counts.L1.txt
	awk '{print $2}' ${SAMPLE}.BC${i}.quickprocessed.L1.txt > ${SAMPLE}.${i}.seq.L1.txt
	
	nl ${SAMPLE}.${i}.seq.L1.txt | awk '{print ">" $1 "\n" $2}' > ${SAMPLE}.BC${i}.fasta2u.L1.txt
	bowtie-build -q ${SAMPLE}.BC${i}.fasta2u.L1.txt indexes/BC${i}fasta2uL1
	bowtie -v 3 -p 10 -f --best -a indexes/BC${i}fasta2uL1 ${SAMPLE}.BC${i}.fasta2u.L1.txt bowtiealign.${i}.2u.L1.txt
	awk '{print $1}' bowtiealign.${i}.2u.L1.txt > bowtie${i}_2u_1L1.txt;awk '{print $3}' bowtiealign.${i}.2u.L1.txt > bowtie.${i}.2u.3.L1.txt
done


