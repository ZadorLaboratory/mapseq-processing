#!/bin/bash -x
#example preprocessing of sequencing libary M205

FQ1="M205_HZ_S1_R1_001"
FQ2="M205_HZ_S1_R2_001"
SAMPLE="M205_HZ_S1"

cd barcodesplitter

BCidx=($(seq 0 1 81; seq 97 1 177; seq 193 1 233; seq 249 1 279)) #the first number should be n-1
BCidx=($(seq 0 1 26)) 
for i in {1..26}; do 
    #this number should be exactly the same as the total RT primer number used in BCidx
	#filter out reads with Ns, cut off indexes and unique datafiles
	awk "NR%2==0" BC${BCidx[$i]} | grep -v N | cut -b 1-44 | sort | uniq -c | sort -nr > BC${BCidx[$i]}.44.txt
	#split output files into two files per index, one that is containing the read counts of each unique sequence, the other the unique sequences themselves.
	awk '{print $1}' BC${BCidx[$i]}.44.txt > BC${BCidx[$i]}.44.counts.txt
	awk '{print $2}' BC${BCidx[$i]}.44.txt > BC${BCidx[$i]}.44.seq.txt
done



#cd barcodesplitter/thresholds

cd thresholds
 
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
    	echo "awk -v j=$j '$1<j {print NR}' ../BC${BCidx[$i]}.44.counts.txt | head -n 1 > t$i"
		awk -v j=$j '$1<j {print NR}' ../BC${BCidx[$i]}.44.counts.txt | head -n 1 > t$i
		thresh=$(cat t$i)
		echo "thresh = ${thresh}"
		echo "head  -n $thresh ../BC${BCidx[$i]}.seq.txt | cut -b 1-32 | sort | uniq -c | sort -nr > BC${i}.32.txt"
		head  -n $thresh ../BC${BCidx[$i]}.44.seq.txt | cut -b 1-32 | sort | uniq -c | sort -nr > BC${i}.32.txt
    else
		grep -nr ^${threshold[$i]}$  ../BC${BCidx[$i]}.44.counts.txt -m 1 | cut -f1 -d":" > t$i
		thresh=$(cat t$i)
		echo "thresh = ${thresh}"
		head ../BC${BCidx[$i]}.44.seq.txt -n $thresh | cut -b 1-32 | sort | uniq -c | sort -nr > BC${i}.32.txt
    fi
done

#ok, thats a lot of preprocessing done. Now we have to error correct these sequences. 
# The all against all mapping of barcodes using bowtie is done on the command line. 
# after this we move into MATLAB.

mkdir indexes

for i in {1..26}
do
	echo $i
	in=BC${i}.32.txt
	echo "handling real ${in} ...}"
	#split off real barcodes from spike-ins

	grep -v 'CGTCAGTC$' $in | grep '[TC][TC]$' > BC${i}.real.txt
	awk '{print $1}' BC${i}.real.txt > BC${i}.counts.txt
	awk '{print $2}' BC${i}.real.txt > BC${i}.seq.txt

	nl BC${i}.seq.txt | awk '{print ">" $1 "\n" $2}' > BC${i}.real.fasta 
	bowtie-build -q BC${i}.real.fasta indexes/BC${i}real 
	bowtie -v 3 -p 10 -f --best -a indexes/BC${i}real BC${i}.real.fasta BC${i}.real.bowtie
	awk '{print $1}' BC${i}.real.bowtie > BC${i}.real.read.txt
	awk '{print $3}' BC${i}.real.bowtie > BC${i}.real.align.txt
done


# now deal with spike ins
for i in {1..26}; do 
	echo $i
	in=BC${i}.32.txt 
	echo "handling spikeins ${in} ...}"	
	
	# CGTCAGTC
	grep 'CGTCAGTC$' $in > BC${i}.spike.txt 
	awk '{print $1}' BC${i}.spike.txt > BC${i}.spike.counts.txt 
	awk '{print $2}' BC${i}.spike.txt > BC${i}.spike.seq.txt  

	nl BC${i}.spike.seq.txt | awk '{print ">" $1 "\n" $2}' > BC${i}.spike.fasta 
	bowtie-build -q BC${i}.spike.fasta indexes/BC${i}spike 
	bowtie -v 3 -p 10 -f --best -a indexes/BC${i}spike BC${i}.spike.fasta BC${i}.spike.bowtie 
	awk '{print $1}' BC${i}.spike.bowtie > BC${i}.spike.read.txt
	awk '{print $3}' BC${i}.spike.bowtie > BC${i}.spike.align.txt 
done

# L1 barcodes
for i in {1..26}; do
	echo $i
	in=BC${i}.32.txt 
	echo "handling L1 ${in} ...}"	
	
	grep -v 'CGTCAGTC$' $in | grep '[AG][AG]$' > BC${i}.real.L1.txt 
	awk '{print $1}'  BC${i}.real.L1.txt > BC${i}.real.L1.counts.txt
	awk '{print $2}'  BC${i}.real.L1.txt > BC${i}.real.L1.seq.txt
	
	nl BC${i}.real.L1.seq.txt | awk '{print ">" $1 "\n" $2}' > BC${i}.real.L1.fasta
	bowtie-build -q BC${i}.real.L1.fasta indexes/BC${i}L1
	bowtie -v 3 -p 10 -f --best -a indexes/BC${i}L1 BC${i}.real.L1.fasta BC${i}.real.L1.bowtie
	awk '{print $1}' BC${i}.real.L1.bowtie > BC${i}.real.L1.read.txt
	awk '{print $3}' BC${i}.real.L1.bowtie > BC${i}.real.L1.align.txt
done

