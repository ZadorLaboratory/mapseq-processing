#	MAPseq Processing
This protocol will guide you through the steps for turning raw sequencing reads into a viral barcode matrix, which you can then analyze yourself or in conjunction with BARseq data.

The pipeline consists of four command-line utilities. Below is a schematic with a summary of what each tool does. Initial input is the raw paired-end sequencing files in FASTQ format. Final output is a normalized viral barcode matrix suitable for further analysis.

## Install software and dependencies

These instructions assume familiarity with running bioinformatics pipelines...

* Install Conda. 
[https://docs.conda.io/projects/miniconda/en/latest/index.html](https://docs.conda.io/projects/miniconda/en/latest/index.html)

* Create an environment for the MAPseq pipeline.

```
conda create -n mapseq python==3.9 
```

* Activate the environment

```
conda activate mapseq
```

* Install additional Conda repositories

```
conda config --add channels conda-forge
conda config --add channels bioconda
```
* Install dependencies, useful tools, confirm update

```
conda install -y pandas ipython jupyter numpy scikit-learn matplotlib seaborn biopython scipy matplotlib openpyxl bowtie bowtie2 natsort kneed git
conda update scipy 
```
* Clone the mapseq-processing software from the repository to the standard location. (This assumes you already have git installed. If not, install it first). 

```
mkdir ~/git
git clone https://github.com/ZadorLaboratory/mapseq-processing.git 
```
All the code is currently organized so it is run directly from the git directory (rather than installed). 
* Create a working directory for your experiment, and copy in a metadata file and the fastq sequencing data files, and the default configuration file. 

```
mkdir \~/mapseq/M205 ; cd ~/mapseq/M205
mkdir fastq ; cp /some/path/<fastqfiles>  fastq/ 
cp /some/path/EXP_sampleinfo.xlsx ./ 
cp \~/git/mapseq-processing/etc/mapseq.conf ./
```

## Experiment Metadata, Initial Configuration
All the commands in the pipeline look up metadata about the experiment from an Excel spreadsheet. We have provided the relevant spreadsheet for our example data. All that is absolutely required for the pipeline to work is that there be a worksheet entitled "Sample Information", and that there be columns titled:
"RT primers for MAPseq", 
"Brain", 
"Region",
"Our Tube #", 
"Site information"

The Site information column should indicate whether the dissected area is a target or injection, and additionally whether it is a negative (biological control) or control (water). 
Region and Brain should be obvious. Note that the region labels will be used (if available) for the default plots.  

An example spreadsheet is included in 
~/git/mapseq-processing/etc/M001_sampleinfo.xlsx

The RT primer barcode table is included in the software, as ~/git/mapseq-processing/etc/barcodes_v2.txt

By default, commands in the pipeline will take their defaults from a single configuration file, included in the distribution ~/git/mapseq-processing/etc/mapseq.conf.
Often, these can be overridden from the command line. You can also copy the default configuration file, edit it, and use it from each command with the -c <configfile> switch. 

For each command, there are typical arguments. Additional arguments may be supported. Run the -h switch to see the full usage help. 

### process_fastq_pairs.py

```
~/git/mapseq-processing/scripts/process_fastq_pairs.py 
	-v  				# give verbose output.  
	-o fastq.out/M211.all.fasta 	# write to FASTA file
	fastq/M211_HZ_S1_R*.fastq	# paired-end input FASTQs. 
```
process_fastq_pairs.py pulls out sequences lines from both read files, trims, and combines them. Some QC filtering is done at this step, with removal of long stretches of homo-polymers and removal of any sequences with low-confidence base calls ('N' bases).

For a standard ~400M read experiment, this takes about 45 minutes.  

### align_collapse.py

```
~/git/mapseq-processing/scripts/align_collapse.py
	-v 				# give verbose output
	-b 30 				# viral barcode length=30
	-m 3 				# maximum Hamming distance of 3 
	-O collapse.out		# output directory
	fastq.out/M211.all.fasta	# input FASTA
```

align_collapse.py separates out the viral barcode part of the data, and partitions them into groups where all the members are within a Hamming distance of 3 of each other. It  then sets all members of the partitions to have the same sequence. (This sequence is a randomly chosen one from the original partition, because it doesn't matter what the exact sequence is, as long as they all share it, and it is unique within the data). 

This is necessary because there is a fair amount of mutation that occurs during viral replication. We don't worry about mismatches in the other components of the sequence (SSI barcode, UMI sequence) because these are added during sample processing and thus have lower error rates.

For a standard ~400M read experiment, this takes about 90 minutes on a high-memory (192GB RAM) node.  

### process_fasta.py

```
~/git/mapseq-processing/scripts/process_fasta.py 
	-v						# give verbose output
	-s M211_sampleinfo.xlsx 			# sample metadata
	-O fasta.out 					# output directory 
	collapse.out/M205.all.collapsed.fasta	# adjusted FASTA
```

Now we break the data into region by separating based on the SSI barcode.

This takes about 50 minutes for 400M reads and 10 SSI barcodes.  

### process_ssifasta.py
This program actually calls a sub-program that handles each SSI  barcode in a separate process, allowing you to leverage all the CPUs on your system and finish faster. 

```
~/git/mapseq-processing/scripts/process_ssifasta.py 
	-v				# give verbose output
	-n 				# do not collapse (already done earlier)
	-t 24				# processes to use 
	-s M211_sampleinfo.xlsx 	# sample metadata 
	-o M211.merged.all.tsv		# output file for all VBCs 
	-O ./ssifasta.out 		# output directory
	fasta.out/BC*.fasta		# input FASTA, one for each SSI barcode 
```

This reads in all the separate SSI barcode FASTAs, and further separates out real VBC (viral barcodes), spike-ins (used for normalization), and L1s (used to detect template switching--a QC measure). These are all then merged and labelled into a single TSV with all the data. This includes molecule counts per VBC

This typically takes about 12 minutes. 

### process_merged.py
```
~/git/mapseq-processing/scripts/process_merged.py 
	-v 
	-s M211_sampleinfo.xlsx  
	-e M211.gac 
	-O merged.out 
	ssifasta.out/M211.merged.tsv
```

This program produces, for each brain, three matrices of viral barcodes by dissected region. This is the step at which (optionally) various threshold can be applied. 
-- Minimum UMI count
-- Filter out VBCs not present in the injection
-- Minimum original read count.   

The three matrices are a real barcode matrix (rbcm), a real barcode matrix normalized by spike-ins (nbcm), and a scaled+normalized barcode matrix (scbcm). The real barcode matrix contains raw numbers of molecules (UMIs) for each VBC. The normalized matrix adjusts the raw numbers by spike-in amounts. This normalization is done such that the real counts in each column are weighted by the total number of spike-in molecules in that column, in such a way that the weights are always 1.0 or greater. A third, normalized+scaled matrix is useful for drawing heatmaps, since all values are scaled back to values from 0.0 to 1.0 so they can be represented by color intensity.      

These matrices, the all.merged.tsv file, and various statistics gathered along the way are the full output of the processing phase, implemented in the mapseq-processing project. Code useful for further investigation and inference is contained in the mapseq-analysis project. 

