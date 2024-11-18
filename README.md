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
# For Linux and Intel MacOS
	conda install -y pandas pyarrow ipython jupyter numpy scikit-learn matplotlib seaborn biopython scipy fastcluster matplotlib openpyxl bowtie bowtie2 natsort git dask  
 
# For MacOS with MX chips, bowtie not made, git part of x-code, and bowtie2 only available via brew ( https://brew.sh/ )
    conda install -y pandas pyarrow ipython jupyter numpy scikit-learn matplotlib seaborn biopython scipy fastcluster matplotlib openpyxl natsort dask
 
```

* Clone the mapseq-processing software from the repository to the standard location. (This assumes you already have git installed. If not, install it first). 

```
mkdir ~/git
git clone https://github.com/ZadorLaboratory/mapseq-processing.git 
```
All the code is currently organized so it is run directly from the git directory (rather than installed). 
* Create a working directory for your experiment, and copy in a metadata file and the fastq sequencing data files, and optionally the default configuration file. 

```
mkdir \~/mapseq/M205 ; cd ~/mapseq/M205
mkdir fastq ; cp /some/path/<fastqfiles>  fastq/ 
cp /some/path/EXP_sampleinfo.xlsx ./ 
cp \~/git/mapseq-processing/etc/mapseq.conf ./
```

## Experiment Metadata, Initial Configuration, Default Behavior
All the commands in the pipeline look up metadata about the experiment from an Excel spreadsheet. We have provided the relevant spreadsheet for our example data. All that is absolutely required for the pipeline to work is that there be a worksheet entitled "Sample Information", and that there be columns titled:
"RT primers for MAPseq", 
"Brain", 
"Region",
"Our Tube #", 
"Site information"

The Site information column should indicate whether the dissected area is a target or injection, and additionally whether it is a negative from a test brain (biological control), negative from a non-test brain, or control (water). 
Region and Brain should be obvious. Note that the region labels will be used (if available) for the default plots.  

An example spreadsheet is included in 
~/git/mapseq-processing/etc/M001_sampleinfo.xlsx

The RT primer barcode table is included in the software, as ~/git/mapseq-processing/etc/barcodes_v2.txt

By default, commands in the pipeline will take their defaults from a single configuration file, included in the distribution ~/git/mapseq-processing/etc/mapseq.conf.
Often, these can be overridden from the command line. You can also copy the default configuration file, edit it, and use it from each command with the -c <configfile> switch. 

For each command, there are typical arguments. Additional arguments may be supported. Run the -h switch to see the full usage help.

Note that for version 2, all the commands automatically generate Parquet file output (in addition to the standard TSV format). For full-sized datasets on moderate-memory systems (less than 128G) TSV files become problematic, causing out of memory errors. In those cases, the parquet file should be used as the input for the next step, as it will load faster and use less memory.  

For most of the commands, specifying the output file (with subdirectory) should be sensible. The programs will make needed sub-directories.  
 

### process_fastq_pairs.py

```
~/git/mapseq-processing/scripts/process_fastq_pairs.py 
	-v  				# give verbose output.  
	-o fastq.out/M253.reads.tsv 	# write to TSV file
	fastq/M253_CR_S1_R*.fastq	# paired-end input FASTQs. 
```
process_fastq_pairs.py pulls out sequences lines from both read files, trims, and combines them. Some QC filtering is done at this step, with removal of long stretches of homo-polymers and removal of any sequences with low-confidence base calls ('N' bases). Note that input files can be fastq or fastq.gz--the program will decompress on the fly if needed.

Program will accept an even-number set of inputs. Just be sure that they are in the correct read1-read2-read1'-reqd2' order.  

For a standard ~400M read experiment, this takes about 45 minutes. (Or 15 minutes on a new Mac M3). 

### aggregate_reads.py

```
~/git/mapseq-processing/scripts/aggregate_reads.py 
	-v  						# give verbose output.
	-t ~/scratch					# Local, fast, roomy place for DASK to put temp files.
	-m 2						# Minimum reads to retain data.    
	-o aggregated.out/M253.aggregated.tsv 	# write to TSV file
	fastq.out/M253.reads.tsv			# all reads merged to single sequence 
```
aggregate_reads.py takes all the assembled paired end reads (typically 52 nucleotides), and aggregates them by exact sequence, and sets a read_count column to the count of that unique sequence. Optionally can threshold by read count. 

This step is isolated from others because it uses Dask, a Python framework to allow working on very large datasets with limited amounts of memory. It leverages disk space for temporary storage instead of loading all data into memory. 

### filter_split.py

```
~/git/mapseq-processing/scripts/filter_split.py 
	-v  					# give verbose output.  
	-o filtered.out/M253.filtered.tsv 	# write to TSV file
	aggregated.out/M253.aggregated.tsv	# aggregated full reads (with read_count). TSV or Parquet 
```
filter_split.py performs sequence-level QC, removing reads that are likely to be erroneous. And it splits the reads into columns that are relevant for the MAPseq protocol.  
  

### align_collapse.py

```
~/git/mapseq-processing/scripts/align_collapse.py
	-v 					# give verbose output
	-m 3 					# maximum Hamming distance of 3 
	-o collapse.out/M253.collapsed.tsv 	# output TSV 
	filtered.out/M253.filtered.tsv	# Filtered and split data. TSV or Parquet 
```

This program takes the viral barcode part of the reads (vbc_read column), and partitions them into groups where all the members are within a Hamming distance of 3 of each other. It does this by running an all-by-all alignment of all unique viral barcode sequences. Self-matches are discarded. The remainder are used to create an edge graph (nodes are sequences, edges are between sequences less than 3 edits apart). This is then given to Tarjan's algorithm to determine "components", sets of sequences connected by edges. It then sets all members of the components to have the same sequence. (This sequence is the one in the original set with the most reads, but it shouldn't matter what the exact sequence is, as long as they all share it, and it is unique within the data). 

This processing step is necessary because there is a fair amount of mutation that occurs during viral replication. We don't worry about mismatches in the other components of the sequence (SSI barcode, UMI sequence) because these are added during sample processing and thus have lower error rates.

For a standard ~400M read experiment, this takes about 70 minutes on a high-memory (192GB RAM) node. (Or 50 minutes on a new Mac M3, 96GB RAM).  


### make_readtable.py

```
~/git/mapseq-processing/scripts/make_readtable.py 
	-v						# give verbose output
	-s M253_sampleinfo.xlsx 			# sample metadata
	-o readtable.out/M253.readtable.tsv 		# output TSV 
	collapse.out/M253.collapsed.tsv		# collapsed reads/counts (or parquet)
```

Now we break the data into useful sequence regions, and fill in addtional relevant information. Note that this is the only command that requires the sample information from the Excel spreadsheet (or derived TSV). 
 
We also create site-specific read count shoulder plots, in order to confirm that our read count thresholds are reasonable. 

Takes about 45 minutes for 400M reads on MacBook Pro M3 
  

### make_vbctable.py

```
~/git/mapseq-processing/scripts/process_ssifasta.py 
	-v					# give verbose output
	-o vbctable.out/M253.vbctable.tsv 	# output file for all VBCs 
	readtable.out/M253.readtable.tsv 	# fully populated read table. 
```
Consumes the readtable, drops all reads with missing/unmatched tags/sequences, drops reads that do not meet the read count threshold for the site type. It then collapses by viral barcode, calculating UMI count for each VBC.
Outputs a VBC-oriented table. 

This typically takes about 12 minutes. (3 minutes on Mac M3)

### make_matrices.py
```
~/git/mapseq-processing/scripts/make_matrices.py 
	-v 
	-e M253 
	-O matrices.out 
	vbctable.out/M253.vbctable.tsv
```

This program produces, for each brain, several matrices of viral barcodes by dissected region. This is the step at which (optionally) various threshold can be applied. 
-- Minimum UMI count in target and/or injection areas (min_target_umi, min_inj_umi)
-- Filter out VBCs not present in the injection site (require_injection)

The matrices produced are: 

A real barcode matrix (rbcm). The real barcode matrix contains raw numbers of molecules (UMIs) for each VBC.

A spike-in barcode matrix (sbcm). Contains raw numbers of spike-in molecules for each site. 

A real barcode matrix normalized by spike-ins (nbcm). The normalized matrix adjusts the raw numbers by spike-in amounts. This normalization is done such that the real counts in each column are weighted by the total number of spike-in molecules in that column, in such a way that the weights are always 1.0 or greater. For most further analysis, this can be considered the canonical processing output file. 

A scaled+normalized barcode matrix (scbcm). This normalized+scaled matrix is useful for drawing heatmaps, since all values are scaled back to values from 0.0 to 1.0 so they can be represented by color intensity.      


Code useful for further investigation,  inference, and visualization is contained in the mapseq-analysis project:
	[https://github.com/ZadorLaboratory/mapseq-analysis](https://github.com/ZadorLaboratory/mapseq-analysis) 

