
# Data Processing Information for {{ conf['project']['project_id'] }}

## Output Data Description

### Matrix information

{{ conf['project']['project_id'] }}.vbctable.tsv:		All viral barcode (VBC) sequences, broken down by type (real, spike-in, L1) and target region.  
			Barcode sequences may appear more than once, if found in multiple target areas. \`\`  
umi\_count is the number of unique UMIs seen, with read\_count representing the total number of reads behind those UMIs (the allocation of reads to UMIs is not included). 

For each brain:  
\<brain\>.rbcm.tsv	raw (real) barcode matrix	  
Raw real barcode UMI counts from filtered VBC table.   

Only includes barcodes where:   
If the experiment has injection samples, then those areas must have a minimum of  \> {{conf['vbctable']['inj_min_umi']}} molecules (UMIs).   
At least one target area had \> {{conf['vbctable']['target_min_umi']}} molecules (UMIs)  

If the experiment had injection samples, then the matrix contains only barcodes present in both injection and target areas. If there 
Rows represent barcodes and columns represent areas. 

\<brain\>.sbcm.tsv	spike-in barcode matrix.   
			Spike in barcode UMI counts.   
			Rows represent barcodes and columns represent areas. 

\<brain\>.nbcm.tsv	normalized barcode matrix  
			Filtered barcode values normalized by spike-in counts per SSI area. 


### Pipeline Overview

1. Raw paired-end sequencing data is parsed and assembled into single blocks consisting of barcode, UMI, and SSI.   
   During this stage all resulting reads with 'N's are filtered, along with any reads containing sections longer than {{stats['fastq_filter']['max_repeats']}} bp containing the same base (since homopolymers introduce errors).    
2. We take the barcodes of all remaining reads, and perform an all-against-all alignment and find all sets of barcodes with Hamming distance of 3 or less from each other (based on 30 bp of barcode length). These sets are then used to 'collapse' the barcode in all the full sequences (VBC+UMI+SSI) to a single unique sequence, essentially 'fixing' replication and sequencing errors in the barcodes.   
3. We then take the full sequences and sort them all according to SSIs (which correspond to unique brain areas), and then determine the molecule number of a barcode in a certain brain area by counting UMIs. As this is done we aggregate the original read counts for each UMI.   
4. Within each SSI we separate out spike-in and real sequences based on their properties.   
5. We create the aggregate table \<experiment\>.vbctable.tsv above.   
6. The aggregate table can now be converted to connection matrices, with the number of each real and spike-in barcode in each brain area. In this matrix, each row is one barcode, each column is a brain area, and the element of the matrix corresponds to the molecule number of the barcode in the brain area.  
7. Lastly, we normalized the number of molecules of each barcode in each brain area to the total number of spike-in molecules in the corresponding brain area.  This is to compensate RT/PCR variations during sequencing library preparation.  This produces the normalized barcode matrix above. 

### Quality Control

 For this MAPseq data set, we have checked:

1. Sequencing depth for projection sites. 
   There are a minimum of {{ conf['vbctable']['target_min_umi'] }} reads for each barcode in target site samples and 2 reads for injection site samples. This sequencing depth is considered enough for target sites.

Spike-in counts are mostly uniform across all projections sites and all injection sites.  
We have seen 20 folds difference between injection sites and projection sites, and 2-3 folds difference within most projection sites, which is considered normal. 

1. We have estimated the fraction of neurons labeled with the same barcodes using barcode diversity (\~20M for the current library) and unique barcodes recovered. 

The estimated amount of reused barcodes within each brain and across different brains are listed in the form below:

|  | B 1 | B 2 |
| :---- | :---- | :---- |
| B 1 | 1 | 3 |
| B 2 | 3 | 4 |

1. We have checked barcode amount in H2O controls. 

The UMI counts from negative controls are listed in the matrix named “ctrl”. This background looks normal. 

1. We have calculated false positive rates. For target sites, the false positive rate is 0 when the UMI threshold is set above 3\.  
1. We have checked the template switching rate induced during sample processing (PCR). The rate is 0.08%, which means 0.08% of the total molecules in the target sites coming from the template switching. This rate is considered low.

### Further Analysis

What we haven’t done but MAPseq users should do is to:

1. Check distribution/mean barcode counts in projection areas: does it fit what's known?  
1. Population-level correlation with ABA or other bulk tracing data  
1. Check area-to-area correlations. They should fit what's known (e.g. cortical neurons mostly don't project to contralateral cortex and brain stem at the same time). Check for potential contaminations/fibers of passage/bad dissections.   
1. Check consistency across multiple brains: Are neurons mixed in t-sne? Are the min distance within a brain similar to min distance across brains?  
1. Clustering: do you get known cell types?

### Acknowledgement Policy

If the services from the MAPseq Core Facility are used to generate data in any publications, please acknowledge the MAPseq Core Facility in CSHL for processing samples for MAPseq sequencing.


