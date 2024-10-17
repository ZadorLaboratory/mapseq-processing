#!/bin/bash
cd ~/mapseq
wget https://labshare.cshl.edu/shares/mbseq/mapseq/htna2024/M205_sampleinfo.xlsx
mkdir fastq ; cd fastq
wget https://labshare.cshl.edu/shares/mbseq/mapseq/htna2024/fastq/M205_HZ_S1_R1_001.fastq.gz
wget https://labshare.cshl.edu/shares/mbseq/mapseq/htna2024/fastq/M205_HZ_S1_R2_001.fastq.gz
