#!/bin/bash
cd ~/mapseq
mkdir fastq ; cd fastq
wget https://labshare.cshl.edu/shares/mbseq/mapseq/htna2024/fastq/M205_HZ_S1_R1_001.fastq.gz -O M205_HZ_S1_R1_001.fastq.gz
wget https://labshare.cshl.edu/shares/mbseq/mapseq/htna2024/fastq/M205_HZ_S1_R2_001.fastq.gz -O M205_HZ_S1_R2_001.fastq.gz
cd ~/mapseq
wget https://raw.githubusercontent.com/ZadorLaboratory/mapseq-processing/refs/heads/main/misc/htna2024/M205.mapseq.conf
wget https://raw.githubusercontent.com/ZadorLaboratory/mapseq-processing/refs/heads/main/misc/htna2024/M205_sampleinfo.xlsx
conda update -y -n base -c defaults conda
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -y -n mapseq python==3.12
echo "activate mapseq env and run setup3.sh"