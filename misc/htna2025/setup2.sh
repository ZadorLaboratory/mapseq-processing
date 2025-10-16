#!/bin/bash
cd ~/mapseq
mkdir M205 ; cd M205
wget https://raw.githubusercontent.com/ZadorLaboratory/mapseq-processing/refs/heads/main/misc/htna2025/M205.mapseq.conf -O M205.mapseq.conf
wget https://raw.githubusercontent.com/ZadorLaboratory/mapseq-processing/refs/heads/main/misc/htna2025/M205.sampleinfo.xlsx -O M205.sampleinfo.xlsx
mkdir fastq ; cd fastq
wget https://labshare.cshl.edu/shares/mbseq/mapseq/workshops/htna2025/fastq/M205_HZ_S1_R1_001.fastq.gz -O M205_HZ_S1_R1_001.fastq.gz
wget https://labshare.cshl.edu/shares/mbseq/mapseq/workshops/htna2025/fastq/M205_HZ_S1_R2_001.fastq.gz -O M205_HZ_S1_R2_001.fastq.gz
cd ../..
mkdir output ; cd output
wget https://labshare.cshl.edu/shares/mbseq/mapseq/workshops/htna2025/M205.htna25.tgz
tar -xvzf M205.htna25.tgz

cd ~/mapseq
conda update -y -n base -c defaults conda
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -y -n mapseq python==3.12
echo "activate mapseq env and run setup3.sh"