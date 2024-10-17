#!/bin/bash
cd ~/mapseq
wget https://labshare.cshl.edu/shares/mbseq/mapseq/htna2024/M205_sampleinfo.xlsx
mkdir fastq ; cd fastq
wget https://labshare.cshl.edu/shares/mbseq/mapseq/htna2024/fastq/M205_HZ_S1_R1_001.fastq.gz
wget https://labshare.cshl.edu/shares/mbseq/mapseq/htna2024/fastq/M205_HZ_S1_R2_001.fastq.gz
conda update -y -n base -c defaults conda
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -y -n mapseq python==3.12
conda activate mapseq
conda install -y pandas pyarrow ipython jupyter numpy scikit-learn matplotlib seaborn biopython scipy fastcluster openpyxl bowtie2 natsort git paramiko
mkdir -p ~/git ; cd ~/git
git clone https://github.com/ZadorLaboratory/mapseq-processing.git
git clone https://github.com/ZadorLaboratory/mapseq-analysis.git
cd ~/mapseq
