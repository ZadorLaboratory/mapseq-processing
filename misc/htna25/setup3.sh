#!/bin/bash
conda install -y pandas pyarrow ipython jupyter numpy scikit-learn matplotlib seaborn biopython scipy fastcluster openpyxl bowtie2 natsort git paramiko h5py dask networkx mergedeep
mkdir -p ~/git ; cd ~/git
git clone https://github.com/ZadorLaboratory/mapseq-processing.git
git clone https://github.com/ZadorLaboratory/mapseq-analysis.git
cd 
echo "run jupyter notebook, copy URL into browser, navigate to git/mapseq-processing/notebooks"
