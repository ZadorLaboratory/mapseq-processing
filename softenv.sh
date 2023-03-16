#!/usr/bin/env -S bash -l
set -e
# install conda
if ! command -v conda&> /dev/null; then
	echo "installing miniconda..."
	wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh
	bash Miniconda3-py39_4.9.2-Linux-x86_64.sh -b
	rm -f Miniconda3-py39_4.9.2-Linux-x86_64.sh
	~/miniconda3/bin/conda init 
	echo "miniconda installed. restart terminal."
	exit 0
else
	echo "miniconda installed already."
fi
# be sure to restart terminal to allow conda to start

conda create -y -n mapseq python=3.9
conda update -n base conda
conda activate mapseq

conda config --add channels conda-forge
conda config --add channels bioconda
conda install bowtie2  
conda install biopython=1.79 ipython=8.2.0 numpy=1.22.3 pandas=1.4.2 scipy=1.9.3 matplotlib jupyter seaborn openpyxl natt
conda install natsort
# conda install -c bioconda fastx_toolkit
