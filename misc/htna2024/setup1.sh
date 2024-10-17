#!/bin/bash
#
# Set up conda
cd ; mkdir mapseq ; cd mapseq
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh -b -u 
~/miniconda/bin/conda init bash
wget https://raw.githubusercontent.com/ZadorLaboratory/mapseq-processing/refs/heads/main/misc/htna2024/setup2.sh
chmod +x setup2.sh
exit


