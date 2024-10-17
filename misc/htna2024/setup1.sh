#!/bin/bash
#
# Set up conda
cd ; mkdir mapseq ; cd mapseq
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh -b -u 
~/miniconda/bin/conda init bash
wget https://raw.githubusercontent.com/ZadorLaboratory/mapseq-processing/refs/heads/main/misc/htna2024/setup2.sh -O setup2.sh
chmod +x setup2.sh
wget https://raw.githubusercontent.com/ZadorLaboratory/mapseq-processing/refs/heads/main/misc/htna2024/setup3.sh -O setup3.sh
chmod +x setup3.sh
echo "restart shell and run setup2.sh"


