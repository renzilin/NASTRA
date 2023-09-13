# NASTRA
Innovative Short Tandem Repeat Analysis through Cluster-Based Structure-Aware Algorithm in Nanopore Sequencing Data


## Install

```bash
git clone https://github.com/renzilin/NASTRA.git # the size of test data is ~67Mb

cd NASTRA

# only test on Ubuntu

conda create -n nastra_env python=3.8 -y  

conda install -c bioconda samtools

pip install tqdm numpy pandas parasail

```


## Uasge

```bash
seq -w 01 08 | xargs -P 8 -L 1 -I {} python NASTRA/nastra.py call -b test_data/alignment_8standards_1h/barcode{}.bam -o test_data/nastra_out/barcode{}.txt;
```