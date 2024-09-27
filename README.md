# NASTRA
This is the code of the article *Innovative Short Tandem Repeat Analysis through Cluster-Based Structure-Aware Algorithm in Nanopore Sequencing Data*. 

In this repository, we provide a STR genotyping method for forensic-level accuracy.

## Table of content
- [NASTRA](#nastra)
  - [Table of Contents](#table-of-contents)
  - [Installation and Test](#installation-and-test)
  - [Usage](#usage)
  - [Citation](#citation)

## Installation and Test

In this part, we introcduce the installation step and test the code for NASTRA. We use `conda` to manage the computing environment, and install samtools v1.9 firstly via the `bioconda` channel. 

Once finish the installtion step, we provide the alignment files of eight standard cell lines for code test (barcode01-08). Additionally, the files are downsampled using [NanoTime](https://github.com/renzilin/NanoTime), which is described in the article.


```bash
# download the source code, the size of test data is ~67Mb, this may be time-consuming
git clone https://github.com/renzilin/NASTRA.git

# enter the directory
cd NASTRA

# install dependencies. Only test on Ubuntu System (samtools can't be install via conda for ARM)
conda env create -n nastra_env -f env.yaml

# code test
conda activate nastra_env

# sample-level parallelization using xargs
seq -w 01 08 | xargs -P 8 -L 1 -I {} python NASTRA/nastra.py call -b test_data/alignment_8standards_1h/barcode{}.bam -o test_data/nastra_out/barcode{}.txt;
```
The sample information of `BAM` files is listed below:
| File name      | Sample   |
|----------------|----------|
| barcode01.bam  | NA12878  |
| barcode02.bam  | NA24143  |
| barcode03.bam  | NA24149  |
| barcode04.bam  | NA24694  |
| barcode05.bam  | NA24695  |
| barcode06.bam  | 9947A    |
| barcode07.bam  | 9948     |
| barcode08.bam  | 2800M    |


## Usage
The `nastra.py` requires several options:
1. **required**: the path of BAM file, usually generated by `minimap2`.
2. **required**: the path of output file.
3. **optional**: the path of factsheet, storing the information of repeat motifs of STRs. We provide a revised file `NASTRA/cfgs/repeat_structure.pat`. Assessible in [NIST STRBase FactSheet](https://strbase.nist.gov/Loci/FactSheet).
4. **optional**: the path of panel file, storing the positional information of STRs. We provide a revised file `NASTRA/cfgs/panel_forenseq.csv`.
5. **optional**: the path of config file of SN, storing the optimal thresholds for NASTRA. We provide one in `NASTRA/cfgs/threshold.cfg`, which was utilized in our manusctipts. 
6. **optional**: the predefined SN threshold used for quality control, default is 25.
7. **optional**: the path of `samtools`, which is required to extract aligned reads.

```bash
usage: nastra.py call [-h] -b BAM [-f FACTSHEET] [-p PANEL] [-c CONFIG] -o OUTPUT [--samtools SAMTOOLS] [--sncutoff SNCUTOFF]

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     the path of bam file
  -f FACTSHEET, --factsheet FACTSHEET
                        the path of STR factsheet
  -p PANEL, --panel PANEL
                        the path of the loci panel file
  -c CONFIG, --config CONFIG
                        the config file of SN cutoff
  -o OUTPUT, --output OUTPUT
                        the path of output
  --samtools SAMTOOLS   the PATH to samtools

```


## Citiing NASTRA
Zilin Ren, Jiarong Zhang, Yixiang Zhang, Tingting Yang, Pingping Sun, Jiguo Xue, Xiaochen Bo, Bo Zhou, Jiangwei Yan, Ming Ni, NASTRA: accurate analysis of short tandem repeat markers by nanopore sequencing with repeat-structure-aware algorithm, Briefings in Bioinformatics, Volume 25, Issue 6, November 2024, bbae472, https://doi.org/10.1093/bib/bbae472
