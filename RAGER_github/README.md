# **RAGER: A User-Friendly Computational Platform for Integrated Analysis of RNA-Seq and ATAC-Seq Data**

## Introduction to RAGER:
**RAGER** is a computational platform, that integrates the popular bioinformatics tools in an automated thread for joint mining of RNA-seq and ATAC-seq data. RAGER facilitates integrative analysis of transcriptome and chromatin accessibility by providing an automated workflow that minimizes the need for bioinformatics expertise and significantly reduces processing time. We demonstrate RAGER's utility for novel biological discovery by characterizing the transcriptome and chromatin accessibility of two recently published datasets.

## Table of Contents
1. [Quick start]()
2. [Preprocess RNAseq data]("https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Preprocess_RNAseq/RNAseq_analysis.md")
3. [Preprocess ATACseq data]() 
4. [Joint analysis]()
5. [Custom analysis]()

# **Quick start**
## System requirements:
Some of the tools that RAGER uses, e.g. Hisat2 and Bowtie2 are very memory intensive programs. Therefore we recommend the following system requirements for VIPER:

### Minimal system requirements:
- **Operating System:** Linux (tested on Ubuntu 20.04)  
- **Memory:** ≥ 16 GB RAM  
- **CPU:** ≥ 8 cores 
- **Storage:** ≥ 500 GB available space  

### Recommended system requirements:
We recommend that you have at least 1T of ram and at least a 20-core CPU if you want to run RAGER in multi-threaded mode (which will speedup the workflow significantly). Our own servers equipped with 192 CPU cores, 1.5 TB of RAM, 80 TB of storage, and running CentOS Stream.


## Getting Started - Install the Conda virtual environment and the necessary data: 

### Install virtual environment.
1. Install conda. Download the conda installer from https://www.anaconda.com/. 
```
bash <conda-installer-name>-latest-Linux-x86_64.sh
conda --version
```
(onda-installer-name)This is a placeholder. During actual use, it needs to be replaced with the specific installation package name.The conda version should appear if it has been installed correctly.

2. install git.
```
git --version
```
If it does not return the installed version of git, follow https://github.com/git-guides/install-git to install git.

3. Clone the branch from the GitHub repository which contains the RAGER code and the necessary files for this pipeline.
```
mkdir PROJECT
cd PROJECT
git clone https://github.com/bioinfo202408/RAGER
```
you can name your PROJECT folder anything

4. Create RAGER virtual environment.
```
conda env create -f RAGER.yml
```
We have provided a complete Conda environment file (E3RC.yml) in the GitHub repository, which includes all the dependencies required for RAGER.Users can download this file, optionally edit the environment name parameter (default "E3RC"), and create the environment using conda. Besides, **MEME Suite AME (v5.5.7)**, which required root permission to install, is also used for analysis. Please follow https://meme-suite.org/meme/meme_5.5.7/ to install MEME Suite AME (v5.5.7).

## Download static reference files
RAGER is dependent on reference files which can be found for the supported species listed below: download link [hg38](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/) [mm10](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/)
1. If you are analyzing human data
```
#Download and unzip reference geneme
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.fa.gz

gunzip GRCh38.p14.genome.fa.gz
#Index reference geneme
bowtie2-build --threads 5 GRCh38.p14.genome.fa /home/yjliu/rager_case2/genomeanno/bowtie2Index/GRCh38

hisat2-build -p 5 GRCh38.p14.genome.fa /home/yjliu/rager_case2/genomeanno/hisat2Index/GRCh38
#Download and unzip annotation file
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

gunzip gencode.v44.annotation.gtf.gz
```
2. If you are analyzing mouse data
```
#Download and unzip reference geneme
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz

gunzip GRCm38.p6.genome.fa.gz
#Index reference geneme
bowtie2-build --threads 5 GRCm38.p6.genome.fa /home/yjliu/rager_case2/genomeanno/bowtie2Index/GRCm38

hisat2-build -p 5 GRCm38.p6.genome.fa /home/yjliu/rager_case2/genomeanno/hisat2Index/GRCm38
#Download and unzip annotation file
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz

gunzip gencode.vM25.annotation.gtf.gz
```

## Download or generate RNA-seq and ATACseq datasets
RAGER has been verified to have excellent functionality in datasets GSE85632 and GSE261119.
Raw sequencing data were downloaded from the GEO database (GSE85632) From this dataset, we selected specific SRR files (SRR4032350, SRR4032351, SRR4032352, SRR4032353,SRR4032269, SRR4032270, SRR4032271, SRR4032272) relevant to our analysis, resulting in a total of 16 fastq files.
Use the SRR4032350 file as an example.
```
prefetch SRR4032350 -O ../../datasets/mappingfile/srafile

fastq-dump --split-files -O ../../datasets/mappingfile/fastqfile ../../datasets/mappingfile/ srafile/SRR4032350.sra
```
**Users who have custom sequencing data only need to unzip the original data and rename it to a unified filename.**
```
gunzip -c test_R1.fq.gz > test_1.fastq

gunzip -c test_R2.fq.gz > test_2.fastq
```

## **The subsequent analysis steps should be referred to the readme files of each process.**
