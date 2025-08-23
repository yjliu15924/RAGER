# **RAGER table of contents**
1. [Quick start]()
2. [Preprocess RNAseq data]("D:\desktop\RAGER\RAGER_github\Scripts\Preprocess_RNAseq\RNAseq_analysis.md")
3. [Preprocess ATACseq data]() 
4. [Joint analysis]()
5. [Custom analysis]()

## **List of processes**
- [RNAseq_quality_control]()
  - [Trim_galore]()
  - [Multiqc]()
- [RANseq_mapping_reads_to_genome]()
  - [Align_RNA-seq_reads_to_the_reference]()
  - [Filter_multi-mapped_reads]()
  - [Convert_sam_files_to_bam_format]()
  - [Sort_bam_files_by_genome_coordinates]()
  - [Mark_and_remove_duplicate_reads]()
  - [Index_bam_files]()
- [RSeQC_gene_body_coverage]()
- [RNAseq_alignment_rates]()
- [StringTie_expression_quantification]()
- [RANseq_extract_expression_data]()
- [RNAseq_cluster_analysis]()
- [RANseq_PCA_analysis]()
- [RNAseq_differential_expression_analysis]()
  - [Construct_samplelist]()
  - [Prepare_differential_expression_analysis]()
  - [Differential_expression_analysis]()

## **Preprocessing**

## RNAseq_quality_control

### **Trim_galore**

**Description**

Trim_galore is used to remove low quality bases and adapter sequences from the raw RNA-seq reads. It performs quality control by trimming low-quality ends from reads in addition to adapter removal.

**Inputs**
- Raw FASTQ files: `*_1_fq`, `*_2.fq` (for paired-end)
- Raw FASTQ files: `*_1_.fq` (for single-end)

**Parameters**
- `--quality`: Quality threshold for trimming low-quality ends (default: 25)
- `--paired33`:  Indicate that the input is paired-end data with Phred+33 quality score encoding.
- `--stringency`: Overlap with adapter sequence required to trim (default: 3)
- `--length`: Discard reads shorter than this length (default: 35)
- `--fastqc`: Quality Control Report

**Outputs**
- Trimmed FASTQ files: `*_1_val_1.fq`, `*_2_val_2.fq` (for paired-end)
- Trimmed FASTQ files: `*_1_trimmed.fq` (for single-end)
- FastqQC reports: `*fastqc.html`

**Output directory** 
- `./datasets/RNAseq/quality_control_file/`

### **Multiqc**
**Description**

MultiQC aggregates results from various bioinformatics analyses into a single report. It collects the output from Trim_galore and other QC tools to generate comprehensive quality control metrics.

**Input**
- FastqQC reports: `*fastqc.html`

**Outputs**  
- HTML report: `multiqc_report.html`
- Data directory: `multiqc_data/`

**Output directory**  
- `./datasets/RNAseq/quality_control_file/multiqc_report.html`

## RANseq mapping reads to genome

### **Align RNA-seq reads to the reference genome using HISAT2**

**Description**  
HISAT2 is used to align RNA-seq reads to a reference genome. It is designed specifically for RNA-seq data and can handle spliced alignments efficiently.

**Input**
- Trimmed FASTQ files: `*_1_val_1.fq`, `*_2_val_2.fq` (for paired-end)
- Trimmed FASTQ files: `*_1_trimmed.fq` (for single-end)

**Parameters**  
- `-x`: Path to HISAT2 index
- `-1`: Path to forward reads
- `-2`: Path to reverse reads (for paired-end)
- `-S`: Output SAM file
- `--dta`: Output alignments tailored for transcript assemblers
- `-p`: Number of threads
- `--summary-file`: Output alignments rate

**Outputs**  
- SAM files: `*accepted_hits.sam`
- Alignment summary: `mapping_summary.txt`

**Output directory**  
- `./datasets/RNAseq/hisat2file/`

### **Filter multi-mapped reads**

**Description**  
Use "grep" to process the alignment results and filter out reads that align to multiple genomic locations,retaining only uniquely mapped reads to improve downstream analysis accuracy.

**Parameters**  
- `-v`: nvert match — select lines that do not match the pattern.
- `-E`: Use extended regular expressions (ERE) for the pattern.
- `-w`: Match only whole words (pattern must match an entire word).

**Outputs**  
- Filtered BAM files: `*accepted_hits_NHi1.sam`

**Output directory**  
- `./datasets/RNAseq/hisat2file/`

### **Convert sam files to bam format using samtools**

**Description**

Samtools view converts SAM files to BAM format, which is a compressed binary version of the SAM format. This step reduces file size and prepares the data for downstream processing. It can also be used to filter reads based on various criteria.

**Input**
- SAM files: `*accepted_hits_NHi1.sam` 

**Parameters**
- `-S`: Input is SAM format
- `-o`: Output file name
- `-@`: Number of threads

**Outputs**
- BAM files: `*accepted_hits_NHi1.bam`

**Output directory**
- `./datasets/RNAseq/hisat2file/`

### **Sort bam files by genome coordinates**

**Description**  
This step sorts the BAM files by genomic coordinates, which is required for many downstream analyses and improves the efficiency of data access.

**Input**
- BAM files: `*accepted_hits_NHi1.bam`

**Parameters**  
- `-o`: Output file name
- `-@`: Number of threads

**Outputs**  
- Sorted BAM files: `*accepted_hits_NHi1.sorted.bam `

**Output directory**  
- `./datasets/RNAseq/hisat2file/`

### **Mark and remove duplicate reads using Picard**

**Description**  
Picard's MarkDuplicates identifies and flags or removes duplicate reads that may have resulted from PCR amplification during library preparation.

**INput**
- Sorted BAM files: `*accepted_hits_NHi1.sorted.bam `

**Parameters**  
- `-Xmx15g`: Set the maximum Java heap memory to 15 GB.
- `I`: Input BAM file
- `O`: Output BAM file
- `METRICS_FILE`: File to write duplicate metrics
- `REMOVE_DUPLICATES`: Whether to remove (true) or just mark (false) duplicates(default: true)
- `ASSUME_SORT_ORDER`: Specify the expected sort order of the input BAM file, skipping sort-order validation (default:coordinate)

**Outputs**  
- Deduplicated BAM files: `*[sampleID].bam`
- Duplicate metrics: `*.metricsFile`

**Output directory**  
- `./datasets/RNAseq/hisat2file/`

### **Index bam files using samtools**

**Description**  
This step creates an index for the BAM files, which allows for quick random access to the data, essential for many downstream analysis tools.
**Input**
- Deduplicated BAM files: `*[sampleID].bam`

**Outputs**  
- BAM index files: `*[sampleID].bam.bai`

**Output directory**  
- `./datasets/RNAseq/hisat2file/`

## RSeQC_gene_body_coverage

**Description**  
RSeQC's gene body coverage analysis evaluates the uniformity of coverage across gene bodies, which can help identify potential biases in the RNA-seq experiment.

**Input**
- All samples BAM files: `*[sampleID].bam`

**Parameters**  
- `-i`: Input BAM file
- `-r`: Reference gene model in BED format
- `-o`: Output prefix

**Outputs**  
- Gene body coverage plot: `*_geneBodyCoverage.curves.pdf`
- Gene body heatMap plot: `*_geneBodyCoverage.heatMap.pdf`

**Output directory**  
- `./datasets/RNAseq/`

## RNAseq_alignment_rates

**Description**  
This step will present the alignment rates of all samples in a stacked bar chart, providing a quick overview of the mapping quality.

**Input**
- Alignment summary: `mapping_summary.txt` for all samples

**Outputs**  
- Alignment rate plot: `RNAseq_alignment_rates_plot.pdf`

**Output directory**  
- `/datasets/RNAseq/hisat2file`

## StringTie_expression_quantification

**Description**
StringTie is used to estimate transcript abundances from RNA-seq alignments. It performs reference-guided quantification to generate gene and transcript expression levels for downstream analysis.

**Input**
All samples BAM files: `*[sampleID].bam`

**Parameters**

- `-p`: Number of processing threads (e.g., 5)
- `-e`: Quantification mode only (no novel transcript assembly)
- `-B`: Enable output for Ballgown downstream analysis
- `-G`: Reference annotation file in GTF format
- `-A`: Output file for gene abundances in tab-delimited format
- `-o`: Output file for assembled transcripts in GTF format

**Outputs**

- Gene abundance table: `gene_abund.tab`
- Assembled transcripts: `transcripts.gtf`

**Output directory**
`./datasets/RNAseq/stringtiefile/`

## RANseq_extract_expression_data

**Description**  
In this step, we use a custom script to extract the gene expression data (TPM) from the results of Stringtie.

**Input**
- Gene abundance table: `gene_abund.tab`

**Outputs**  
- Expression matrix: `gene_TPM.txt`


**Output directory**  
- `/datasets/RNAseq/stringtiefile`

## RNAseq_cluster_analysis

**Description**  
This step performs hierarchical clustering of samples based on their gene expression matrix to identify similarities and differences between samples.

**Input**
- Expression matrix: `gene_TPM.txt`

**Outputs**  
- Clustering dendrogram: `RNAseq_sample_clustering.pdf`

**Output directory**  
- `/datasets/RNAseq/stringtiefile`

## RANseq_PCA_analysis

**Description**  
Principal Component Analysis (PCA) is performed to visualize the overall structure of the data.We employ a custom code that will produce both 2D and 3D PCA maps.

**Input**  
- Expression matrix: `gene_TPM.txt`

**Outputs**  
- 2D PCA plot: `2DPCA_PC1_PC2.pdf` `2DPCA_PC1_PC3.pdf` `2DPCA_PC2_PC3.pdf`
- 2D PCA plot: `3DPCA.pdf.pdf`

**Output directory**  
- `/datasets/RNAseq/stringtiefile`

## RNAseq_differential_expression_analysis

### **Construct_samplelist**

**Description**  
DEseq2 requires defining sample names and experimental conditions in advance. This step creates a sample list file that defines the experimental design and the input file path.

**Input**
- Assembled transcripts: `transcripts.gtf`

**Outputs**  
- Sample list file: `samplelist.txt`

**Output directory**  
- `/datasets/RNAseq/stringtiefile`

### **Prepare_differential_expression_analysis**

**Description**  
This step prepares the data for differential expression analysis by normalizing the counts and creating the necessary data structures.

**Inputs**  
- Sample list file: `samplelist.txt`

**Outputs**  
- Normalized counts: `gene_count_matrix.count`

**Output directory**  
- `/datasets/RNAseq/stringtiefile`

### **Differential_expression_analysis**

**Description**  
This step identifies differentially expressed genes between conditions using statistical methods DESeq2 

**Inputs**  
- Normalized counts: `gene_count_matrix.count`

**Outputs**  
- Differential expression results: `Exper_vs_ctr_DEG.csv.txt`

**Output directory**  
- `/datasets/RNAseq/stringtiefile`