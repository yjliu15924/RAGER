# **RAGER table of contents**
1. [Quick start](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/README.md#quick-start)
2. [Preprocess RNAseq data](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Preprocess_RNAseq/RNAseq_analysis.md)
3. [Preprocess ATACseq data](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Preprocess_ATACseq/ATACseq_analysis.md) 
4. [Joint analysis](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Joint_analysis/Joint_analysis.md)
5. [Custom analysis](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Custom_analysis/Custom_analysis.md)

## **List of processes**
- [ATAC_quality_control](#)
  - [Trim_galore](#)
  - [Multiqc](#)
- [ATACseq_mapping_reads_to_genome](#)
  - [Align_ATACseq_reads_to_the_reference](#)
  - [Filter_multi-mapped_reads](#)
  - [Convert_sam_to_bam_format](#)
  - [Sort_bam_files_by_genome_coordinates](#)
  - [Mark_and_remove_duplicate_reads](#)
  - [Index_bam_files](#)
- [ATACseq_alignment_rates](#)
- [ATACseqQC](#)
- [Trans_gene_anno_to_bed](#)
- [Bam_to_bw](#)
- [Quantify TPM signals](#)
- [ATACseq_cluster_analysis](#)
- [ATACseq_PCA_analysis](#)
- [Merge_bam_file](#)
- [ATACseq_peak_calling](#)
- [ATACseq_peak_annotation](#)
- [ATACseq_peak_bar_plot](#)

## **Preprocessing**

## ATAC_quality_control

### **Trim_galore**

**Description**

Trim_galore is used to remove low quality bases and adapter sequences from the raw ATAC-seq reads. It performs quality control by trimming low-quality ends from reads in addition to adapter removal, which is particularly important for ATAC-seq data to ensure accurate mapping and peak calling.

**Inputs**
- Raw FASTQ files: `*_1.fq`, `*_2.fq` (for paired-end)
- Raw FASTQ files: `*_1.fq` (for single-end)

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
- `./datasets/ATACseq/quality_control_file/`

### **Multiqc**

**Description**

MultiQC aggregates results from various bioinformatics analyses into a single report. It collects the output from Trim_galore and other QC tools to generate comprehensive quality control metrics for ATAC-seq data.

**Input**
- FastqQC reports: `*fastqc.html`

**Outputs**  
- HTML report: `multiqc_report.html`
- Data directory: `multiqc_data/`

**Output directory**  
- `./datasets/ATACseq/quality_control_file/`

## ATACseq_mapping_reads_to_genome

### **Align_ATACseq_reads_to_the_reference**

**Description**  
Bowtie2 is used to align ATAC-seq reads to a reference genome. It is well-suited for short read alignment and is commonly used for ATAC-seq data to map open chromatin regions.

**Input**
- Trimmed FASTQ files: `*_1_val_1.fq`, `*_2_val_2.fq` (for paired-end)
- Trimmed FASTQ files: `*_1_trimmed.fq` (for single-end)

**Parameters**  
- `-x`: Reference genome index prefix
- `-1`: Forward/left reads file (R1)
- `-2`: Reverse/right reads file (R2)
- `-t`: Print wall-clock time taken by search phases
- `-q`: Input reads are in FASTQ format
- `-N`: Max mismatches in seed alignment (default: 1)
- `-L`: Length of seed substrings (default: 25)
- `-S`: Output alignment results in SAM format
- `-p`: Launch specified number of parallel search threads
- `--no-mixed`: Suppress unpaired alignments for paired reads
- `--no-discordant`: Suppress discordant alignments for paired reads

**Outputs**  
- SAM files: `*accepted_hits.sam`
- Alignment summary: `mapping_summary.txt`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

### **Filter_multi-mapped_reads**

**Description**  
This step filters out reads that align to multiple genomic locations, retaining only uniquely mapped reads to improve the accuracy of peak calling and other downstream analyses.

**Parameters**  
- `-H`: Extract header lines only
- `grep 'AS:'`: Find reads with alignment score (AS:) tag
- `grep -v 'XS:'`: Exclude reads with suboptimal alignment score (XS:) tag
**Outputs**  
- Filtered SAM files: `*accepted_hits_NHi1.sam`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

### **Convert_sam_to_bam_format**

**Description**

Samtools view converts SAM files to BAM format, which is a compressed binary version of the SAM format. This step reduces file size and prepares the data for downstream processing.

**Input**
- Filtered SAM files: `*accepted_hits_NHi1.sam`

**Parameters**
- `-S`: Input is SAM format
- `-b`: Output BAM format
- `-o`: Output file name
- `-@`: Number of threads

**Outputs**
- BAM files: `*accepted_hits_NHi1.bam`

**Output directory**
- `./datasets/ATACseq/bowtie2file/`

### **Sort.bam_files_by_genome_coordinates**

**Description**  
This step sorts the BAM files by genomic coordinates, which is required for many downstream analyses including peak calling and visualization.

**Input**
- BAM files: `*accepted_hits_NHi1.bam`

**Parameters**  
- `-o`: Output file name
- `-@`: Number of threads

**Outputs**  
- Sorted BAM files: `*accepted_hits_NHi1.sorted.bam`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

### **Mark_and_remove_duplicate_reads_using_Picard**

**Description**  
Picard's MarkDuplicates identifies and flags or removes duplicate reads that may have resulted from PCR amplification during library preparation, which is particularly important for ATAC-seq data to prevent bias in peak calling.

**Input**
- Sorted BAM files: `*accepted_hits_NHi1_sorted.bam`

**Parameters**  
- `-Xmx15g`: Set the maximum Java heap memory to 15 GB
- `I`: Input BAM file
- `O`: Output BAM file
- `METRICS_FILE`: File to write duplicate metrics
- `VALIDATION_STRINGENCY`: Controls validation strictness level(default: LENIENT)
- `REMOVE_DUPLICATES`: Whether to remove (true) or just mark (false) duplicates(default: true)
- `ASSUME_SORT_ORDER`: Specify the expected sort order of the input BAM file, skipping sort-order validation (default:coordinate)

**Outputs**  
- Deduplicated BAM files: `*[sampleID].bam`
- Duplicate metrics: `*.metricsFile`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

### **Index_bam_files_using_samtools**

**Description**  
This step creates an index for the BAM files, which allows for quick random access to the data, essential for visualization and downstream analyses.

**Input**
- Deduplicated BAM files: `*[sampleID].bam`

**Outputs**  
- BAM index files: `*[sampleID].bam.bai`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## ATACseq_alignment_rates

**Description**  
This step will present the alignment rates of all samples in a stacked bar chart, providing a quick overview of the mapping quality for ATAC-seq data.

**Input**
- Alignment summary: `mapping_summary.txt` for all samples

**Outputs**  
- Alignment rate plot: `ATACseq_alignment_rates_plot.pdf`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## ATACseqQC

**Description**  
ATACseqQC performs quality control analyses specific to ATAC-seq data, including fragment size distribution, TSS enrichment, and heatmap.

**Input**
- Indexed BAM files: `*[sampleID].bam` and `*[sampleID].bam.bai`

**Outputs**  
- Fragment size distribution plot: `fragmentSizeDistribution.pdf`
- Heatmap of signal around genomic features: `heatmap.pdf`
- Coverage curve plot: `coverage_curve_plot.pdf`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## Trans_gene_anno_to_bed

**Description**  
This step converts gene annotation from GTF format to BED format, which is required for many downstream analyses including peak annotation.

**Input**
- Gene annotation file: `gencode.v44.annotation.gtf`

**Outputs**  
- BED format gene annotation: `gencode.v44.annotation.bed`

**Output directory**  
- `./reference_annotation_file/geneanno/`

## Bam_to_bw

**Description**  
This step uses the **Deeptool** software bamCoverage tool to convert BAM files into bigWig format, which is more efficient for subsequent calculation of TPM signals and visualization of the genome browser.

**Input**
- Indexed BAM files: `*[sampleID].bam` and `*[sampleID].bam.bai`

**Parameters**
- `--binsize`: Bin size for signal calculation (default:10)
- `--ignoreDuplicates`: Skip PCR duplicate reads marked by Picard
- `--normalizeUsing` : Coverage normalization strategy (default:BPM)
- `-p`: Number of threads
- `-of`: Output format (bigwig)

**Outputs**  
- BigWig files: `*[sampleID].bw`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## Quantify TPM signals

**Description**  
This step uses the **Deeptool** software computeMatrix tool to quantifies the normalized signal (TPM - Tags Per Million) from ATAC-seq data, which allows for comparison between samples.

**Input**
- Indexed BAM files: `*[sampleID].bam`

**Parameters**
- `--reference-point`: Reference-point centered analysis mode
- `-b`: Distance upstream of reference point (bp) (default:1000)
- `-a`: Distance downstream of reference point (bp) (default:1000)
- `--binsize`: Bin size for signal calculation (default:10)

**Outputs**  
- Normalized signal table: `all_scaled.tab`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## ATACseq_cluster_analysis

**Description**  
This step performs hierarchical clustering of samples based on their ATAC-seq signals to identify similarities and differences between samples.

**Input**
- Normalized signal table: `all_scaled.tab`

**Outputs**  
- Clustering dendrogram: `sample_clustering.pdf`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## ATACseq_PCA_analysis

**Description**  
Principal Component Analysis (PCA) is performed to visualize the overall structure of the ATAC-seq data and identify major sources of variation.We employ a custom code that will produce both 2D and 3D PCA maps.

**Input**  
- Normalized signal table: `all_scaled.tab`

**Outputs**  
- 2D PCA plot: `2DPCA_PC1_PC2.pdf` `2DPCA_PC1_PC3.pdf` `2DPCA_PC2_PC3.pdf`
- 2D PCA plot: `3DPCA.pdf.pdf`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## Merge_bam_file

**Description**  
This step merges BAM files from replicates or related samples to increase the signal-to-noise ratio for peak calling.

**Input**
- Indexed BAM files: `*[sampleID].bam`

**Outputs**  
- Merged BAM files: `merged.bam`
- Merged BAM index files: `merged.bam.bai`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## ATACseq_peak_calling

**Description**  
MACS2 is used to call peaks from ATAC-seq data, which represent regions of open chromatin. This step identifies differential peaks between conditions.

**Input**
- Merged BAM files: `experiment.bam`, `ctr.bam`

**Parameters**
- `--treatment`: Treatment BAM file
- `--control`: Control BAM file
- `--genome`: Species of this data
- `--name`: Output file prefix
- `--outdir`: Output directory
- `-q`: Q-value cutoff for peak detection

**Outputs**  
- Peak files: `ATACexperiment_vs_ATACctr_summits.bed`, `ATACctr_vs_ATACexperiment_summits.bed`

**Output directory**  
- `./datasets/ATACseq/macs2file/`

## ATACseq_peak_annotation

**Description**  
This step annotates the called peaks with genomic features like promoters, enhancers, and gene bodies to provide functional context.

**Input**
- Peak files: `ATACexperiment_vs_ATACctr_summits.bed`, `ATACctr_vs_ATACexperiment_summits.bed`

**Outputs**  
- Annotated peak tables: `Experiment_vs_ctr_peak.csv`, `Ctr_vs_experiment_peak.csv`

**Output directory**  
- `./datasets/ATACseq/macs2file/`

## ATACseq_peak_bar_plot

**Description**  
This step creates bar plots to visualize the number and distribution of differential peaks between conditions.

**Input**
- Annotated peak tables: `GDC_vs_ctr_peak.csv`, `ctr_vs_GDC_peak.csv`

**Outputs**  
- Diverging bar plot: `peaks_diverging_bar.pdf`
- All peaks bar plot: `all_peaks_bar.pdf`

**Output directory**  
- `./datasets/ATACseq/macs2file/`

