# Genome Analysis Project – *Enterococcus faecium* E745

This repository contains the analysis files and supplementary material for the re-analysis of genomic and transcriptomic data from Zhang et al. (2017).

The project includes genome assembly and evaluation, genome annotation, genome comparison, RNA-seq quality control and trimming, RNA-seq mapping, read counting, differential expression analysis, and additional Tn-seq analysis.

## Analysis Workflow

The complete analysis workflow, including methods, results, figures and discussion, is presented in the project Wiki:

[Full Analysis Workflow](https://github.com/oriettazhan-collab/fake/wiki/Full-Analysis-Workflow)

## Code

Most short analysis commands are shown directly in the relevant Methods sections of the Wiki.

The complete R script used for the DESeq2 differential expression analysis, PCA and volcano plot generation is stored separately because of its length:

[`code/DESeq2_analysis.R`](https://github.com/oriettazhan-collab/fake/blob/main/code/DESeq2_analysis.R)

## Data and Output Files

Important supplementary files available in this repository include:

- `final_count_matrix_for_R.txt` – RNA-seq count matrix used as input for DESeq2
- `TN_SEQ_FINAL_MATRIX.txt` – Tn-seq insertion-site count matrix
- `ERR1797972_1_trim_fastqc.html` – representative post-trimming FastQC report

## Main Tools

The main tools used in the workflow were:

- Canu – genome assembly
- QUAST – assembly quality evaluation
- Prokka – genome annotation
- MUMmer – genome comparison and synteny analysis
- FastQC and fastp – RNA-seq quality control and preprocessing
- HISAT2 and SAMtools – RNA-seq alignment and BAM processing
- HTSeq-count – gene-level read counting
- DESeq2 – differential expression analysis
- Bowtie2 – Tn-seq alignment
