# Project: Enterococcus faecium E745 Genome & Transcriptome Analysis

This project reproduces the analysis of the multidrug-resistant clinical isolate *Enterococcus faecium* E745, focusing on its adaptation to human serum using long-read genome assembly, RNA-seq, and Tn-seq essential gene screening.

---

## Project Overview
- **Genome assembly**: PacBio long-read assembly with Canu, polished and annotated with Prokka.
- **RNA-seq differential expression**: Compare transcriptomes in BHI vs. human serum to identify genes up/downregulated under serum stress.
- **Tn-seq essential genes**: Use mariner transposon mutagenesis to identify genes required for survival in human serum.

---

## Repository Structure
- `raw_data/`: Sequencing data (PacBio, Illumina, RNA-seq, Tn-seq)
- `assembly/`: Genome assembly and quality control results
- `annotation/`: Prokka annotation files (GFF, FAA, etc.)
- `rna_seq/`: HISAT2 alignment, read counts, and DESeq2 differential expression results
- `tn_seq/`: Bowtie2 alignment, insertion site counts, and final matrix
- `wiki/`: GitHub Wiki with full analysis pipeline, commands, and figures

---

## Key Results
1.  **Genome**: 3.14 Mb complete genome, with a single 2.76 Mb contig (N50=2.76 Mb) and 3118 predicted CDS.
2.  **RNA-seq**: ~95% alignment rate across all samples; PCA shows clear separation between BHI and serum conditions.
3.  **Tn-seq**: Insertion counts at TA sites identified genes essential for serum survival.

---

## Access the Full Analysis
For step-by-step commands, visualizations, and methods, please visit the [GitHub Wiki](https://github.com/oriettazhan-collab/fake/wiki).
