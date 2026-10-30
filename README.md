# RNA-seq Analysis | Differential Gene Expression Profiling

## Project Overview

This project presents an RNA sequencing (RNA-seq) analysis workflow for investigating differential gene expression using bioinformatics tools and R-based statistical analysis.

The workflow integrates raw sequencing quality control, read trimming, reference genome alignment, gene quantification, differential expression analysis, and biological data visualization.

The differential expression analysis evaluates transcriptional changes associated with dexamethasone treatment by comparing treated and untreated samples.

## Project Objectives

- Assess the quality of RNA-seq sequencing reads.
- Perform quality-based read trimming and preprocessing.
- Align sequencing reads to the human GRCh38 reference genome.
- Quantify gene-level expression from aligned reads.
- Identify differentially expressed genes between treated and untreated conditions.
- Visualize gene expression patterns and statistical results.
- Annotate genes and explore their biological functions.

## Bioinformatics Workflow

### 1. Quality Control — FastQC

Evaluated raw sequencing data using FastQC to assess sequencing quality and identify potential technical issues.

### 2. Read Trimming — Trimmomatic

Performed quality-based trimming to remove low-quality bases and improve sequencing read quality.

### 3. Read Alignment — HISAT2 and SAMtools

Aligned sequencing reads against the human GRCh38 reference genome using HISAT2 and generated sorted BAM files using SAMtools.

### 4. Gene Quantification — featureCounts

Quantified gene-level expression using featureCounts and generated a count matrix for downstream analysis.

### 5. Differential Gene Expression — DESeq2

Performed differential expression analysis in R using DESeq2, including:

- Low-count gene filtering
- Size factor estimation and normalization
- Dispersion estimation
- Negative binomial generalized linear model fitting
- Wald statistical testing
- Multiple-testing correction using the Benjamini–Hochberg method
- Identification of statistically significant differentially expressed genes

### 6. Quality Assessment and Visualization

Generated visualizations to examine gene expression patterns and differential expression results:

- Principal Component Analysis (PCA)
- Dispersion plot
- MA plot
- Volcano plot
- Heatmap of selected differentially expressed genes

### 7. Gene Annotation and Functional Analysis

Mapped Ensembl gene identifiers to gene symbols using the `org.Hs.eg.db` annotation package and explored gene functional annotations using DAVID.

## Tools and Technologies

| Category | Tools |
|---|---|
| Operating environment | Ubuntu Linux (WSL) |
| Quality control | FastQC |
| Read preprocessing | Trimmomatic |
| Alignment | HISAT2 |
| BAM processing | SAMtools |
| Gene quantification | featureCounts |
| Statistical analysis | R, DESeq2 |
| Data manipulation | dplyr |
| Visualization | ggplot2, ComplexHeatmap |
| Gene annotation | org.Hs.eg.db |
| Functional annotation | DAVID |
| Reference genome | Homo sapiens GRCh38 |

## Repository Structure

```text
rna-seq-analysis/
├── README.md
├── raw_data/
│   ├── script1.sh
│   ├── script2.sh
│   ├── script3.sh
│   ├── script4.sh
│   ├── output.csv
│   └── FastQC reports
└── DEG/
    ├── RNAseq_analysis.R
    ├── counts.csv
    ├── metadata2.csv
    ├── Finaltest.csv
    ├── best_genes.csv
    ├── pcaplot.png
    ├── Dispersionplot.png
    ├── MAplot.png
    ├── Volcanoplot.png
    └── heatmapplot.png
```

## Statistical Analysis

Differential gene expression was evaluated using DESeq2, which models RNA-seq count data using a negative binomial distribution.

Adjusted p-values were used to account for multiple hypothesis testing. Genes were examined using an FDR threshold of 0.05, with log2 fold changes used to assess the magnitude and direction of expression differences.

## Project Outputs

The analysis produced:

- Quality control reports
- Trimmed and aligned sequencing data
- Gene-level count matrices
- Differential gene expression result tables
- Annotated gene lists
- PCA, dispersion, MA, volcano, and heatmap visualizations
- Functional annotation results

## Data Availability

Large FASTQ, BAM, reference genome, and genome-index files are excluded from this repository to maintain a lightweight project structure.

The repository contains analysis scripts, selected processed data, result tables, and visualizations. Original datasets and third-party materials remain subject to their respective access and redistribution permissions.
