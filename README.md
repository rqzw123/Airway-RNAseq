# Airway RNA-Seq Analysis: Transcriptional Response to Dexamethasone

This repository contains code (in R) to evaluate the effects of synthetic glucocorticoid exposure on airway smooth muscle cells.

## Project Overview

- **Objective**: The goal of the analysis is to assess the transcriptional response of human airway smooth muscle (ASM) cells to synthetic glucocorticoid exposure.

- **Experimental Design**: Bulk RNA-seq was performed on four human ASM cell lines under 2 conditions: [untreated]{.underline} and [treated with dexamethasone]{.underline}.

- **Data Source**: GEO accession number [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)

## Workflow

The core analysis pipeline is documented in Airway_Seq-Script.md and performs the following steps:

1.  **Data Importing**: Loading read counts, sample metadata, and gene annotations
2.  **Quality Control**: Filter low-count genes (min. 10), and evaluates sample distance using clustering heatmaps and principal component analysis (PCA).
3.  **Differential Expression**: Utilises the `DESeq2` framework, with a paired design formula for both cell line and dexamethasone treatment.
4.  **Visualisation**: Generates MA plots, an annotated heatmap of the top 30 differentially expressed genes, and a volcano plot.
5.  **Pathway Enrichment**: Performed Gene Set Enrichment Analysis (GSEA) on GO Biological Processes using `clusterProfiler` to identify functional changes.

## Prerequisites & Packages

To run the script successfully, ensure that the following R packages are installed:

- `DESeq2`

- `tidyverse`

- `clusterProfiler` and `enrichplot`

- `EnhancedVolcano`

- `pheatmap` and `RColorBrewer`

- `org.Hs.eg.db` and `AnnotationDbi`

## Key Conclusions 

- Exploratory analyses confirm robust data quality with high replicate concordance and no outlier samples.

- The primary source of variation is driven largely by treatment status (*dexamethasone vs untreated*).
