# Airway RNA-Seq Analysis: Transcriptional Response to Dexamethasone

This repository contains code (in R) to evaluate the effects of synthetic glucocorticoid exposure on airway smooth muscle cells.

## Project Overview

- **Objective**: The goal of the analysis is to assess the transcriptional response of human airway smooth muscle (ASM) cells to synthetic glucocorticoid exposure.

- **Experimental Design**: Bulk RNA-seq was performed on four human ASM cell lines under 2 conditions: *untreated* and treated with *dexamethasone*.

- **Data Source**: GEO accession number [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)

(note that the full GEO dataset contains other treatment arms such as `albuterol` and `dexamethasone+albuterol`, but dexamethasone will be the focus for this analysis.)

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

- Secondary clustering reflects donor-specific differences, confirming the necessity of a paired design formula (`~ cellLine + dexamethasone`) to effectively control for biological background variation.
- Applying log2 fold change shrinkage (`apeglm`) successfully mitigated statistical noise from low-count transcripts, ensuring that downstream pathway analysis was driven by true biological signals rather than mathematical artifacts.
- Over-Representation Analysis (ORA) revealed that dexamethasone strongly upregulates networks associated with cytoskeletal organisation, cell-substrate adhesion, and metabolic shifts (e.g., insulin response).
- Conversely, downregulated gene signatures were heavily enriched in broad developmental and extracellular matrix-related pathways, capturing the pervasive suppressive effects of glucocorticoids on cellular proliferation and off-target signaling.
