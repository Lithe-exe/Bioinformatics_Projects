# End-to-End scRNA-seq Analysis Pipeline

## Overview
This program is meant to run an automated pipeline for analyzing single-cell RNA sequencing (scRNA-seq) data using **Scanpy**.The pipeline processes raw count matrices and performs quality control, doublet detection to cell-type annotation and differential expression analysis. Eventually it should be modular enough to run multiple different datasets automatically but is currently based only on one library.

## Key Features
*   **Doublet Detection:** Integrated `Scrublet` workflow to identify and remove multiplet artifacts.
*   **Quality Control (QC):** Automated filtering based on mitochondrial gene content, gene counts, and total UMI counts.
*   **Dimensionality Reduction:** Highly Variable Gene selection, PCA, and UMAP embedding.
*   **Clustering & Annotation:** Clustering using the Leiden algorithm and automated mapping to known PBMC cell types.
*   **Visuals:** Generation of UMAP plots and marker gene dotplots.

## Tech Stack
*   **Language:** Python 3.x
*   **Core Libraries:** `Scanpy`, `AnnData`
*   **Statistics/ML:** `Scrublet`, `NumPy`, `Pandas`
*   **Visualization:** `Matplotlib`, `Seaborn`

## Pipeline Workflow
1.  **Data Ingestion:** Load 10x Genomics PBMC 3k dataset.
2.  **Preprocessing:** 
    *   Scrublet for doublet removal.
    *   Filtering cells with high mitochondrial counts (>5%).
3.  **Normalization:** Total count normalization and log-transformation.
4.  **Feature Selection:** Identification of top highly variable genes.
5.  **Clustering:** Neighborhood graph construction followed by Leiden clustering.
6.  **Annotation:** Identifying cluster-specific markers (t-test) and assigning biological identities.

## Results
The pipeline outputs a final object containing the processed data as well as visualizations 