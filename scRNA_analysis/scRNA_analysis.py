import scanpy as sc
import scrublet as scr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def run_scrna_pipeline():
    # Loads the data 
    print("Loading data...")
    adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique() # Removes duplicate gene names

    # Scrublet for doublet detection
    print("Running Scrublet for doublet detection...")
    scrub = scr.Scrublet(adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    adata.obs['doublet_scores'] = doublet_scores
    adata.obs['predicted_doublets'] = predicted_doublets

    # This filters out the predicted doublets from the dataset
    adata = adata[~adata.obs['predicted_doublets']].copy()

    # Quality Control 
    print("Performing quality control...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Pre filtered QC metrics Visualization
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show=False)
    plt.savefig('qc_metrics_violin.png')

    #Saves raw data
    adata.raw= adata
    adata = adata[:, adata.var.highly_variable].copy()