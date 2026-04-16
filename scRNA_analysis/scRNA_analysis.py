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
    print("Normalizing and identifying variable genes...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    #Saves raw data
    adata.raw= adata
    adata = adata[:, adata.var.highly_variable].copy()

    sc.tl.pca(adata, svd_solver='arpack')

    sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_pca')
    sc.tl.leiden(adata, resolution=0.5)
    sc.tl.umap(adata)
    print("Differential expression analysis...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

    marker_map = {
        "0": "CD4 T-cells", "1": "B-cells", '2':' CD14+ Monocytes', '3':'NK cells',
        '4':'CD8 T-cells', '5':'FCGR3A+ Monocytes', '6':'Dendritic cells', 
        '7':'platelets'
    }
    #edge case for when leiden clustering does not yield 8 clusters
    adata.obs['cell_type']=adata.obs['leiden'].map(lambda x: marker_map.get(x, f'cluster {x}'))


    print('Creating visuals')
    sc.pl.umap(adata, color=['cell_type', 'doublet_scores'],
        title=['Cell Type annotations', 'Doublet Scores'],
        frameon=False, legend_loc='on data', save='umap_celltype_doublets.png')  

    sc.pl.rank_genes_groups_dotplot(adata, n_genes=3, values_to_plot='logfoldchanges',
                                    min_logfoldchange=3, vmax=7, vmin=7, cmap ='bwr', save='_markers.png')

if __name__ == "__main__":
    processed_adata = run_scrna_pipeline()          