# Title: scrublet_estimate_doublecell.py
# Date: 2024-04-22
# Attention: how to rationally get a multi-matrix anndata including FilterMatrix, SpliceMatrix and UnspliceMatrix.

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
from matplotlib.pyplot import savefig
from pathlib import Path
import shutil
import gzip
import os
import sys
import scrublet
import leidenalg
import argparse

#get outdoor parameter
species = sys.argv[1]
input_mingenes = int(sys.argv[2])
input_mincells = int(sys.argv[3])
group_key = sys.argv[4]
matrix_filepath_txt = sys.argv[5]
splice_filepath_txt = sys.argv[6]
unsplice_filepath_txt = sys.argv[7]
samples_txt_path = sys.argv[8]
print(len(sys.argv))

try:
    txt_rho = sys.argv[9]
except IndexError:
    print("Warning: txt_rho parameter not provided. Using un-soupx deal.")
    
try:
    if len(sys.argv) in [11, 12]:
        mito_genes = sys.argv[10] if len(sys.argv) == 12 else sys.argv[9]
        mito_threshold = float(sys.argv[11]) if len(sys.argv) == 12 else float(sys.argv[10])
    else:
        raise ValueError("Invalid number of arguments")
except (IndexError, ValueError) as e:
    print(f"Warning: {e}. Using default values for mito_genes and mito_threshold.")
    mito_genes = "None_mito_genes.csv"
    mito_threshold = 0.05
print(f"mito_genes: {mito_genes}")
print(f"mito_threshold: {mito_threshold}")

#from text transform to array
with open(matrix_filepath_txt, 'r') as file:
    matrix_files = file.read().strip().split(',')
with open(splice_filepath_txt, 'r') as file:
    splice_files = file.read().strip().split(',')
with open(unsplice_filepath_txt, 'r') as file:
    unsplice_files = file.read().strip().split(',')
with open(samples_txt_path, 'r') as filen:
    sample_names = filen.read().strip().split(',')

#check soupxmatrix
if len(matrix_files) == 1 and len(matrix_files) != len(sample_names) and matrix_files[0].endswith("result"):
    base_path = matrix_files[0]
    matrix_files = [f"{base_path}/{sample}" for sample in sample_names]
    #matrix_files = new_matrix_files
else:
    matrix_files = matrix_files

print(matrix_files)

def copy_and_process(matrixfile, featuresfile, barcodesfile, target_folder):
    original_dir = os.getcwd()
    os.chdir(target_folder)
    shutil.copy(matrixfile, "matrix.mtx.gz")
    shutil.copy(featuresfile, "features.tsv.gz")
    shutil.copy(barcodesfile, "barcodes.tsv.gz")
    with gzip.open('matrix.mtx.gz', 'rb') as g_file1, open("matrix.mtx", "wb") as f_out:
        f_out.write(g_file1.read())
    with gzip.open('features.tsv.gz', 'rb') as g_file2, open("features.tsv", "wb") as f_out:
        f_out.write(g_file2.read())
    with gzip.open('barcodes.tsv.gz', 'rb') as g_file3, open("barcodes.tsv", "wb") as f_out:
        f_out.write(g_file3.read())
    with open('features.tsv', 'r') as f_in, open('genes.tsv', 'w') as f_out:
        for line in f_in:
            f_out.write(line.strip() + '\t' + line.strip() + '\n')
    os.chdir(original_dir)

def process_missing_genes(adata, all_genes):
    adata = adata[:, adata.var['gene_symbols'].isin(all_genes)].copy()
    missing_genes = np.setdiff1d(all_genes, adata.var['gene_symbols'])
    if len(missing_genes) > 0:
        missing_data = ad.AnnData(
            X=np.zeros((adata.n_obs, len(missing_genes))),
            var=pd.DataFrame(index=missing_genes, data=missing_genes, columns=['gene_symbols']),
            obs=adata.obs.copy()
        )
        missing_data.var.set_index('gene_symbols', inplace=True)
        adata = ad.concat([adata, missing_data], axis=1, join='outer')
    adata.var.set_index('gene_symbols', inplace=True)
    return adata

def run_scrublet_view(species, input_mingenes, input_mincells, group_key, sample_names, trans_matrix_list, trans_splice_list, trans_unsplice_list, mito_genes, mito_threshold):
    adatas = {}
    for i in range(len(sample_names)): 
        key = sample_names[i]
        adata_filter = sc.read_10x_mtx(trans_matrix_list[i], var_names='gene_ids')
        adata_splice = sc.read_10x_mtx(trans_splice_list[i], var_names='gene_ids')
        adata_unsplice = sc.read_10x_mtx(trans_unsplice_list[i], var_names='gene_ids')
        all_genes = adata_filter.var_names
        adata_filter = process_missing_genes(adata_filter, all_genes)
        adata_splice = process_missing_genes(adata_splice, all_genes)
        adata_unsplice = process_missing_genes(adata_unsplice, all_genes)
        adata = adata_filter.copy()
        adata.layers['splice'] = adata_splice.X
        adata.layers['unsplice'] = adata_unsplice.X
        # rename cells to include sample key
        adata.obs_names = [f"{cell_name}_{key}" for cell_name in adata.obs_names]
        # store in dictionary
        adatas[key] = adata
    adata = ad.concat(adatas, label=group_key, join="inner")  # the variable [join] is key, could select "outer" or "inner".
    print(adata.obs[group_key].value_counts())

    # Set parameters for figures
    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80, facecolor='white')
    
    # Check mitogenes and filter
    if os.path.exists(mito_genes):
        mt_genes = pd.read_csv(mito_genes, header=None, names=["gene_name"])
        mt_genes_list = mt_genes["gene_name"].tolist()
        print(mt_genes_list[:10])
        adata.var["mt"] = adata.var_names.isin(mt_genes)
        print("calculate mt genes")
        sc.pp.calculate_qc_metrics(adata,qc_vars=["mt"],inplace=True,log1p=True)
        sc.pl.violin(adata,["n_genes_by_counts", "total_counts", "pct_counts_mt"],jitter=0.4,multi_panel=True,save="_mitogene.pdf")
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save="_mitogenes.pdf")
        adata = adata[adata.obs.pct_counts_mt < mito_threshold].copy()
    else:
        print("mitochondrial list not exist")
        sc.pp.calculate_qc_metrics(adata, inplace=True, log1p=True)
    sns.jointplot(data=adata.obs, x="log1p_total_counts", y="log1p_n_genes_by_counts", kind="hex")
    savefig("qc.pdf")

    # Pre-process, control quality, and Scrublet
    sc.pp.filter_cells(adata, min_genes=input_mingenes)
    sc.pp.filter_genes(adata, min_cells=input_mincells)
    sc.external.pp.scrublet(adata, batch_key=group_key)

    adata.layers["counts"] = adata.X.copy()

    # Visualization
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=group_key)
    sc.tl.pca(adata)
    sc.pl.pca(adata, color=[group_key, "doublet_score", "log1p_total_counts", "log1p_n_genes_by_counts"],
              dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)], ncols=2, size=2, save='_potentially_undesired_features.pdf')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=group_key, size=2, save="_batch.pdf")
    sc.tl.leiden(adata, resolution=1)
    adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype('category')
    sc.pl.umap(adata, color=["leiden", "log1p_n_genes_by_counts", "predicted_doublet", "doublet_score"], ncols=2, save="_quality.pdf")
    for res in [0.02, 0.2, 0.5, 0.8, 1.0, 1.3, 1.6, 2.0]:
        sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res)
    sc.pl.umap(adata, color=["leiden_res_0.02", "leiden_res_0.20", "leiden_res_0.50", "leiden_res_0.80", "leiden_res_1.00", "leiden_res_1.30", "leiden_res_1.60", "leiden_res_2.00"], legend_loc="on data", save="_leiden_clus.pdf")
    sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")
    sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=5, save="marker.pdf")
    marker = sc.get.rank_genes_groups_df(adata, group=None)
    marker.to_csv("leiden_res_0.50.markers.csv")
    with open('summary.txt', 'w') as f:
        f.write(species + ' data summary' + '\n')
        f.write('Total cells: ' + str(adata.n_obs) + '\n')
        f.write('Total genes: ' + str(adata.n_vars) + '\n')
        f.write('Average genes per cell: ' + str(adata.obs['n_genes'].mean()) + '\n')
        f.write('Median genes per cell: ' + str(adata.obs['n_genes'].median()) + '\n')
        f.write('Average counts per cell: ' + str(adata.obs['total_counts'].mean()) + '\n')
        f.write('Median counts per cell: ' + str(adata.obs['total_counts'].median()) + '\n')
    adata.write_h5ad(filename=species + '.h5ad', compression="gzip")

trans_matrix_list = []
trans_splice_list = []
trans_unsplice_list = []

process_types = [
    ("filter", matrix_files),
    ("splice", splice_files),
    ("unsplice", unsplice_files)
]

if len(matrix_files) > 0:
    for i in range(len(sample_names)):
        sample = sample_names[i]
        for process_name, file_list in process_types:
            directory_path = Path(f"./{sample}/{process_name}")
            directory_path.mkdir(parents=True, exist_ok=True)
            folder_path = os.path.abspath(directory_path)
            if process_name == "filter":
                trans_matrix_list.append(folder_path)
                matrixfile = file_list[i] + '/matrix.mtx.gz'
                featuresfile = file_list[i] + '/features.tsv.gz'
                barcodesfile = file_list[i] + '/barcodes.tsv.gz'
            elif process_name == "splice":
                trans_splice_list.append(folder_path)
                matrixfile = file_list[i] + '/matrix.mtx.gz'
                featuresfile = file_list[i] + '/features.tsv.gz'
                barcodesfile = file_list[i] + '/barcodes.tsv.gz'
            elif process_name == "unsplice":
                trans_unsplice_list.append(folder_path)
                matrixfile = file_list[i] + '/unspliced.mtx.gz'
                featuresfile = file_list[i] + '/features.tsv.gz'
                barcodesfile = file_list[i] + '/barcodes.tsv.gz'
            copy_and_process(matrixfile, featuresfile, barcodesfile, folder_path)
    run_scrublet_view(species, input_mingenes, input_mincells, group_key, sample_names, trans_matrix_list, trans_splice_list, trans_unsplice_list, mito_genes, mito_threshold)
else:
    print("No samples to process")
    with open('summary.txt', 'w') as f:
        f.write(species + ' data summary' + '\n')
        f.write('No samples to process' + '\n')
