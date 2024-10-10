import numpy as np
from scipy.stats import median_abs_deviation
from scipy.sparse import issparse
import os
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

#this function is to determine outliers. It can only be used after using sc.pp.calculate_qc_metrics().
#it takes a specific qc metric, like total count or total gene count and extracts the column corresponding to this metric form the adata object.
#then it computes wether the MAD (median absolute deviation) deviates nmads (e.g. 5) from the median. The output is an array of F/T.
def is_mad_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

#this function is also to determine outliers, but its related to finding the .xth quantile and filters based on this.
def is_qntl_outlier(adata, metric: str, qntl_threshhold: float):
    M = adata.obs[metric]
    lower_lim = np.quantile(M, qntl_threshhold)
    upper_lim = np.quantile(M, (1 - qntl_threshhold))
    outlier = (M < lower_lim) | (M > upper_lim)
    return outlier

#This function is only to be used when utilising hashtag oligos. It is meant to generate a viable input for the hashsolo demultiplexing function.
#The output of cellranger contains the HTOs as var_names withing the adata object, however hashsolo expects the counts for HTOs to be in the .obs dtaframe.
#This function takes the HTO feature types (assigned as Antibody Capture) and subsets the adata object to contain only this specific feature type.
#This is saved into a pandas dataframe, which can finally be joined with the .obs dataframe.
#Since no longer needed, the feature type is eliminated from the adtada.X matrix and only remains in the .obs df.
def hashing_columns(data, rename_to = None, rm_var_cols=False):   
    #copy the andata object
    dat_cpy = data.copy() 
    #find those andata.var rows matching antibody capture (which only works if you do not combine HTO with real antibody capture experiments)
    htos = dat_cpy.var['feature_types']  == 'Antibody Capture'
    #create a dataframe containing the counts for all HTOs
    htos_df = dat_cpy[:, htos].copy().to_df()
    #rename the dataframes columnames to whatever you find suitable (if you want to rename)
    if rename_to:
        #check if len of the list of names matches the len of the dataframe colnames
        if len(rename_to) == len(htos_df.columns):
            #if so, replace the old names with the new names
            htos_df.columns = rename_to
    #join the obs columns with the HTO columns
    dat_cpy.obs = dat_cpy.obs.join(htos_df)
    #if wanting to remove the var columns completely, remove them
    if rm_var_cols: dat_cpy = dat_cpy[:, ~htos].copy()
    return dat_cpy

#function for finding empty droplets in the raw data
def find_empty_drops(rawdata, range = [0, 100]):
    cell_sums = np.array(rawdata.X.sum(axis=1)).flatten()
    condition = (cell_sums > range[0]) & (cell_sums < range[1]) #why greater then 0 and not greater or equal? well you are aiming to compute an empty droplet profile, hence only droplets that do contain some genes are relevant
    empty_drops = rawdata[condition, :].X.T
    return empty_drops

#this function is to convert the andata.X.T sparse matrix into a csc format (if its below 32 bits)
#this is due to scran deconvolved normalization functions only being able to deal with csc type
def convert_to_csc(matrix):
    if issparse(matrix):
        if matrix.nnz > 2**31 - 1:
            matrix = matrix.tocoo()
        else:
            matrix = matrix.tocsc()
    return matrix

#######
# Now to the functions for processing all pools
#######

# This is to read all Pools within the given folder
# takes the path to the folder where all file subfolders are stored + where within these subfolders you can find the actual count matrix as input
def read_all_pools(path_to_folders: str, path_to_count_matrix_within_folders: str):
    #initialize a list for all andata objects
    adatas = []
    #read all folders in the given path
    subfolders = sorted(os.listdir(path_to_folders))
    #iterate through all subfolder, within each subfolder, go to the folder containing the count matrix
    for folder in subfolders:
        count_matrix_path = os.path.join(path_to_folders, folder, path_to_count_matrix_within_folders)
        #if the the given path is a directory, import the count matrix using scanpy and append to the list of andata objects
        if os.path.isdir(count_matrix_path):
            print(f'Reading from: {count_matrix_path}')
            adata = sc.read_10x_mtx(count_matrix_path, gex_only=False)
            adatas.append(adata)
        #if the given path is no directory, inform the user about it
        else:
            print(f'Could not read {count_matrix_path}')
    #finaly return the lenght of the given elements and the list of andata objects
    print(f'Read a total of {len(adatas)}')
    return adatas


def quality_control(adatas: list, method = 'mad'):
    #copy adatas list
    adatas_cp = []
    #iterate through copy and apply qc
    for i, adata in enumerate(adatas):
        print(f'computing qc for Pool {i+1}')
        #create a colum for mitochondrial genes and include it for qc
        adata.var['mt'] = adata.var_names.str.startswith('mt-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
        #determine mad, quantile or absolute outliers and filter them out
        if method.lower() == 'mad':
            adata.obs["outlier"] = (
                is_mad_outlier(adata, "total_counts", 5) #log1p_
                |  is_mad_outlier(adata, "n_genes_by_counts", 5) #log1p_
                |  is_mad_outlier(adata, "pct_counts_mt", 20) 
                # |  is_mad_outlier(adata, "pct_counts_in_top_20_genes", 5)
                |  (adata.obs["pct_counts_mt"] > 20)
            )
        elif method.lower() == 'qntl':
            adata.obs["outlier"] = (
                is_qntl_outlier(adata, "total_counts", .01)
                |  is_qntl_outlier(adata, "n_genes_by_counts", .01)
                |  is_qntl_outlier(adata, "pct_counts_mt", .01)
                |  (adata.obs["pct_counts_mt"] > 20)
                )
        elif method.lower() == 'abs':
            adata.obs['outlier'] = (
                (adata.obs['n_genes_by_counts'] > 6000)
                | (adata.obs['n_genes_by_counts'] < 200)
                | (adata.obs['pct_counts_mt'] > 20)
            )
        else:
            print('Filtering must be either mad, qntl or abs')
            return 0
        #finally subset according to previous filtering
        adata_subset = adata[(~adata.obs.outlier), :].copy()
        adatas_cp.append(adata_subset)
    return adatas_cp


#function for plotting the qc metrics
def plot_qc_metrics(adatas: list, adatas_qc: list):
    for i, (adata, adata_qc) in enumerate(zip(adatas, adatas_qc)):
        fig, axs = plt.subplots(1, 3, figsize=(16, 6))
        cells_before = adata.obs.shape[0]
        cells_after = adata_qc.obs.shape[0]
        fig.suptitle(f'QC Metrics Before (Cells: {cells_before}) and After (Cells: {cells_after}) for pool: {i+1}')

        metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
        before_qc = adata.obs[metrics].copy()
        after_qc = adata_qc.obs[metrics].copy()
        before_qc['Condition'] = 'Before QC'
        after_qc['Condition'] = 'After QC'
        combined_df = pd.concat([before_qc, after_qc])

        for j, metric in enumerate(metrics):
            # Melt the DataFrame to use seaborn's factor plot style
            melted_df = pd.melt(combined_df, id_vars='Condition', value_vars= metric, var_name='Metric')
            sns.violinplot(x='Metric', y='value', hue='Condition', data=melted_df, dodge=True, split=False, ax = axs[j], cut=0)
            axs[j].set_xlabel('')
            axs[j].set_ylabel('')
        axs[0].set_ylabel('Value')
        handles, labels = axs[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.9, 1))
        for ax in axs:
            ax.legend_.remove()
        plt.show()

#transfer HTOs from actual matrix to adata.obs (to use hashsolo later on in the workflow and not falsify quality control metrics)
def transfer_htos(adatas: list, adatas_raw: list):
    adatas_new = []
    adatas_raw_new = []
    for (adata, adata_raw) in zip(adatas, adatas_raw):
        conditions = list(adata.var[adata.var.feature_types == 'Antibody Capture'].gene_ids)
        adata_new = hashing_columns(adata, rename_to=conditions, rm_var_cols=True) 
        adata_raw_new = hashing_columns(adata_raw, rm_var_cols=True)
        adatas_new.append(adata_new)
        adatas_raw_new.append(adata_raw_new)
    return adatas_new, adatas_raw_new