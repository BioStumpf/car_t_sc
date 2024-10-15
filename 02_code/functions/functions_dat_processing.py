import numpy as np
from scipy.stats import median_abs_deviation
from scipy.sparse import issparse
import os
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.sparse import csr_matrix

import anndata2ri
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import logging
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
ro.numpy2ri.activate()
anndata2ri.activate()


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
    conditions = []
    for (adata, adata_raw) in zip(adatas, adatas_raw):
        condition = list(adata.var[adata.var.feature_types == 'Antibody Capture'].gene_ids)
        #move HTO columns from variables to observables and demultiplex using hashsolo
        adata_new = hashing_columns(adata, rename_to=condition, rm_var_cols=True) 
        adata_raw_new = hashing_columns(adata_raw, rm_var_cols=True)
        #append to list of andata objects
        adatas_new.append(adata_new)
        adatas_raw_new.append(adata_raw_new)
        conditions.append(condition)
    return adatas_new, adatas_raw_new, conditions

def demultiplex(adatas, conditions):
    for (adata, condition) in zip(adatas, conditions):
        sce.pp.hashsolo(adata, condition)

#do the preclustering on each adata object, which is necessary for e.g. SoupX and scran normalization
# def pregroup(adatas: list, resolution = None):
#     for adata in adatas:
#         scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
#         adata.layers['log1p'] = sc.pp.log1p(scales_counts["X"], copy=True)
#         sc.pp.pca(adata, layer='log1p', n_comps=resolution)
#         sc.pp.neighbors(adata)
#         sc.tl.leiden(adata, key_added="groups", flavor="igraph", n_iterations=2)

def pregroup(adata, resolution = None):
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    adata.layers['log1p'] = sc.pp.log1p(scales_counts["X"], copy=True)
    sc.pp.pca(adata, layer='log1p', n_comps=resolution)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, key_added="groups", flavor="igraph", n_iterations=2)
    return adata

#define soupX function to apply soupX to all the pools
def cook_soup(adata, adata_raw):
    r_code = """
    make_soup <- function(data, data_raw, genes, cells, soupx_groups, empty_drops)
    {
        rownames(data) = genes
        colnames(data) = cells
        data <- as(data, "sparseMatrix")
        data_raw <- as(data_raw, "sparseMatrix")

        sc = SoupChannel(data_raw, data, calcSoupProfile = FALSE) #technically you do not even need the raw data since the empty droplets were computed manually, you may also use 2x data
        soupProf = data.frame(row.names = genes, est = rowSums(empty_drops)/sum(empty_drops), counts = rowSums(empty_drops))
        sc = setSoupProfile(sc, soupProf)
        sc = setClusters(sc, soupx_groups)
        sc  = autoEstCont(sc, doPlot=FALSE)
        if (is.null(sc$metaData$rho)) {
            stop("Contamination fractions were not calculated")
        }
        out = adjustCounts(sc, roundToInt = TRUE)
        return(out)
    }
    """
    groups = adata.obs['groups']
    empty_drops = find_empty_drops(adata_raw)
    cells = adata.obs_names
    genes = adata.var_names
    data = adata.X.T
    data_raw = adata_raw.X.T

    ro.r(r_code) #run the code above in the rpy2 environment to define the function in R.
    r_cook_soup = ro.globalenv['make_soup'] #take the function from the rpy2 so the R environment and make it globally available
    res = r_cook_soup(data, data_raw, genes, cells, groups, empty_drops) #apply the function
    adata.layers["counts"] = adata.X
    adata.layers["soupX_counts"] = res.T
    adata.X = adata.layers["soupX_counts"]
    return adata

#do the log1p/shifted log normalization
# def log1p_norm(adatas: list):
#     for adata in adatas:
#         scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
#         adata.layers['log1p'] = sc.pp.log1p(scales_counts["X"], copy=True)
def log1p_norm(adata):
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    adata.layers['log1p'] = sc.pp.log1p(scales_counts["X"], copy=True)
    return adata

#write function for scran normalization (optional)
def scran_norm(adata):
    r_code = """
    scran <- function(dat_mat, input_groups)
    {
        dat_mat <- SingleCellExperiment(list(counts=dat_mat))
        norm <- computeSumFactors(dat_mat, clusters=input_groups,  min.mean = 0.1, BPPARAM = MulticoreParam()) 
        sf <- sizeFactors(norm)
        return (sf)
    }
    """
    dat_mat = convert_to_csc(adata.layers['log1p'].T)
    input_groups = adata.obs['groups']

    ro.r(r_code)
    r_scran = ro.globalenv['scran']
    sf = r_scran(dat_mat, input_groups)

    adata.obs['size_factors'] = sf
    scran_norm = adata.X /adata.obs["size_factors"].values[:, None]
    adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran_norm))
    return adata


#write function to plot the normalization
def plot_normalization(adatas, norm_layer, title):
    for adata in adatas:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
        axes[0].set_title("Total counts")
        sns.histplot(
            adata.layers[norm_layer].sum(1), bins=100, kde=False, ax=axes[1]
        )
        axes[1].set_title(f"log1p with {title} estimated size factors")
        plt.show()



#write function for deviance feature selection
# def deviance_feature_selection(adata, n_top_genes):
#     r_code = """
#     get_deviance <- function(adata)
#     {
#         sce <- devianceFeatureSelection(adata, assay = "X")
#         return rowData(sce)$binomial_deviance
#     }
#     """
#     ro.r(r_code)
#     r_get_deviance = ro.globalenv['get_deviance']
#     binomial_deviance = r_get_deviance(adata)

#     idx = binomial_deviance.argsort()[-n_top_genes:]
#     mask = np.zeros(adata.var_names.shape, dtype=bool)
#     mask[idx] = True


#write function to plot deviance feature selection
