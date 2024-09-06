import numpy as np
from scipy.stats import median_abs_deviation
from scipy.sparse import csr_matrix, issparse

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
    condition = (cell_sums > range[0]) & (cell_sums < range[1])
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
