import anndata as ad
import numpy as np
import pandas as pd
import os
import regex as re

#this is to read the cellbarcodes that matched to the cart receptor
def read_and_merge_CAR_annotation(Folder):
    files = sorted(os.listdir(Folder))
    VDJ_GEX_list = []
    for regex in ['.*GEX', '.*VDJ']:
        r = re.compile(regex)
        file_subset = list(filter(r.match, files)) #r.match generates a function that takes any string as an input an returns true or false depending on whether the string could be matched
        pools = [pd.read_csv(os.path.join(Folder, file)) for file in file_subset]
        VDJ_GEX_list.append(pools)
        # print(file_subset)
    return VDJ_GEX_list

#merge also identical columns (both VDJ and GEX contain identical columns, however pandas does not properly merge them)
def merging_identical_cols(merged):
    for col in merged.columns:
        if col.endswith('_x'):  # Find columns with the `_x` suffix
            base_col = col[:-2]  # Remove `_x` to get the base name
            merged[base_col] = merged[col] + merged.get(f"{base_col}_y", 0)
            merged.drop([col, f"{base_col}_y"], axis=1, inplace=True)

#this is to merge the matched VDJ and GEX reads into one dataframe
def merge_VDJ_and_GEX(GEX, VDJ):
    if VDJ.empty:
        VDJ = pd.DataFrame(columns=['Cellbarcode', 'ReadIDs'])
    if GEX.empty:
        GEX = pd.DataFrame(columns=['Cellbarcode', 'ReadIDs'])
    VDJ = VDJ.drop(columns=['ReadIDs'])
    GEX = GEX.drop(columns=['ReadIDs'])
    merged = pd.merge(GEX, VDJ, on="Cellbarcode", how='outer')
    merged.fillna(0, inplace=True)
    merging_identical_cols(merged)
    return merged

#this is to add a column containing information as to whether a cell is a carT cell to the adata object
def annotate_mapped_cars(andata_obj, list_of_annotated_pools):
    if "dataset" in andata_obj.obs:
        new_obs_df = pd.DataFrame()
        for dataset in np.unique(andata_obj.obs.dataset):
            subset = andata_obj.obs[andata_obj.obs['dataset'] == dataset] #note: dataset is a string for some reason
            subset.index = [s[0:16] for s in subset.index]
            annotated_pool = list_of_annotated_pools[int(dataset)].copy()
            annotated_pool.index = annotated_pool.Cellbarcode
            annotated_pool.drop(columns=['Cellbarcode'], inplace=True)
            merged_obs = subset.join(annotated_pool, how='left')
            new_obs_df = pd.concat([new_obs_df, merged_obs])
    else:
        print('Column "dataset" non existent in andata object')
    andata_obj_cpy = andata_obj.copy()
    andata_obj_cpy.obs = new_obs_df
    return andata_obj_cpy

#this is to select only carTs
def isCAR(adata):
    return adata[(adata.obs['CD19_trunc'] > 0) | (adata.obs['R11_ScFV'] > 0)]

# Function to calculate how many cars have a mapped cd19, ScFV, neither or both to allow plotting
def calculate_categories(df, datasets):
    categories = {'CD19_only': [], 'R11_only': [], 'both': [], 'none': []} #, 'none': []
    for dataset in datasets:
        df_subset = df[df.dataset == dataset]
        cd19_only = ((df_subset['CD19_trunc'] > 0) & (df_subset['R11_ScFV'] < 1)).sum()
        r11_only = ((df_subset['R11_ScFV'] > 0) & (df_subset['CD19_trunc'] < 1)).sum()
        both = ((df_subset['CD19_trunc'] > 0) & (df_subset['R11_ScFV'] > 0)).sum()
        none = ((df_subset['CD19_trunc'] == 0) & (df_subset['R11_ScFV'] == 0)).sum()
        
        categories['CD19_only'].append(cd19_only)
        categories['R11_only'].append(r11_only)
        categories['both'].append(both)
        categories['none'].append(none)
        indx = np.array(datasets) + 1
    return pd.DataFrame(categories, index=indx)