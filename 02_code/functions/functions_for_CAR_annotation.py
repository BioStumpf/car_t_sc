import anndata as ad
import numpy as np
import pandas as pd
import os
import regex as re
import matplotlib.pyplot as plt

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


#this function reads all files containing tsv dictionaries with information about the raw cellbarcodes
#and what barcode they have become within the actual anndata object
def read_cellbarcode_tsvs(folder):
    files = sorted(os.listdir(folder))
    col_names = ['raw', 'corrected']
    pool_corrections_list = []
    for file in files:
        full_path = os.path.join(folder, file)
        pool_correction = pd.read_csv(full_path, names = col_names, sep=" ", header = None)
        pool_correction[col_names[1]] = [brcd[:16] for brcd in pool_correction[col_names[1]]]
        pool_corrections_list.append(pool_correction)
    return pool_corrections_list


#this function will merge and format the dfs for containing the corrected barcode info
#together with the dfs containing the information about CAR-T counts, allowing to join the originial cellranger object based on corrected barcodes
def merge_CAR_counts_with_corrected_cellbarcodes(barcodes, counts):
    merged = pd.merge(barcodes, counts, left_on='raw', right_on='Cellbarcode', how='inner')
    merged.drop(columns=['Cellbarcode', 'raw'], inplace=True)
    merged.rename(columns={"corrected": "Cellbarcode"}, inplace=True)
    merged = merged.groupby("Cellbarcode", as_index=False).sum()
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
    categories = {'CD19_only': [], 'scFV_only': [], 'both': [], 'none': []} #, 'none': []
    for dataset in datasets:
        df_subset = df[df.dataset == dataset]
        cd19_only = ((df_subset['CD19_trunc'] > 0) & (df_subset['R11_ScFV'] < 1)).sum()
        scFV_only = ((df_subset['R11_ScFV'] > 0) & (df_subset['CD19_trunc'] < 1)).sum()
        both = ((df_subset['CD19_trunc'] > 0) & (df_subset['R11_ScFV'] > 0)).sum()
        none = ((df_subset['CD19_trunc'] == 0) & (df_subset['R11_ScFV'] == 0)).sum()
        
        categories['CD19_only'].append(cd19_only)
        categories['scFV_only'].append(scFV_only)
        categories['both'].append(both)
        categories['none'].append(none)
        indx = np.array(datasets) + 1
    return pd.DataFrame(categories, index=indx)

#write function extracting information about receptor count from adata object
def extract_receptor_count(adata):
    to_extract = ['condition', 'day', 'CD19_trunc', 'R11_ScFV']
    df = adata.obs[to_extract].copy()
    df_tmp = adata.obs[to_extract[:2]].copy()

    # Convert TCR counts to binary (1 if present, 0 if absent)
    df_tmp['CD19_trunc'] = ((df['CD19_trunc'] >= 1) & (df['R11_ScFV'] <  1)).astype(int)
    df_tmp['scFV'] = ((df['CD19_trunc'] <  1) & (df['R11_ScFV'] >= 1)).astype(int)
    df_tmp['both'] = ((df['CD19_trunc'] >= 1) & (df['R11_ScFV'] >= 1)).astype(int)
    df_tmp['none'] = ((df['CD19_trunc'] < 1) & (df['R11_ScFV'] < 1)).astype(int)

    # Count number of cells where TCRa or TCRb is present per condition and day
    df_counts = df_tmp.groupby(['condition', 'day'])[['CD19_trunc', 'scFV', 'both', 'none']].sum()
    # df_counts['TCRab_tot'] = df_counts.TCRa_only + df_counts.TCRb_only
    df_counts_reset = df_counts.reset_index()
    df_counts_reset['Total'] = df_counts_reset[['CD19_trunc', 'scFV', 'both', 'none']].sum(axis=1)
    return df_counts_reset

#plotting of receptor count per condition
def plot_receptor_count(df_counts_reset, xmax = 60, counts = 'Absolute Counts'):
    conditions = np.unique(df_counts_reset.condition)
    # TCRs = df_counts_reset.columns[2:6]

    fig, axs = plt.subplots(len(conditions), figsize=(8, 6))

    for idx, condition in enumerate(conditions):
        ax = axs[idx] if len(conditions) > 1 else axs  # Handles case with only one condition
        cond_subset = df_counts_reset[df_counts_reset.condition == condition]
        cond_subset = cond_subset.sort_values(by='day', ascending=False)

        colors = ['#56B4E9', '#E69F00', '#009E73']
        labels = ['CD19_only', 'scFV_only', 'both']
        categories = df_counts_reset.columns[2:5]

        left = np.zeros(len(cond_subset))  # Initialize left offsets as zeros

        for i, category in enumerate(categories):
            ax.barh(cond_subset.day, cond_subset[category], color=colors[i], label=labels[i], left=left)
            left += cond_subset[category].values  # Accumulate left offsets

        ax.set_xlim(0, xmax)
        ax.text(1.1, 0.5, f'Condition: {condition}', transform=ax.transAxes, ha='center', va='center', fontsize=12)

        for key, spine in ax.spines.items():
            spine.set_visible(False)

        if idx == len(conditions) - 1:
            ax.set_xlabel(counts)
        else:
            ax.set_xticks([])

    fig.legend(labels, loc='center left', bbox_to_anchor=(1.1, 0.81), title="")
    plt.show()