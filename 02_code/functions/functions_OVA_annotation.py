import anndata as ad
import numpy as np
import pandas as pd
import os
import regex as re
import matplotlib.pyplot as plt

#this is to ectract evalues from hmmer output files
def extract_evals(hmmer_output):
    with open(hmmer_output, 'r') as hmmer:
        evals = []
        for line in hmmer:
            if line.startswith('#'):
                continue
            # Split the line by tab and take the first entry == id
            line_elements = line.split()
            eval = float(line_elements[4])
            evals.append(eval)
        return evals
    
#this extracts all evals from a given folder
def extract_evals_from_folder(path_to_folder):
    pool_evals = {}
    for hmmer_file in sorted(os.listdir(path_to_folder)):
        full_path = os.path.join(path_to_folder, hmmer_file)
        if not os.path.isfile(full_path):
            continue
        evals = extract_evals(full_path)
        key = os.path.basename(full_path)
        pool_evals[key] = evals
    return pool_evals

#this is to read the cellbarcodes that matched to the cart receptor
def read_and_merge_annotation(Folder):
    files = sorted(os.listdir(Folder))
    VDJ_GEX_list = []
    for regex in ['.*GEX', '.*VDJ']:
        r = re.compile(regex)
        file_subset = list(filter(r.match, files)) #r.match generates a function that takes any string as an input an returns true or false depending on whether the string could be matched
        pools = [pd.read_csv(os.path.join(Folder, file)) for file in file_subset]
        VDJ_GEX_list.append(pools)
        # print(file_subset)
    return VDJ_GEX_list

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
def merge_OVA_counts_with_corrected_cellbarcodes(barcodes, counts):
    merged = pd.merge(barcodes, counts, left_on='raw', right_on='Cellbarcode', how='inner')
    merged.drop(columns=['Cellbarcode', 'raw'], inplace=True)
    merged.rename(columns={"corrected": "Cellbarcode"}, inplace=True)
    merged = merged.groupby("Cellbarcode", as_index=False).sum()
    return merged


#this is to add a column containing information as to whether a cell is a carT cell to the adata object
def annotate_mapped_OVA(andata_obj, list_of_annotated_pools):
    if "dataset" in andata_obj.obs:
        new_obs_df = pd.DataFrame()
        for dataset in np.unique(andata_obj.obs.dataset):
            subset = andata_obj.obs[andata_obj.obs['dataset'] == dataset].copy() #note: dataset is a string for some reason
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
def isOVA(adata):
    return adata[(adata.obs['TCRa'] > 0) | (adata.obs['TCRb'] > 0)]

#write function extracting information about receptor count from adata object
def extract_receptor_count(adata, to_extract, colnames, nr_groups):
    df = adata.obs[to_extract].copy()
    df_tmp = adata.obs[to_extract[:nr_groups]].copy()
    argc = len(to_extract)
    groups = to_extract[:nr_groups]
    targetA = to_extract[argc - 2]
    targetB = to_extract[argc - 1]

    # Convert TCR counts to binary (1 if present, 0 if absent)
    df_tmp[colnames[0]] = ((df[targetA] > 0) & (df[targetB] == 0)).astype(int)
    df_tmp[colnames[1]] = ((df[targetA] == 0) & (df[targetB] > 0)).astype(int)
    df_tmp[colnames[2]] = ((df[targetA] > 0) & (df[targetB] > 0)).astype(int)
    df_tmp[colnames[3]] = ((df[targetA] == 0)  & (df[targetB] == 0)).astype(int)

    # Count number of cells where TCRa or TCRb is present per condition and day
    df_counts = df_tmp.groupby(groups)[colnames].sum()
    # df_counts['TCRab_tot'] = df_counts.TCRa_only + df_counts.TCRb_only
    df_counts_reset = df_counts.reset_index()
    df_counts_reset['Total'] = df_counts_reset[colnames].sum(axis=1)
    return df_counts_reset

#plotting of receptor count per condition
def plot_receptor_counth(df_counts_reset, grouping, hue, xmax = 60, counts = 'Absolute Counts', figrsize = (8,6), label=slice(2, 5)):
    groups = np.unique(df_counts_reset[grouping])
    # TCRs = df_counts_reset.columns[2:6]
    fig, axs = plt.subplots(len(groups), figsize=figrsize)

    for idx, group in enumerate(groups):
        ax = axs[idx] if len(groups) > 1 else axs  # Handles case with only one condition
        group_subset = df_counts_reset[df_counts_reset[grouping] == group]
        # group_subset = group_subset.sort_values(by=hue, ascending=False) if sort==True else 
        group_subset = group_subset[::-1]

        colors = ['#E69F00', '#56B4E9', '#009E73', '#808080']
        categories = df_counts_reset.columns[label]

        left = np.zeros(len(group_subset))  # Initialize left offsets as zeros
        for i, category in enumerate(categories):
            ax.barh(group_subset[hue], group_subset[category], color=colors[i], label=categories[i], left=left)
            left += group_subset[category].values  # Accumulate left offsets

        ax.set_xlim(0, xmax)
        ax.text(1.1, 0.5, group, transform=ax.transAxes, ha='center', va='center', fontsize=12)

        for key, spine in ax.spines.items():
            spine.set_visible(False)

        if idx == len(groups) - 1:
            ax.set_xlabel(counts)
        else:
            ax.set_xticks([])

    fig.legend(categories, loc='center left', bbox_to_anchor=(1.03, 0.81), title="")
    plt.show()

def plot_receptor_countv(df_counts_reset, grouping, hue, ymax = 60, counts = 'Absolute Counts', figrsize = (8,6), label=slice(2, 5)):
    groups = np.unique(df_counts_reset[grouping])
    # TCRs = df_counts_reset.columns[2:6]
    fig, axs = plt.subplots(ncols=len(groups), figsize=figrsize)

    for idx, group in enumerate(groups):
        ax = axs[idx] if len(groups) > 1 else axs  # Handles case with only one condition
        group_subset = df_counts_reset[df_counts_reset[grouping] == group]
        # group_subset = group_subset.sort_values(by=hue, ascending=False) if sort==True else 
        # group_subset = group_subset[::-1]

        # colors = ['#E69F00', '#56B4E9', '#009E73', ]
        colors = ['#E69F00', '#56B4E9', '#009E73', '#808080']  # Add a fourth color

        # labels = ['TCRa only', 'TCRb only', 'TCRb + TCRa']
        # labels = df_counts_reset.columns[2:5]
        categories = df_counts_reset.columns[label]

        bottom = np.zeros(len(group_subset))  # Initialize left offsets as zeros

        for i, category in enumerate(categories):
            # ax.barh(group_subset[hue], group_subset[category], color=colors[i], label=categories[i], left=left)
            ax.bar(group_subset[hue], group_subset[category], color=colors[i], label=categories[i], bottom=bottom)
            bottom += group_subset[category].values  # Accumulate left offsets

        ax.set_ylim(0, ymax)
        ax.text(0.4, -0.1, group, transform=ax.transAxes, ha='center', va='center', fontsize=12)
        ax.set_xticks(group_subset[hue])  # Set tick positions
        ax.set_xticklabels(group_subset[hue], rotation=45, ha="right")  # Rotate labels

        for key, spine in ax.spines.items():
            spine.set_visible(False)

        if idx == 0:
            ax.set_ylabel(counts)
        else:
            ax.set_yticks([])
            

    fig.legend(categories, loc='center left', bbox_to_anchor=(0.8, 0.9), title="")
    plt.show()

#to generate a .obs column for each replicate, usefull for statistical tesing
def extract_replicate(adata, column):
    rep_info = adata.obs[column]
    r = re.compile('(C|DM|P)(\d)')
    replicates = [int(re.search(r, cell).group(2)) if re.search(r, cell) else 1 for cell in rep_info]
    adata.obs['replicate'] = replicates
