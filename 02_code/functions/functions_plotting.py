############
#Functions for Plotting
############

from matplotlib import pyplot as plt
import anndata as ad
import numpy as np
import pandas as pd
import os
import regex as re

#This function extracts a specific day and condition from the adata objecs classification
def extract_DayCond(days, conditions, classification):
    day = ''.join(filter(lambda day: day in classification, days))
    condition = ''.join(filter(lambda condition: condition in classification, conditions))
    return day, condition
#Tis function applies the day and condition extraction for all given days and conditions
def extract_DaysConds(days, conditions, classifications):
    return [extract_DayCond(days, conditions, classification) for classification in classifications]
#this function extracts all given days and conditions from the classification column and adds them to a distinct .obs column of the adata object
def add_DaysConds_to_adata(adata, days, conditions):
    classifications = adata.obs.Classification.values #specific_class, Classification
    res = extract_DaysConds(days, conditions, classifications)
    cond_days, cond_conds = zip(*res)
    adata.obs['condition'] = cond_conds
    adata.obs['day'] = cond_days

#This function creates a stacked barplot for the adata object, based on 
def stacked_barplot(adata, cell_subtype, obs_column):
    conditions = np.unique(adata.obs.condition.values)
    days = sorted(np.unique(adata.obs.day), key=lambda x: int(x), reverse=True)
    colors = ['#FF69B4', '#4B0082', '#228B22', '#FFD700', '#4682B4', '#FF4500', '#8B4513']

    # Filter for only selected subtype
    subtype_only = list(map(lambda subtype: cell_subtype in subtype, adata.obs[obs_column]))
    adata = adata[subtype_only, :].copy()
    unq_celltypes = np.unique(adata.obs[obs_column])

    # Set up subplots for each condition
    fig, axs = plt.subplots(len(conditions), figsize=(8, 6)) #, sharey=True

    # Plot each condition in a separate subplot
    for idx, condition in enumerate(conditions):
        ax = axs[idx] if len(conditions) > 1 else axs  # Handles case with only one condition
        subset = adata[adata.obs.condition == condition, :]
        
        # Initialize dictionary for each celltype's counts by day
        celltype_counts = {celltype: [] for celltype in unq_celltypes}
        
        for celltype in unq_celltypes:
            celltype_subset = subset[subset.obs[obs_column] == celltype, :]
            for day in days:
                day_subset = celltype_subset[celltype_subset.obs.day == day, :]
                day_count = len(day_subset)  # Count the entries for each day and cell type
                celltype_counts[celltype].append(day_count)

        # Normalize cell type counts by day
        norm_factors = [sum(day_counts) for day_counts in zip(*celltype_counts.values())]
        norm_factors = [nf if nf > 0 else 1 for nf in norm_factors] 
        for key, values in celltype_counts.items():
            celltype_counts[key] = (np.array(values) / norm_factors) * 100

        # Stacked bar plot for the current condition
        left = np.zeros(len(days))  # Reset the 'left' offset for each condition
        for i, (celltype, celltype_count) in enumerate(celltype_counts.items()):
            p = ax.barh(days, celltype_count, 0.8, label=celltype, left=left, color=colors[i % len(colors)])
            left += celltype_count
            # ax.bar_label(p, label_type='center', labels=[f'{val:.1f}' for val in celltype_count])
        ax.text(1.05, 0.5, f'Condition: {condition}', transform=ax.transAxes, ha='center', va='center', fontsize=12) #rotation=270,
        if idx == len(conditions) - 1:
            ax.set_xlabel('Percentage')
        else:
            ax.set_xticks([])

        for key, spine in ax.spines.items():
            spine.set_visible(False)

    fig.legend(unq_celltypes, loc='center left', bbox_to_anchor=(1.005, 0.75), title="Cell Types")
    # plt.subplots_adjust(left=0.05, right=0.8, wspace=0.15)
    fig.suptitle(f'Fraction of {cell_subtype} T cells by Condition')
    # plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout to fit suptitle
    plt.show()

#this is to read the cellbarcodes that matched to the cart receptor
def read_and_merge_CAR_annotation(Folder):
    files = sorted(os.listdir(Folder))
    VDJ_GEX_list = []
    for regex in ['.*GEX', '.*VDJ']:
        r = re.compile(regex)
        file_subset = list(filter(r.match, files))
        pools = [pd.read_csv(os.path.join(Folder, file)) for file in file_subset]
        VDJ_GEX_list.append(pools)
        # print(file_subset)
    return VDJ_GEX_list

#this is to merge the matched VDJ and GEX reads into one dataframe
def merge_VDJ_and_GEX(GEX, VDJ):
    VDJ = VDJ.drop(columns=['ReadIDs'])
    GEX = GEX.drop(columns=['ReadIDs'])
    merged = pd.merge(GEX, VDJ, on="Cellbarcode", how='outer')
    merged.fillna(0, inplace=True)
    return merged

#this is to add a column containing information as to whether a cell is a carT cell to the adata object
def annotate_mapped_cars(andata_obj, list_of_annotated_pools):
    if "dataset" in andata_obj.obs:
        intersection_all = []
        for dataset in np.unique(andata_obj.obs.dataset):
            subset = andata_obj.obs[andata_obj.obs['dataset'] == dataset] #note: dataset is a string for some reason
            # print(subset.index)
            subset_barcodes = [s[:16] for s in subset.index]
            # print(subset_barcodes)
            matches = set(list_of_annotated_pools[int(dataset)].Cellbarcode)
            # print(matches)
            intersection_subset = [barcode in matches for barcode in subset_barcodes]
            # print(sum(intersection_subset))
            intersection_all.extend(intersection_subset)
    else:
        print('Column "dataset" non existent in andata object')
    andata_obj_cpy = andata_obj.copy()
    andata_obj_cpy.obs['isCAR'] = intersection_all
    return andata_obj_cpy