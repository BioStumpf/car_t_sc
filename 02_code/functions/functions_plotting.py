############
#Functions for Plotting
############

from matplotlib import pyplot as plt
import anndata as ad
import numpy as np
import pandas as pd
import os 
import warnings
import seaborn as sns
from adjustText import adjust_text
import matplotlib.patheffects as PathEffects
import random
import io
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import sys

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
# def stacked_barplot(adata, cell_subtype, obs_column):
#     conditions = np.unique(adata.obs.condition.values)
#     days = sorted(np.unique(adata.obs.day), key=lambda x: int(x), reverse=True)
#     colors = ['#FF69B4', '#4B0082', '#228B22', '#FFD700', '#4682B4', '#FF4500', '#8B4513']

#     # Filter for only selected subtype
#     subtype_only = list(map(lambda subtype: cell_subtype in subtype, adata.obs[obs_column]))
#     adata = adata[subtype_only, :].copy()
#     unq_celltypes = np.unique(adata.obs[obs_column])

#     # Set up subplots for each condition
#     fig, axs = plt.subplots(len(conditions), figsize=(8, 6)) #, sharey=True

#     # Plot each condition in a separate subplot
#     for idx, condition in enumerate(conditions):
#         ax = axs[idx] if len(conditions) > 1 else axs  # Handles case with only one condition
#         subset = adata[adata.obs.condition == condition, :]
        
#         # Initialize dictionary for each celltype's counts by day
#         celltype_counts = {celltype: [] for celltype in unq_celltypes}
        
#         for celltype in unq_celltypes:
#             celltype_subset = subset[subset.obs[obs_column] == celltype, :]
#             for day in days:
#                 day_subset = celltype_subset[celltype_subset.obs.day == day, :]
#                 day_count = len(day_subset)  # Count the entries for each day and cell type
#                 celltype_counts[celltype].append(day_count)

#         # Normalize cell type counts by day
#         norm_factors = [sum(day_counts) for day_counts in zip(*celltype_counts.values())]
#         norm_factors = [nf if nf > 0 else 1 for nf in norm_factors] 
#         for key, values in celltype_counts.items():
#             celltype_counts[key] = (np.array(values) / norm_factors) * 100

#         # Stacked bar plot for the current condition
#         left = np.zeros(len(days))  # Reset the 'left' offset for each condition
#         for i, (celltype, celltype_count) in enumerate(celltype_counts.items()):
#             p = ax.barh(days, celltype_count, 0.8, label=celltype, left=left, color=colors[i % len(colors)])
#             left += celltype_count
#             # ax.bar_label(p, label_type='center', labels=[f'{val:.1f}' for val in celltype_count])
#         ax.text(1.05, 0.5, f'Condition: {condition}', transform=ax.transAxes, ha='center', va='center', fontsize=12) #rotation=270,
#         if idx == len(conditions) - 1:
#             ax.set_xlabel('Percentage')
#         else:
#             ax.set_xticks([])

#         for key, spine in ax.spines.items():
#             spine.set_visible(False)

#     fig.legend(unq_celltypes, loc='center left', bbox_to_anchor=(1.005, 0.75), title="Cell Types")
#     # plt.subplots_adjust(left=0.05, right=0.8, wspace=0.15)
#     fig.suptitle(f'Fraction of {cell_subtype} T cells by Condition')
#     # plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout to fit suptitle
#     plt.show()


def stacked_barplot(adata, obs_column, xmax, xlabel, colors, common_cell_subtype = None, norm = True, save_dir=None):
    conditions = np.unique(adata.obs.condition.values)
    days = sorted(np.unique(adata.obs.day), key=lambda x: int(x), reverse=True)
    #to avoid some annoying infos
    plt.set_loglevel('WARNING') 

    # colors = ['#FF69B4', '#4B0082', '#228B22', '#FFD700', '#4682B4', '#FF4500', '#8B4513']
    # Filter for only selected common subtype 
    if common_cell_subtype:
        subtype_only = list(map(lambda subtype: common_cell_subtype in subtype, adata.obs[obs_column]))
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

        if norm:
            # Normalize cell type counts by day
            norm_factors = [sum(day_counts) for day_counts in zip(*celltype_counts.values())]
            norm_factors = [nf if nf > 0 else 1 for nf in norm_factors] 
            for key, values in celltype_counts.items():
                celltype_counts[key] = (np.array(values) / norm_factors) * 100

        # Stacked bar plot for the current condition
        left = np.zeros(len(days))  # Reset the 'left' offset for each condition
        for i, (celltype, celltype_count) in enumerate(celltype_counts.items()):
            p = ax.barh(days, celltype_count, 0.8, label=celltype, left=left, color=colors[celltype]) #colors[i % len(colors)]   colors[celltype]
            left += celltype_count
            # ax.bar_label(p, label_type='center', labels=[f'{val:.1f}' for val in celltype_count])
        ax.text(1.05, 0.5, condition, transform=ax.transAxes, ha='center', va='center', fontsize=12) #rotation=270,
        if idx == len(conditions) - 1:
            ax.set_xlabel(xlabel)
        else:
            ax.set_xticks([])

        for key, spine in ax.spines.items():
            spine.set_visible(False)
        ax.set_xlim(0, xmax)

    fig.legend(unq_celltypes, loc='center left', bbox_to_anchor=(1.005, 0.75), title="Cell Types")
    # plt.subplots_adjust(left=0.05, right=0.8, wspace=0.15)
    # fig.suptitle(f'Fraction of {common_cell_subtype} T cells by Condition')
    # plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout to fit suptitle
    if save_dir:
        plt.savefig(save_dir, bbox_inches='tight', dpi=300)
    plt.show()


#for counting tils and dlns 
def extract_count(adata):
    to_extract = ['condition', 'day', 'Location']
    df_tmp = adata.obs[to_extract].copy()
    df_counts = df_tmp.groupby(['condition', 'day'])['Location'].value_counts().unstack(fill_value=0)   
    df_counts_reset = df_counts.reset_index()
    return df_counts_reset

#plotting of til and dln count per condition
def plot_count(df_counts_reset, xmax = 60, counts = 'Absolute Counts', save_dir=None):
    conditions = np.unique(df_counts_reset.condition)

    fig, axs = plt.subplots(len(conditions), figsize=(8, 6))

    for idx, condition in enumerate(conditions):
        ax = axs[idx] if len(conditions) > 1 else axs  # Handles case with only one condition
        cond_subset = df_counts_reset[df_counts_reset.condition == condition]
        cond_subset = cond_subset.sort_values(by='day', ascending=False)

        colors = ['#E69F00', '#56B4E9', '#009E73']
        labels = ['TIL', 'dLN', 'in-vitro']
        categories = ['TIL', 'dLN', 'in-vitro']

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
    # file_path = os.path.join(save_dir, f'qc_metrics_p{i+1}.png')
    if save_dir:
        plt.savefig(save_dir, bbox_inches='tight', dpi=300)
    plt.show()


###aggregate to pseudobulk experiment
#create a pseudobulk RNA-seq object from the scRNA-seq object, ergo group by celltype, day and replicate and sum up gene counts of all cells within the given group
def aggregate_and_sum(diff_adata, obs_column, target_cell):
    pbs = []
    #subset the adata object for a celltype
    adata_subcell = diff_adata[diff_adata.obs[obs_column] == target_cell]
    #iterate trough each condition
    for condition in adata_subcell.obs.condition.unique():
        #subset the adata object containig only one celltype for the specific condition
        adata_subcell_cond = adata_subcell[adata_subcell.obs.condition == condition] 
        #iterate trough each day
        for day in adata_subcell_cond.obs.day.unique():
            #subset acording to each day
            adata_subcell_cond_day = adata_subcell_cond[adata_subcell_cond.obs.day == day]
            #iterate trough each replicate
            for rep in adata_subcell_cond_day.obs.replicate.unique():
                #subset according to each replicate
                adata_subcell_cond_day_rep = adata_subcell_cond_day[adata_subcell_cond_day.obs.replicate == rep]
                #for pool1/day0 we only have one replicate, split that up into randomely generated pseudo replicates of 3 (minimum, we do not have many cells so dont split too much)
                indices = list(adata_subcell_cond_day_rep.obs_names)
                random.shuffle(indices)
                if day == "0":
                    indices = np.array_split(np.array(indices), 3)
                else:
                    indices = np.array_split(np.array(indices), 1)
                for i, pseudo_rep_indx in enumerate(indices):
                    adata_subcell_cond_day_rep_pseudo = adata_subcell_cond_day_rep[pseudo_rep_indx]
                    #sum the gene counts for each gene together for all cells of all replicates or pseudo replicates
                    rep_adata = sc.AnnData(X = adata_subcell_cond_day_rep_pseudo.X.sum(axis=0), var = adata_subcell_cond_day_rep_pseudo.var[[]])
                    rep_adata.obs_names = [str(len(pbs))]
                    rep_adata.obs['condition'] = [condition]
                    rep_adata.obs['day'] = [day]
                    rep_adata.obs['replicate'] = [f'{rep}_{i+1}']
                    rep_adata.obs['nr_cells'] = [adata_subcell_cond_day_rep_pseudo.shape[0]]
                    pbs.append(rep_adata)
    pb = sc.concat(pbs)
    return pb


##do pydeseq2 for pseudobulk statistics
def compute_statistics(pb, contrast_vecs, equation = '~condition * day'):
    #create result dictionary, storing all result statistic tables
    res_dict = {key: None for key in contrast_vecs}
    #create counts table from pseudobulk anndata object
    counts = pd.DataFrame(pb.X, columns = pb.var_names)
    #create DeseqDataSet object from counts and pseudobulk object
    dds = DeseqDataSet(counts=counts, metadata=pb.obs, design=equation) #includes interaction terms
    #filter out low quality genes which are not present in at least one cell
    sc.pp.filter_genes(dds, min_cells=1)
    #run deseq2
    dds.deseq2() 
    #iterate through all given contrast vectors that shall be analysed
    for key, contrast_vec in contrast_vecs.items():
        #compute statistics for each given contrast
        stat_res = DeseqStats(dds, contrast=contrast_vec)
        #keep annoying messages from popping up
        original_stdout = sys.stdout
        sys.stdout = io.StringIO()
        #compute summary statistic
        stat_res.summary()
        sys.stdout = original_stdout
        #get statistic table
        de = stat_res.results_df
        de = de.reset_index()
        #ad a new metrix weighting pvalues and log2fc
        de['nlog10'] = -np.log10(de['padj']) #make nlog10 column
        de['sorter'] = de['nlog10']*de['log2FoldChange'] #make a column to pick top genes
        #add table to result dictionary
        res_dict[key] = de
    return res_dict



###############
#taken from sanbomics
###############
def volcano(data, log2fc='log2FoldChange', pvalue='padj', symbol='symbol',
           baseMean=None, pval_thresh=0.05, log2fc_thresh=0.75, to_label=5,
           color_dict=None, shape_dict=None, fontsize=10, colors={'not DE': 'lightgrey', 'DE': 'dimgrey', 'picked': 'black'}, #           ['dimgrey', 'lightgrey', 'black']
           top_right_frame=False, figsize=(5,5), legend_pos=(1.4,1),
           point_sizes=(20,200), save=False, shapes=None, shape_order=None,
           ax=None):  # NEW PARAMETER

    # Load data
    if isinstance(data, str):
        df = pd.read_csv(data)
    else:
        df = data.copy(deep=True)

    # Clean and replace zero p-values
    df = df.dropna()
    if df[pvalue].min() == 0:
        df[pvalue] = df[pvalue].replace(0, 1e-323)

    # Convert p-value threshold
    pval_thresh = -np.log10(pval_thresh)
    df['nlog10'] = -np.log10(df[pvalue])
    df['sorter'] = df['nlog10'] * df[log2fc]

    # Size the dots if baseMean is provided
    if baseMean is not None:
        df['logBaseMean'] = np.log(df[baseMean])
        baseMean = 'logBaseMean'
    else:
        point_sizes = None

    # Determine labeling genes
    if isinstance(to_label, int):
        label_df = pd.concat((df.sort_values('sorter')[-to_label:], df.sort_values('sorter')[:to_label]))
    else:
        label_df = df[df[symbol].isin(to_label)].copy()
        sbst = [abs(log) >= log2fc_thresh and pval >= pval_thresh for log,pval in zip(label_df[log2fc], label_df['nlog10'])]
        label_df = label_df[sbst]


    # Define color mapping
    def map_color(a):
        log2FoldChange, zymbol, nlog10 = a
        if zymbol in label_df[symbol].tolist():
            return 'picked'
        elif abs(log2FoldChange) < log2fc_thresh or nlog10 < pval_thresh:
            return 'not DE'
        # else:
        #     return 'DE'
        elif abs(log2FoldChange) >= log2fc_thresh and nlog10 >= pval_thresh:
            return 'DE'


    df['color'] = df[[log2fc, symbol, 'nlog10']].apply(map_color, axis=1)
    # hues = ['DE', 'not DE', 'picked'][:len(df.color.unique())]
    # std_hue = ['DE', 'not DE', 'picked']
    # hues = df.color.unique()
    

    # Set up figure and axis (MODIFICATION)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)  # Create new figure only if ax not provided

    # Plot scatter plot
    sns.scatterplot(data=df, x=log2fc, y='nlog10', hue='color', palette=colors, #, hue_order=hues
                    size=baseMean, sizes=point_sizes, ax=ax)

    # Add gene labels
    if label_df.shape[0] > 0:
        texts = []
        for i in range(len(label_df)):
            # if label_df.loc[label_df.index[i], 'color'] == 'not DE':
            # # if label_df['color'].iloc[i] == 'not DE':
            #     continue
            txt = ax.text(x=label_df.iloc[i][log2fc], y=label_df.iloc[i].nlog10, s=label_df.iloc[i][symbol],
                          fontsize=fontsize, weight='bold')
            txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='w')])
            texts.append(txt)
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='k', zorder=5), ax=ax)

    # Plot threshold lines
    ax.axhline(pval_thresh, zorder=0, c='k', lw=2, ls='--')
    ax.axvline(log2fc_thresh, zorder=0, c='k', lw=2, ls='--')
    ax.axvline(-log2fc_thresh, zorder=0, c='k', lw=2, ls='--')

    # Beautify plot
    for axis in ['bottom', 'left', 'top', 'right']:
        ax.spines[axis].set_linewidth(2)
    if not top_right_frame:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    ax.tick_params(width=2)
    ax.set_xlabel("$log_{2}$ fold change", size=15)
    ax.set_ylabel("-$log_{10}$(pval)", size=15)

    # Legend
    ax.legend(loc=1, bbox_to_anchor=legend_pos, frameon=False, prop={'weight': 'bold'})

    # Save if requested
    if save:
        plt.savefig(save + '.png', dpi=300, bbox_inches='tight')
        plt.savefig(save + '.svg', bbox_inches='tight')

    # Show only if we created the figure (i.e., not using subplot mode)
    if ax is None:
        plt.show()
    return df

