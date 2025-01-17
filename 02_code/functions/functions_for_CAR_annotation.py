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