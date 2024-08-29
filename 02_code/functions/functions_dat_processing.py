import numpy as np
from scipy.stats import median_abs_deviation

#this function is to determine outliers. It can only be used after using sc.pp.calculate_qc_metrics().
#it takes a specific qc metric, like total count or total gene count and extracts the column corresponding to this metric form the adata object.
#then it computes wether the MAD (median absolute deviation) deviates nmads (e.g. 5) from the median. The output is an array of F/T.
def is_mad_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def is_pctl_outlier(adata, metric: str, pctl_threshhold: float):
    M = adata.obs[metric]
    lower_lim = np.quantile(M, pctl_threshhold)
    upper_lim = np.quantile(M, (1 - pctl_threshhold))
    outlier = (M < lower_lim) | (M > upper_lim)
    return outlier