library(BiocManager)
library("Seurat") #5.1.0
library('ProjecTILs') #3.3.1
# library('scRepertoire') #1.12.0
library('STACAS') #2.2.2
library('scGate') #1.6.2
library('dittoSeq') #1.14.3
# library('Biostrings') #2.70.3
library('gprofiler2')
# library('SeuratDisk')
# library('SeuratData')
library('stringr')
library('sceasy')
library('remotes')
library('scImmunuCC')
set.seed(1234)


if (!requireNamespace("GSVA", quietly = TRUE))
    BiocManager::install("GSVA")

remotes::install_github("wuaipinglab/scImmuCC")

file_path <- "~/car_t_sc/01_data/processed/merged_and_processed/XXXCAR_genome/XXXCAR_genome_after_qc_TIL_only_pure_TC_annotation_non_TC_filtered.RData"
load(file_path)

