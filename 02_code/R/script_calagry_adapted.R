
library('Seurat') #5.1.0
library('ProjecTILs') #3.3.1
# library('scRepertoire') #1.12.0
library('STACAS') #2.2.2
library('scGate') #1.6.2
library('dittoSeq') #1.14.3
library('Biostrings') #2.70.3
library('gprofiler2')
library('SeuratDisk')
library('SeuratData')
library('stringr')
set.seed(1234)


# install.packages(c("scGate", "BiocManager"))
# BiocManager::install(c("Biostrings","dittoSeq", "scRepertoire"))
# BiocManager::install("dittoSeq")
# BiocManager::install("scRepertoire")


######################################################
#functions for processing
######################################################


##################
#createSeuratObject
##################
import_data_to_seurat <- function(path){
    data <- list()
    data[[1]] <- Read10X(path)[['Gene Expression']]
    data[[2]] <- Read10X(path, gene.column = 1)[['Antibody Capture']]
    seurat_obj1 <- CreateSeuratObject(counts=data[[1]])
    seurat_obj1[['HTO']] <- CreateAssayObject(counts=data[[2]])
    return(seurat_obj1)
}


##################
#KeepNo-zeroHTOCounts
##################
eliminate_zero_counts <- function(seurat_obj) {
    counts <- GetAssayData(seurat_obj, assay = "HTO", layer = "counts")
    zero_count_cells <- colnames(counts[, colSums(counts) == 0])
    seurat_obj <- subset(seurat_obj, cells = setdiff(colnames(seurat_obj), zero_count_cells))
    return(seurat_obj)
}


##################
#Perform HTODemux
##################
demultiplex <- function(seurat_obj, pool_nr) {
    norm <- c(3,4,5,6,8)
    if (pool_nr %in% norm) {
        seurat_obj <- NormalizeData(seurat_obj, assay= "HTO", normalization.method= "CLR")
        print("normalizing")
    }
    seurat_obj <- eliminate_zero_counts(seurat_obj)
    if (pool_nr == 9) {
        seurat_obj <- HTODemux(seurat_obj, assay = "HTO", kfunc="clara", positive.quantile = 0.99)
        print("clara")
    }
    else {
        seurat_obj <- HTODemux(seurat_obj, assay = "HTO", kfunc="kmeans", positive.quantile = 0.99)
        print("kmeans")
    }
    return(seurat_obj)
}


##################
#Delete Singlets
##################
delete_singlets <- function(seurat_obj) {
    singlets <- subset(seurat_obj, subset=HTO_classification.global == "Singlet")
    hto_classification <- singlets@meta.data$HTO_classification
    return(singlets)
}


##################
#Perform Quality Control
##################
quality_control <- function(singlets) {
    singlets[["percent.mt"]] <-PercentageFeatureSet(singlets, pattern = "^mt-")
    singlets$percent.mt[is.nan(singlets$percent.mt)] <- 0
    singlets <-subset(singlets, subset=nFeature_RNA >200 & nFeature_RNA <6000 & percent.mt <20)
    return(singlets)
}


##################
#Perform cell cycle correction and normalize
##################
further_processing <- function(singlets) {
    singlets <- FindVariableFeatures(singlets)
    mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
    mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
    singlets <- NormalizeData(singlets, assay = "RNA")
    common_s <- intersect(mmus_s, rownames(singlets))
    common_g2m <- intersect(mmus_g2m, rownames(singlets))
    singlets <- CellCycleScoring(singlets, s.features = common_s, g2m.features = common_g2m)
    singlets <- ScaleData(singlets, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
    singlets <- RunPCA(singlets, features = VariableFeatures(object = singlets))
    singlets <- FindNeighbors(singlets, dims = 1:15)
    singlets <- FindClusters(singlets, resolution = 0.5)
    singlets <- RunUMAP(singlets, dims = 1:15)
    return(singlets)
}


##################
# ###ProjecTILs on each sample respectively
##################

project_TILs_single_pool <- function(singlets) {
    ref.cd8 <- load.reference.map("~/car_t_sc/01_data/reference_datasets_project_TILs/CD8T_human_ref_v1.rds")
    ref.cd4 <- load.reference.map("~/car_t_sc/01_data/reference_datasets_project_TILs/CD4T_human_ref_v2.rds")
    ncores = 8
    DefaultAssay(ref.cd4) <- "integrated"
    query.projected <- ProjecTILs.classifier(query= singlets, ref= ref.cd4, reduction= "umap", ncores = ncores, split.by="hash.ID")
    # table(query.projected$functional.cluster, useNA="ifany")
    DefaultAssay(ref.cd8) <- "integrated"
    query.projected <- ProjecTILs.classifier(query= query.projected, ref= ref.cd8, reduction= "umap", ncores = ncores, split.by="hash.ID", overwrite=FALSE)
    # table(query.projected$functional.cluster, useNA="ifany")
    genes4radar <- c("Mki67", "Foxp3", "Cd4", "Cd8a", "Tcf7", "Ccr7", "Gzmb", "Gzmk", "Pdcd1", "Havcr2", "Tox")
    # Check and remove cells with NA in functional.cluster
    query.projected <- query.projected[, !is.na(query.projected@meta.data$functional.cluster)]
    # print(sum(is.na(query.projected@meta.data$functional.cluster)))
    # Extract genes_interested seurat object for RadarPlot
    new_obj <- subset(query.projected, features=genes4radar)
    gene_names <- rownames(new_obj)
    uppercase_gene_names <- toupper(gene_names)
    rownames(new_obj@assays$RNA) <- uppercase_gene_names
    rownames(new_obj)
    # plot.states.radar(ref.cd4, new_obj, genes4radar=uppercase_gene_names, min.cells=1)
    # plot.states.radar(ref.cd8, new_obj, genes4radar=uppercase_gene_names, min.cells=1)
}



##################
#Apply functions for all objects/pools
##################
path_to_folders <- '~/car_t_sc/01_data/raw/cellranger_multi'
path_to_count_matrices_in_folders <- 'per_sample_outs/count/sample_filtered_feature_bc_matrix'
output <- "~/car_t_sc/01_data/processed/preprocessed_pools_R/RData_files/processed_pseudocount"
# path_to_raw_count_matrices_in_folders <- 'count/raw_feature_bc_matrix'

# pools <- list()
subfolders <- list.dirs(path_to_folders, recursive = FALSE)
for (i in 1:length(subfolders)) {
    folder <- subfolders[i]
    count_matrix_path <- file.path(folder, path_to_count_matrices_in_folders)
    seurat_obj <- import_data_to_seurat(count_matrix_path)
    seurat_obj <- demultiplex(seurat_obj, i)
    singlets <- subset(seurat_obj, subset=HTO_classification.global == "Singlet")
    singlets <- quality_control(singlets)
    singlets <- further_processing(singlets)
    filename <- paste0("P", i, "_singlets_processed.RData")
    file_path <- file.path(output, filename)
    save(singlets, file = file_path)
    # SaveH5Seurat(singlets, filename = file_path)
}


##################
#for a single object, e.g. pool 1 (P1)
##################
#load data
path_p1 <- "~/car_t_sc/01_data/raw/cellranger_multi/2024-06-07_24054SC_Luu_P1_cellranger/per_sample_outs/count/sample_filtered_feature_bc_matrix"
seurat_obj <- import_data_to_seurat(path_p1)
#demultiplex
seurat_obj <- demultiplex(seurat_obj, 1)
print(table(seurat_obj$HTO_classification.global))
print(table(seurat_obj$HTO_classification))
# FeatureScatter(seurat_obj, feature1= "d0-C", feature2="d0-P")
#delete singlets
singlets <- subset(seurat_obj, subset=HTO_classification.global == "Singlet")
# table(singlets@meta.data$HTO_classification)
#do quality control
singlets <- quality_control(singlets)
# VlnPlot(singlets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#normalize and correct the cell cycle
singlets <- further_processing(singlets)
# DimPlot(singlets, reduction = "umap",group.by= "hash.ID", label = TRUE)
#safe the processed data
output <- "~/car_t_sc/01_data/processed/preprocessed_pools_R"
filename <- "P1_singlets_processed.RData"
file_path <- file.path(output, filename)
save(singlets, file = file_path)
# SaveH5Seurat(singlets, filename = file_path, overwrite = TRUE)


############################################################################################################################################
#end of processing single pools
############################################################################################################################################

############################################################################################################################################
#beginning with merging all pools into one object
############################################################################################################################################

##################
#continue with importing safed pools
##################
#import all the data
pools = list()
folder <- "~/car_t_sc/01_data/processed/preprocessed_pools_R/RData_files"
files <- list.files(folder)
# file <- file.path(folder, "P1_singlets_processed.RData")
for (i in 1:length(files)) {
    path_to_file <- file.path(folder, files[i])
    load(path_to_file)
    pools[[i]] <- singlets
}

#this is for h5 files, which apparently does not properly work in R
# load(file)
# p1 <- singlets
# p1 <- LoadH5Seurat(file, assays = list(RNA = "assay.orig"))
# h5_file <- H5File$new(file, mode = "r")
# rna_group <- h5_file[["assays/RNA"]]
# rna_group$ls()


##################
# ###IntegrationOfNineObjects
##################
AAA <- lapply(seq_along(pools), function(i) {
  x <- pools[[i]]
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 6000)
  pool_name <- paste0("P", i)  # Create a name for the pool, e.g., "Pool1", "Pool2", etc.
  x <- AddMetaData(x, metadata = pool_name, col.name = "pool")
  return(x)
})


integrated_data <- merge(AAA[[1]], AAA[2:9])
integrated_data <- JoinLayers(integrated_data, assay = "RNA")

anchors <- FindIntegrationAnchors(object.list = AAA, dims=1:30)
integrated_data <- IntegrateData(anchorset = anchors, dims=1:30)


###safe the merged object
file_path <- "~/car_t_sc/01_data/processed/merged_and_processed/just_merged_calagry_exact.RData"
save(singlets, file = file_path)

load(file_path)

##################
# #pre-processing of merged data
##################
integrated_data <- NormalizeData(integrated_data)
integrated_data <- ScaleData(integrated_data, verbose = FALSE)
integrated_data <- RunPCA(integrated_data, npcs = 30, verbose = FALSE)
integrated_data <- FindNeighbors(integrated_data, reduction = "pca", dims = 1:30)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)
integrated_data <- RunUMAP(integrated_data, dims = 1:30)

##################
#add meta data to merged object from object list, since it is not kept when integrating
##################
metadata_combined <- do.call(rbind, lapply(AAA, function(x) x@meta.data))
# rownames(metadata_combined) <- metadata_combined$cell
integrated_data <- AddMetaData(integrated_data, metadata_combined)
# DimPlot(integrated_data, reduction="umap", group.by="pool")

##################
#extract tne days, conditions and tissue information. store in a new meta.data object
##################
integrated_data@meta.data$condition <- str_extract(integrated_data@meta.data$hash.ID, "C|P|DM")
integrated_data@meta.data$day <- str_extract(integrated_data@meta.data$hash.ID, "0|7|14")
integrated_data@meta.data$Location <- str_extract(integrated_data@meta.data$hash.ID, "TIL|dLN")
integrated_data@meta.data$Location[is.na(integrated_data@meta.data$Location)] <- "in-vitro"

##################
# #RemovingAll"dLN"
##################
Integrated <- subset(integrated_data, subset = Location != "dLN")
# Integrated <- subset(integrated_data, cells = rownames(integrated_data@meta.data)[!grepl("dLN", integrated_data@meta.data$HTO_classification)])


##################
# #now using the merged RNA assay for further STACAS integration, not the intgrated assay, so in theory, the whole process is just in order to merge the seurat objects 
##################
# DefaultAssay(Integrated) <- "RNA"
Integrated <- JoinLayers(Integrated, assay = "RNA")
Integrated <- NormalizeData(Integrated)
Integrated <- FindVariableFeatures(Integrated)
Integrated <- ScaleData(Integrated, verbose=FALSE)
Integrated <- RunPCA(Integrated, npcs=30, verbose=FALSE)
Integrated <- FindNeighbors(Integrated, reduction = "pca", dims = 1:20)
Integrated <- FindClusters(Integrated, resolution = 0.5)
Integrated <- RunUMAP(Integrated, dims = 1:20)
DimPlot(Integrated, reduction="umap", group.by="Location")


##################
# ###STACASIntegration
##################
stacas <- NormalizeData(Integrated) |>
 SplitObject(split.by="pool") |>
 Run.STACAS()
stacas <-RunUMAP(stacas, dims=1:30)
DimPlot(stacas, reduction="umap", group.by="pool")


##################
# ###ScgateForCellSelection
##################
scgate_DB <- Integrated
scgate_DB <- get_scGateDB()
stacas_scgate <- scGate(stacas, model=scgate_DB$mouse$generic)
DimPlot(stacas_scgate, group.by="is.pure_Tcell")
DimPlot(stacas_scgate, group.by="scGate_multi")
FeaturePlot(stacas_scgate, features=c("Tcell_UCell", "CD4T_UCell", "CD8T_UCell"))

##################
# #KeepPureTcell
##################
stacas_Tcells <- subset(stacas_scgate, subset = is.pure_Tcell == "Pure")
DefaultAssay(stacas_Tcells) <- "RNA"
stacas_Tcells <- NormalizeData(stacas_Tcells)

#new
stacas_Tcells <- RunPCA(stacas_Tcells)
stacas_Tcells <- FindNeighbors(stacas_Tcells, dims = 1:15)
stacas_Tcells <- FindClusters(stacas_Tcells, resolution = 0.5)
stacas_Tcells <- RunUMAP(stacas_Tcells, dims = 1:15)
DimPlot(stacas_Tcells, reduction = "umap",group.by= "pool")

##################
# #RegressOutCellCycle&MitochondrialGenes
##################
stacas_Tcells <- ScaleData(stacas_Tcells, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
stacas_Tcells <- FindVariableFeatures(stacas_Tcells)

##################
# ###STACASIntegration again
##################
stacas_Tcells2 <- NormalizeData(stacas_Tcells) |>
 SplitObject(split.by="pool") |>
 Run.STACAS()
stacas_Tcells2 <- RunPCA(stacas_Tcells2)
# ElbowPlot(stacas_Tcells, ndims=50)
stacas_Tcells2 <- FindNeighbors(stacas_Tcells2, dims = 1:15)
stacas_Tcells2 <- FindClusters(stacas_Tcells2, resolution = 0.5)
stacas_Tcells2 <- RunUMAP(stacas_Tcells2, dims = 1:15)
DimPlot(stacas_Tcells2, reduction = "umap",group.by= "pool")


file_path <- "~/car_t_sc/01_data/processed/merged_and_processed/merged_integrated_xTcell_subtypes_calagry_exact.RData"
save(singlets, file = file_path)


load(file_path)

##################
# ###ProjecTILsDefaultParameters
##################
# #Mapping2HumanRef
DefaultAssay(stacas_Tcells) <- "RNA"
ref.cd8 <- load.reference.map("~/car_t_sc/01_data/reference_datasets_project_TILs/CD8T_human_ref_v1.rds")
ref.cd4 <- load.reference.map("~/car_t_sc/01_data/reference_datasets_project_TILs/CD4T_human_ref_v2.rds")
ncores = 8
DefaultAssay(ref.cd4) <- "integrated"
scgate.projected <- ProjecTILs.classifier(query= stacas_Tcells, ref= ref.cd4, reduction= "umap") 
DefaultAssay(ref.cd8) <- "integrated"
scgate.projected <- ProjecTILs.classifier(query= scgate.projected, ref= ref.cd8, reduction= "umap", overwrite=FALSE)
scgate.projected$Day.Condition <- paste(scgate.projected$day, scgate.projected$condition, sep="_")
table(scgate.projected$pool, scgate.projected$Day.Condition)
table(scgate.projected$condition, scgate.projected$day)
table(scgate.projected$functional.cluster, useNA = "ifany")
DimPlot(scgate.projected, group.by="functional.cluster", cols=c("CD4.CTL_Exh"= "#E69F00", "CD4.CTL_GNLY"="#56B4E9", "CD4.NaiveLike"= "#009E73", "CD4.Tfh"= "#F0E442", "CD4.Th17"= "#0072B2", "CD4.Treg"="#D55E00", "CD4.CTL_EOMES"="#CC79A7", "CD8.CM"="#800000", "CD8.EM"= "#9ACD32", "CD8.MAIT"="#2F4F4F", "CD8.NaiveLike"="#FF4500", "CD8.TEMRA"= "#FF6347", "CD8.TEX"= "#9400D3", "CD8.TPEX"= "#7FFFD4"))
# #RemoveNAfunction
new_obj <- scgate.projected[, !is.na(scgate.projected@meta.data$functional.cluster)]
new_obj <- ScaleData(new_obj)
# #DittoBarPlot
# dittoBarPlot(new_obj, x.reorder=c(1,3,2), group.by="day", var="NewAnnotation", split.by="condition")
# dittoBarPlot(new_obj, x.reorder=c(1,3,2), group.by="day", var="functional.cluster", split.by="condition")
cd4_cells <- rownames(new_obj@meta.data)[grepl("CD4", new_obj@meta.data$functional.cluster)]
new_cd4 <- subset(new_obj, cells = cd4_cells)
cd8_cells <- rownames(new_obj@meta.data)[grepl("CD8", new_obj@meta.data$functional.cluster)]
new_cd8 <- subset(new_obj, cells = cd8_cells)
VlnPlot(new_obj, features=c("Tnf", "Ifng", "Il2", "Gzmb", "Tbx21", "Foxp3"), group.by="day", split.by="condition", pt.size=0)
# VlnPlot(new_cd4, features=c("Tnf", "Ifng", "Il2", "Gzmb", "Tbx21", "Foxp3"), group.by="day", split.by="condition", pt.size=0)
# VlnPlot(new_cd8, features=c("Tnf", "Ifng", "Il2", "Gzmb", "Tbx21", "Foxp3"), group.by="day", split.by="condition", pt.size=0)
dittoBarPlot(new_cd4, x.reorder=c(1,3,2), group.by="day", var="functional.cluster", split.by="condition") | dittoBarPlot(new_cd8, x.reorder=c(1,3,2), group.by="day", var="functional.cluster", split.by="condition")

# save.image("/mnt/raw-seq/Maik-scRNAMouse/Analysis/Integrated.RData")




DefaultAssay(stacas_Tcells2) <- "RNA"
scgate.projected2 <- ProjecTILs.classifier(query= stacas_Tcells2, ref= ref.cd4, reduction= "umap", skip.normalize=TRUE) 
scgate.projected2 <- ProjecTILs.classifier(query= scgate.projected2, ref= ref.cd8, reduction= "umap", overwrite=FALSE)
scgate.projected2$Day.Condition <- paste(scgate.projected2$day, scgate.projected2$condition, sep="_")
table(scgate.projected2$functional.cluster, useNA = "ifany")
DimPlot(scgate.projected2, group.by="functional.cluster", cols=c("CD4.CTL_Exh"= "#E69F00", "CD4.CTL_GNLY"="#56B4E9", "CD4.NaiveLike"= "#009E73", "CD4.Tfh"= "#F0E442", "CD4.Th17"= "#0072B2", "CD4.Treg"="#D55E00", "CD4.CTL_EOMES"="#CC79A7", "CD8.CM"="#800000", "CD8.EM"= "#9ACD32", "CD8.MAIT"="#2F4F4F", "CD8.NaiveLike"="#FF4500", "CD8.TEMRA"= "#FF6347", "CD8.TEX"= "#9400D3", "CD8.TPEX"= "#7FFFD4"))
