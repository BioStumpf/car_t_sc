
library('Seurat') #5.1.0
library('ProjecTILs') #3.3.1
library('scRepertoire') #1.12.0
library('STACAS') #2.2.2
library('scGate') #1.6.2
library('dittoSeq') #1.14.3
library('Biostrings') #2.70.3
library('gprofiler2')
set.seed(1234)

linters: with_defaults(line_length_linter = line_length_linter(5000))

##################
###LoadFileFrom10X with correct HTO names (sample instead of HTO label)
##################
# "~/car_t_sc/01_data/raw/cellranger_multi/2024-06-07_24054SC_Luu_P1_cellranger/per_sample_outs/count/sample_filtered_feature_bc_matrix"
path_to_folders <- '~/car_t_sc/01_data/raw/cellranger_multi'
path_to_count_matrices_in_folders <- 'per_sample_outs/count/sample_filtered_feature_bc_matrix'
# path_to_raw_count_matrices_in_folders <- 'count/raw_feature_bc_matrix'

pools <- list()

subfolders <- list.dirs(path_to_folders)
for (i in 1:length(subfolders)) {
    folder <- subfolders[i]
    count_matrix_path <- path.join(folder, path_to_count_matrices_in_folders)
    seurat_obj <- import_data_to_seurat(count_matrix_path)
    seurat_obj <- demultiplex(seurat_obj, i)
    singlets <- delete_singlets(seurat_obj)
    singlets <- quality_control(singlets)
}

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
    norm =c(3,4,5,6,8)
    if (pool_nr in norm) {
        seurat_obj <- NormalizeData(seurat_obj, assay= "HTO", normalization.method= "CLR")
    }
    seurat_obj <- eliminate_zero_counts(seurat_obj)
    if (i == 9) {
        seurat_obj <- HTODemux(seurat_obj, assay = "HTO", kfunc="clara", positive.quantile = 0.99)
    }
    else {
        seurat_obj <- HTODemux(seurat_obj, assay = "HTO", kfunc="kmeans", positive.quantile = 0.99)
    }
    return(seurat_obj)
}

# seurat_obj <- HTODemux(seurat_obj, assay = "HTO", kfunc="kmeans", positive.quantile = 0.99)
# print(table(seurat_obj$HTO_classification.global))
# print(table(seurat_obj$HTO_classification))
# FeatureScatter(seurat_obj, feature1= "Ms.Hashtag-1", feature2="Ms.Hashtag-2")


##################
#Delete Singlets
##################
delete_singlets <- function(seurat_obj) {
    singlets <- subset(seurat_obj, subset=HTO_classification.global == "Singlet")
    hto_classification <- singlets@meta.data$HTO_classification
    return(singlets)
}
# table(hto_classification)


##################
#Perform Quality Control
##################
quality_control <- function(singlets) {
    singlets[["percent.mt"]] <-PercentageFeatureSet(singlets, pattern = "^mt-")
    singlets$percent.mt[is.nan(singlets$percent.mt)] <- 0
    singlets <-subset(singlets, subset=nFeature_RNA >200 & nFeature_RNA <6000 & percent.mt <20)
    return(singlets)
}
# VlnPlot(singlets_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)


##################
#Perform cell cycle correction and normalize
##################
singlets_P1 <- FindVariableFeatures(singlets_P1)
#ConvertHumanCellCycleGenestoMouse
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
singlets_P1 <- NormalizeData(singlets_P1, assay = "RNA")
common_s <- intersect(mmus_s, rownames(singlets_P1))
common_g2m <- intersect(mmus_g2m, rownames(singlets_P1))
singlets_P1 <- CellCycleScoring(singlets_P1, s.features = common_s, g2m.features = common_g2m)
singlets_P1 <- ScaleData(singlets_P1, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
singlets_P1 <- RunPCA(singlets_P1, features = VariableFeatures(object = singlets_P1))
# ElbowPlot(singlets_P1, ndims=50)
singlets_P1 <- FindNeighbors(singlets_P1, dims = 1:15)
singlets_P1 <- FindClusters(singlets_P1, resolution = 0.5)
singlets_P1 <- RunUMAP(singlets_P1, dims = 1:15)
# DimPlot(singlets_P1, reduction = "umap",group.by= "hash.ID", label = TRUE)


##################
# ###ProjecTILs on each sample respectively
##################
ref.cd8 <- load.reference.map("01_data/reference_datasets_project_TILs/CD8T_human_ref_v1.rds")
ref.cd4 <- load.reference.map("01_data/reference_datasets_project_TILs/CD4T_human_ref_v2.rds")
ncores = 8
DefaultAssay(ref.cd4) <- "integrated"
query.projected <- ProjecTILs.classifier(query= singlets_P1, ref= ref.cd4, reduction= "umap", ncores = ncores, split.by="hash.ID")
table(query.projected$functional.cluster, useNA="ifany")
DefaultAssay(ref.cd8) <- "integrated"
query.projected <- ProjecTILs.classifier(query= query.projected, ref= ref.cd8, reduction= "umap", ncores = ncores, split.by="hash.ID", overwrite=FALSE)
table(query.projected$functional.cluster, useNA="ifany")
genes4radar <- c("Mki67", "Foxp3", "Cd4", "Cd8a", "Tcf7", "Ccr7", "Gzmb", "Gzmk", "Pdcd1", "Havcr2", "Tox")
# Check and remove cells with NA in functional.cluster
query.projected <- query.projected[, !is.na(query.projected@meta.data$functional.cluster)]
print(sum(is.na(query.projected@meta.data$functional.cluster)))
# Extract genes_interested seurat object for RadarPlot
new_obj <- subset(query.projected, features=genes4radar)
gene_names <- rownames(new_obj)
uppercase_gene_names <- toupper(gene_names)
rownames(new_obj@assays$RNA) <- uppercase_gene_names
rownames(new_obj)
plot.states.radar(ref.cd4, new_obj, genes4radar=uppercase_gene_names, min.cells=1)
plot.states.radar(ref.cd8, new_obj, genes4radar=uppercase_gene_names, min.cells=1)


###################
#safe data
###################
# save.image("/mnt/raw-seq/Maik-scRNAMouse/2024-06-07_24054SC_Luu_P1_cellranger/P1_workspace.RData")


###################
#continue with integrating all pools
###################
#load("/path/to/saved.RData")

#AddMetaInfo
# ###IntegrationOfNineObjects
# AAA <- list(P1, P2, P3, P4, P5, P6, P7, P8, P9)
# AAA <- lapply(AAA, function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 6000)
#   return(x)
# })
# anchors <- FindIntegrationAnchors(object.list = AAA, dims=1:30)
# integrated_data <- IntegrateData(anchorset = anchors, dims=1:30)
# #pre-processing
# integrated_data <- ScaleData(integrated_data, verbose = FALSE)
# integrated_data <- RunPCA(integrated_data, npcs = 30, verbose = FALSE)
# integrated_data <- FindNeighbors(integrated_data, reduction = "pca", dims = 1:30)
# integrated_data <- FindClusters(integrated_data, resolution = 0.5)
# integrated_data <- RunUMAP(integrated_data, dims = 1:30)
# metadata_combined <- do.call(rbind, lapply(AAA, function(x) x@meta.data))
# rownames(metadata_combined) <- metadata_combined$cell
# integrated_data <- AddMetaData(integrated_data, metadata_combined)
# DimPlot(integrated_data, reduction="umap", group.by="pool")
# #RemovingAll"dLN"
# Integrated <- subset(integrated_data, subset = Location != "dLN")
# DefaultAssay(Integrated) <- "RNA"
# Integrated <- JoinLayers(Integrated, assay = "RNA")
# Integrated <- NormalizeData(Integrated)
# Integrated <- FindVariableFeatures(Integrated)
# Integrated <- ScaleData(Integrated, verbose=FALSE)
# Integrated <- RunPCA(Integrated, npcs=30, verbose=FALSE)
# Integrated <- FindNeighbors(Integrated, reduction = "pca", dims = 1:20)
# Integrated <- FindClusters(Integrated, resolution = 0.5)
# Integrated <- RunUMAP(Integrated, dims = 1:20)
# DimPlot(Integrated, reduction="umap", group.by="Location"

# ###STACASIntegration
# stacas <- NormalizeData(Integrated) |>
#  SplitObject(split.by="pool") |>
#  Run.STACAS()
# stacas <-RunUMAP(stacas, dims=1:30)
# DimPlot(stacas, reduction="umap", group.by="pool")

# ###ScgateForCellSelection
# scgate_DB <- get_scGateDB()
# stacas_scgate <- scGate(stacas, model=scgate_DB$mouse$generic)
# DimPlot(stacas_scgate, group.by="is.pure_Tcell")
# DimPlot(stacas_scgate, group.by="scGate_multi")
# FeaturePlot(stacas_scgate, features=c("Tcell_UCell", "CD4T_UCell", "CD8T_UCell"))
# #KeepPureTcell
# stacas_Tcells <- subset(stacas_scgate, subset = is.pure_Tcell == "Pure")
# DefaultAssay(stacas_Tcells) <- "RNA"
# stacas_Tcells <- NormalizeData(stacas_Tcells)
# #RegressOutCellCycle&MitochondrialGenes
# stacas_Tcells <- CellCycleScoring(stacas_Tcells, s.features = common_s, g2m.features = common_g2m)
# stacas_Tcells <- ScaleData(stacas_Tcells, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
# stacas_Tcells <- FindVariableFeatures(stacas_Tcells)
# stacas_Tcells <- NormalizeData(stacas_Tcells) |>
#  SplitObject(split.by="pool") |>
#  Run.STACAS()
# stacas_Tcells <- RunPCA(stacas_Tcells)
# ElbowPlot(stacas_Tcells, ndims=50)
# stacas_Tcells <- FindNeighbors(stacas_Tcells, dims = 1:15)
# stacas_Tcells <- FindClusters(stacas_Tcells, resolution = 0.5)
# stacas_Tcells <- RunUMAP(stacas_Tcells, dims = 1:15)
# DimPlot(stacas_Tcells, reduction = "umap",group.by= "pool")

# ###ProjecTILsDefaultParameters
# #Mapping2HumanRef
# DefaultAssay(ref.cd4) <- "integrated"
# scgate.projected <- ProjecTILs.classifier(query= stacas_Tcells, ref= ref.cd4, reduction= "umap") 
# DefaultAssay(ref.cd8) <- "integrated"
# scgate.projected <- ProjecTILs.classifier(query= scgate.projected, ref= ref.cd8, reduction= "umap", overwrite=FALSE)
# scgate.projected$Day.Condition <- paste(scgate.projected$Day, scgate.projected$Condition, sep="_")
# table(scgate.projected$pool, scgate.projected$Day.Condition)
# table(scgate.projected$Condition, scgate.projected$Day)
# DimPlot(scgate.projected, group.by="functional.cluster", cols=c("CD4.CTL_Exh"= "#E69F00", "CD4.CTL_GNLY"="#56B4E9", "CD4.NaiveLike"= "#009E73", "CD4.Tfh"= "#F0E442", "CD4.Th17"= "#0072B2", "CD4.Treg"="#D55E00", "CD4.CTL_EOMES"="#CC79A7", "CD8.CM"="#800000", "CD8.EM"= "#9ACD32", "CD8.MAIT"="#2F4F4F", "CD8.NaiveLike"="#FF4500", "CD8.TEMRA"= "#FF6347", "CD8.TEX"= "#9400D3", "CD8.TPEX"= "#7FFFD4"))
# #RemoveNAfunction
# new_obj <- scgate.projected[, !is.na(scgate.projected@meta.data$functional.clust
# er)]
# new_obj <- ScaleData(new_obj)
# #DittoBarPlot
# dittoBarPlot(new_obj, x.reorder=c(1,3,2), group.by="Day", var="NewAnnotation", split.by="Condition")
# dittoBarPlot(new_obj, group.by="Day", var="functional.cluster", split.by="Condition")
# dittoBarPlot(new_obj, x.reorder=c(1,3,2), group.by="Day", var="functional.cluster", split.by="Condition")
# new_cd4 <- subset(new_obj, subset= NewAnnotation =="CD4")
# new_cd8 <- subset(new_obj, subset= NewAnnotation =="CD8")
# VlnPlot(new_obj, features=c("Tnf", "Ifng", "Il2", "Gzmb", "Tbx21", "Foxp3"), group.by="Day", split.by="Condition", pt.size=0)
# VlnPlot(new_cd4, features=c("Tnf", "Ifng", "Il2", "Gzmb", "Tbx21", "Foxp3"), group.by="Day", split.by="Condition", pt.size=0)
# VlnPlot(new_cd8, features=c("Tnf", "Ifng", "Il2", "Gzmb", "Tbx21", "Foxp3"), group.by="Day", split.by="Condition", pt.size=0)
# barplot <- dittoBarPlot(new_cd4, x.reorder=c(1,3,2), group.by="Day", var="functional.cluster", split.by="Condition") | dittoBarPlot(new_cd8, x.reorder=c(1,3,2), group.by="Day", var="functional.cluster", split.by="Condition")

# save.image("/mnt/raw-seq/Maik-scRNAMouse/Analysis/Integrated.RData")
