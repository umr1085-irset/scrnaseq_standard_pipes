#########################
# Pre-processing pipeline
# 2021 IRSET
#########################

args = commandArgs(trailingOnly=TRUE)

input_seurat_obj_path = args[1]
save_seurat_obj_path = args[2]

library(Seurat)
library(SingleCellExperiment)
library(scater)

# read pre_pipe output file
obj <- readRDS(file=input_seurat_obj_path)

# Remove outliers and doublets
obj <- subset(x = obj, subset = scateroutlier != 'TRUE') # keep non-outlier cells
obj <- subset(x = obj, subset = DoubletFinder != 'Doublet') # keep singlets

# cleanup step
obj[['umap']] <- NULL
obj[['pca']] <- NULL
DefaultAssay(object = obj) <- "RNA"
obj[['SCT']] <- NULL
obj@meta.data[['nCount_SCT']] <- NULL
obj@meta.data[['nFeature_SCT']] <- NULL
obj@meta.data[['S.Score']] <- NULL
obj@meta.data[['G2M.Score']] <- NULL
obj@meta.data[['old.ident']] <- NULL
obj@meta.data[['seurat_clusters']] <- NULL
obj@meta.data[['Phase']] <- NULL
lbl <- names(obj@meta.data)[startsWith(names(obj@meta.data),'SCT_snn_res')]
obj@meta.data[[lbl]] <- NULL

# Keep genes with more than 10 expressing cells
numgenes <- nexprs(GetAssayData(object = obj, slot = "counts"), byrow=TRUE)
obj <- obj[numgenes >= 10,]

# SCTransform
obj = SCTransform(obj, vars.to.regress = "percent.mito", verbose = FALSE) # run Seurat sctransform method
s.genes <- cc.genes$s.genes # extract genes associated to S cycle
g2m.genes <- cc.genes$g2m.genes # extract genes associated to G2M cycle
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident=TRUE) # compute cell cyle scores for all cells

# Dim Red
obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs=100)
nDims = 30
obj <- RunUMAP(obj, dims = 1:nDims) # run UMAP on PCA
obj <- FindNeighbors(obj, dims = 1:nDims) # build knn graph then snn graph
obj <- FindClusters(obj) # Louvain clustering

# Save Seurat object
saveRDS(obj, file = save_seurat_obj_path)
