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
library(DoubletFinder)

obj <- readRDS(file=input_seurat_obj_path)

# Add ribo and mito metadata

mito.genes <- rownames(obj)[grep("^[Mm][Tt]-",rownames(obj))]
percent.mito <- colSums(GetAssayData(object = obj, slot = "counts")[mito.genes,])/Matrix::colSums(GetAssayData(object = obj, slot = "counts"))*100
obj <- AddMetaData(obj, percent.mito, col.name = "percent.mito")

ribo.genes <- rownames(obj)[grep("^R[Pp][SLsl]",rownames(obj))]
percent.ribo <- colSums(GetAssayData(object = obj, slot = "counts")[ribo.genes,])/Matrix::colSums(GetAssayData(object = obj, slot = "counts"))*100
obj <- AddMetaData(obj, percent.ribo, col.name = "percent.ribo")

# Keep cells with more than 200 features

obj <- obj[,obj@meta.data$nFeature_RNA>=200]
obj

# Keep genes with more than 10 expressing cells

numgenes <- nexprs(GetAssayData(object = obj, slot = "counts"), byrow=TRUE)
obj <- obj[numgenes >= 10,] 
obj

# Scater

obj.sce <- as.SingleCellExperiment(obj)
obj.sce  <- runColDataPCA(obj.sce,variables=c("percent.mito","percent.ribo","nCount_RNA","nFeature_RNA"),outliers=TRUE,name='PCA_coldata') # run pca on col data
obj@meta.data$scateroutlier <- colData(obj.sce)$outlier
obj@meta.data$scateroutlierPC1 <- as.vector(reducedDim(obj.sce, "PCA_coldata")[,1])
obj@meta.data$scateroutlierPC2 <- as.vector(reducedDim(obj.sce, "PCA_coldata")[,2])

# DoubletFinder

dfcols <- data.frame(matrix(ncol = 1, nrow = 0))
colnames(dfcols) <- "DoubletFinder"

for (sample in unique(obj@meta.data$orig.ident)){
    print(sample)
    #obj.sub <- subset(x = obj, subset = orig.ident == sample)
    obj.sub = obj[,obj@meta.data$orig.ident==sample]
    obj.sub <- NormalizeData(obj.sub)
    obj.sub <- ScaleData(obj.sub)
    obj.sub <- FindVariableFeatures(obj.sub, selection.method = "vst", nfeatures = 2000)
    obj.sub <- RunPCA(obj.sub)
    sweep.res.list_sub.seurat <- paramSweep_v3(obj.sub, PCs = 1:10, sct = FALSE)
    sweep.stats_sub.seurat<- summarizeSweep(sweep.res.list_sub.seurat, GT = FALSE)
    bcmvn_sub.seurat <- find.pK(sweep.stats_sub.seurat)
    pk = 0.09
    nExp_poi <- round(0.075 * dim(obj.sub)[2])
    obj.sub <- doubletFinder_v3(obj.sub, PCs = 1:10, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) # run DoubletFinder
    
    indexDFInfo <- which(names(obj.sub@meta.data) == sprintf("DF.classifications_0.25_%s_%s",pk,nExp_poi)) # extract DF results column index
    tmpcol <- obj.sub@meta.data[indexDFInfo]
    colnames(tmpcol) <- "DoubletFinder"
    dfcols <- rbind(dfcols,tmpcol)
}

obj <- AddMetaData(obj,dfcols,col.name = 'DoubletFinder') # add metadata to main object

# SCTransform

obj = SCTransform(obj, vars.to.regress = "percent.mito", verbose = FALSE) # run Seurat sctransform method
#s.genes <- cc.genes$s.genes # extract genes associated to S cycle
#g2m.genes <- cc.genes$g2m.genes # extract genes associated to G2M cycle
s.genes <- rownames(obj)[tolower(rownames(obj)) %in% tolower(cc.genes$s.genes)]
g2m.genes <- rownames(obj)[tolower(rownames(obj)) %in% tolower(cc.genes$g2m.genes)]
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident=TRUE) # compute cell cyle scores for all cells

# Dim Red

obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs=100)

nDims = 30
obj <- RunUMAP(obj, dims = 1:nDims) # run UMAP on PCA
obj <- FindNeighbors(obj, dims = 1:nDims) # build knn graph then snn graph
obj <- FindClusters(obj) # Louvain clustering

# Save Seurat object

saveRDS(obj, file = save_seurat_obj_path)

# FindAllMarkers

DefaultAssay(obj) <- "RNA"
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000)
markers.obj <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, features=VariableFeatures(obj))
saveRDS(markers.obj, file= sprintf("%s.markers.rds",gsub(".rds$","",save_seurat_obj_path)))