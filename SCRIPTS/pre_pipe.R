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
s.genes <- rownames(obj)[tolower(rownames(obj)) %in% tolower(cc.genes$s.genes)]
g2m.genes <- rownames(obj)[tolower(rownames(obj)) %in% tolower(cc.genes$g2m.genes)]
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident=TRUE) # compute cell cyle scores for all cells
obj <- SCTransform(obj, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('percent.mito', 'S.Score', 'G2M.Score')) # normalize again but this time including also the cell cycle scores

# Dim Red
obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs=100)
nDims = 30
obj <- RunUMAP(obj, dims = 1:nDims) # run UMAP on PCA
obj <- FindNeighbors(obj, dims = 1:nDims) # build knn graph then snn graph
obj <- FindClusters(obj,resolution = 1) # Louvain clustering

# Save Seurat object

saveRDS(obj, file = save_seurat_obj_path)

# FindAllMarkers

DefaultAssay(obj) <- "RNA"
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000)
markers.obj <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, features=VariableFeatures(obj))
saveRDS(markers.obj, file= sprintf("%s.markers.rds",gsub(".rds$","",save_seurat_obj_path)))


############
# QC figures
############

input_seurat_obj_path <- save_seurat_obj_path
#------------------------------------------#
input_markers_obj_path <- sprintf("%s.markers.rds",gsub(".rds$","",input_seurat_obj_path))
save_pdf_path          <- sprintf("%s.figures.pdf",gsub(".rds$","",input_seurat_obj_path))
#------------------------------------------#

#------------------------------------------#
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(scales)
#------------------------------------------#

#------------------------------------------#
cluster.contingency.table <- function(df, var.name, main.var.name = "cluster") {
    df2            <- data.frame(table(df[, main.var.name ], df[, var.name ]))
    colnames(df2)  <- c("cluster","metadata","n")
    df2            <- df2[order(as.numeric(gsub("^c","",df2$cluster)),decreasing=TRUE),]
    df2$cluster    <- factor( df2$cluster, levels= rev(levels(df2$cluster)) )
    return(df2)
}
cluster.coordinates <- function(df, var.name, var.dim1, var.dim2) {
    clust.coordinates <- matrix(ncol=2,nrow=0)
    colnames(clust.coordinates) <- c("dim1","dim2")
    for (clust in unique(df[,var.name])) {
        indexes <- which( df[,var.name] == clust )
    
        clust.coordinates <- rbind(clust.coordinates , matrix( apply( df[indexes, c(var.dim1, var.dim2) , drop=FALSE ],2, median),ncol=2,dimnames=list(clust,c("dim.1","dim.2"))))
    }
    clust.coordinates         <- as.data.frame(clust.coordinates)
    clust.coordinates$cluster <- rownames(clust.coordinates)
    clust.coordinates         <- clust.coordinates[ levels(df[,var.name]), ]

    return(clust.coordinates)
}
inc <- function(x)
{
 eval.parent(substitute(x <- x + 1))
}
scatter.plot <- function(df, var.name, var.dim1, var.dim2, label.coord.data = NULL , title = NULL, gradient.colors = NULL , colors = NULL , with.legend = "no") {
    #--------------------------------------#
    new.df <- df[,c(var.dim1,var.dim2,var.name)]
    colnames(new.df) <- c("dim1","dim2","metadata")
    #--------------------------------------#

    #--------------------------------------#
    p <- ggplot( new.df , aes( x=dim1,y=dim2,color=metadata ) ) +
        geom_point( shape=20 , size=0.5 , alpha=0.5 ) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 10, hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank() ,
            panel.grid.major.x = element_line(colour = "white"),
            panel.grid.minor.x = element_blank(),       
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),       
            panel.background = element_rect(fill = "white", colour="white"),        
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border = element_blank()
        ) 
    if ( with.legend == "no" ) {
        p <- p+ theme( legend.position = 'none' )
    }       
    if ( is.null(title) == FALSE ) {
        p <- p + ggtitle(title)
    }
    if ( is.null(gradient.colors) == FALSE ) {
        p <- p + scale_color_gradientn(colors=gradient.colors)
    }
    if ( is.null(colors) == FALSE ) {
        p <- p + scale_color_manual(values=colors)
    }
    p <- p + xlab("Dim. 1") + ylab("Dim. 2")

    if ( is.null(label.coord.data) == FALSE ) {
        p <- p + geom_text(data=label.coord.data, aes(x=dim1,y=dim2,label=cluster), color="black", size=3)
    }
    #--------------------------------------#
    return( p )
    #--------------------------------------#
}
spot.plot <- function(df, title = NULL) {
    p <- ggplot(data=df, aes(gene, cluster)) +
        geom_point(aes(size = pct, color = scale.expr )) +
        scale_color_gradientn( colors= brewer.pal(n = 8, name = "Reds")) +
        theme_bw() +
        theme(
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6),
            axis.text.y=element_text(size=10),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            legend.text = element_text(size = 10),
            legend.title=element_blank()
        ) 

    if ( is.null(title) == FALSE ) {
        p <- p + ggtitle(title)
    }
    return(p)
}
#------------------------------------------#

#------------------------------------------#
seurat.obj  <- readRDS(input_seurat_obj_path)
markers.obj <- readRDS(input_markers_obj_path)
#------------------------------------------#

#------------------------------------------#
Idents(seurat.obj) <- seurat.obj@meta.data$seurat_clusters
DefaultAssay(seurat.obj)      <- "RNA"
#------------------------------------------#

#------------------------------------------#
df       <- seurat.obj@meta.data
df$umap1 <- seurat.obj@reductions$umap@cell.embeddings[rownames(seurat.obj@meta.data),1]
df$umap2 <- seurat.obj@reductions$umap@cell.embeddings[rownames(seurat.obj@meta.data),2]
df$cells <- rownames(seurat.obj@meta.data)
#------------------------------------------#

#------------------------------------------#
df$outlier <- TRUE
df$outlier[ which(df$DoubletFinder == "Singlet" & df$scateroutlier == FALSE) ] <- "FALSE"

df$scateroutlier <- factor(df$scateroutlier, levels=c(FALSE,TRUE))
df$DoubletFinder <- factor(df$DoubletFinder, levels=c("Singlet","Doublet"))
df$outlier       <- factor(df$outlier, levels=c(FALSE,TRUE))
df$Phase         <- factor(df$Phase, levels=c("G1","S","G2M"))
#------------------------------------------#

#------------------------------------------#
clust.coord  <- cluster.coordinates(df, "seurat_clusters", "umap1", "umap2")
sample.coord <- cluster.coordinates(df, "orig.ident"     , "umap1", "umap2")
#------------------------------------------#

#-----------------------------------------#
#--- Identify top marker genes -----------#
#-----------------------------------------#
sel.markers.obj1 <- markers.obj[ !markers.obj$gene %in% unique(markers.obj$gene[grep("^[Mm][Tt]-|^R[Pp][SsLl]|^HIST|^Hist",markers.obj$gene)]), , drop=FALSE]
sel.markers.obj1 <- sel.markers.obj1[ !tolower(sel.markers.obj1$gene) %in% tolower(c(cc.genes$s.genes, cc.genes$g2m.genes)), , drop=FALSE]
sel.markers.obj1 <- sel.markers.obj1[ which( sel.markers.obj1$pct.1 >= 0.25 & sel.markers.obj1$pct.1 > sel.markers.obj1$pct.2 & sel.markers.obj1$avg_log2FC > 0 ), , drop=FALSE]
sel.markers.obj1 <- sel.markers.obj1[ which( p.adjust(sel.markers.obj1$p_val,method="BH") <= 0.05 ), , drop=FALSE]
sel.markers.obj1 <- sel.markers.obj1[ order(sel.markers.obj1$pct.1 - sel.markers.obj1$pct.2,decreasing=TRUE),, drop=FALSE]

sel.markers.obj2 <- sel.markers.obj1[ sel.markers.obj1$gene %in% names(which( table(sel.markers.obj1$gene) <= 3)), , drop=FALSE]
#-----------------------------------------#

#-----------------------------------------#
n.max <- 5
sel.genes1 <- c()
sel.genes2 <- c()
for (clust in unique(markers.obj$cluster)) {
    clust.indexes1 <- which(sel.markers.obj1$cluster == clust)
    clust.indexes2 <- which(sel.markers.obj2$cluster == clust)

    n1 <- length(clust.indexes1)
    if ( n1 > 0) {
        if ( n1 > n.max) {n1 <- n.max}
        sel.genes1 <- c(sel.genes1, sel.markers.obj1[clust.indexes1[1:n1],,drop=FALSE]$gene)
    }
    n2 <- length(clust.indexes2)
    if ( n2 > 0) {
        if ( n2 > n.max ) {n2 <- n.max}
        sel.genes2 <- c(sel.genes2, sel.markers.obj2[clust.indexes2[1:n2],,drop=FALSE]$gene)
    }
}
sel.genes1 <- intersect(unique(sel.genes1[which( is.na(sel.genes1) == FALSE)]), rownames(seurat.obj))
sel.genes2 <- intersect(unique(sel.genes2[which( is.na(sel.genes2) == FALSE)]), rownames(seurat.obj))
sel.genes3 <- sort(rownames(seurat.obj)[ toupper(rownames(seurat.obj)) %in% c("HBE1","HBB","CD34","PECAM1","CD68","CD14","PTPRC","POU5F1","NANOG","DAZL","SYCP3","SRY","SOX9","AMH","CLU","INSL3","CYP17A1","CYP11A1","UPK3B","KRT19","WT1","GATA4","PDGFRA","PDGFRB","TCF21","SST","PLAU","CXCL14","MARCH3","PAX8","IFI44","IFIT3","OAS3","DDX58","PNOC","CDH6","CD24","FOXL2","LHX2","IRX3","CA9")])
sel.genes  <- unique(c(sel.genes1,sel.genes2, sel.genes3))
sel.genes  <- sort(sel.genes[ which( is.na(sel.genes) == FALSE ) ])
#--------------------------------------#
    
        
#--------------------------------------#
clusts <- levels(seurat.obj@meta.data$seurat_clusters)  
#--------------------------------------#
pct.m  <- matrix(0,nrow=length(sel.genes),ncol=length(clusts), dimnames=list(sel.genes, clusts))
expr.m <- matrix(0,nrow=length(sel.genes),ncol=length(clusts), dimnames=list(sel.genes, clusts))
for (clust in clusts) {
    clust.cells      <- rownames(seurat.obj@meta.data)[which(seurat.obj@meta.data$seurat_clusters == clust)]
    expr.m[,clust] <- apply(seurat.obj@assays$RNA@data[sel.genes,clust.cells,drop=FALSE],1, mean)
    pct.m[,clust]  <- apply(seurat.obj@assays$RNA@data[sel.genes,clust.cells,drop=FALSE],1, function(x) length(which(x > 0))) * 100 / length(clust.cells)
}
expr.m <- round( expr.m, digit = 1)
pct.m  <- round( pct.m , digit = 1)

indexes <- which( apply( expr.m, 1, sd ) > 0 )
expr.m <- expr.m[indexes,,drop=FALSE]
pct.m  <- pct.m[indexes,,drop=FALSE]
#--------------------------------------#
scale.expr.m <- t(apply(expr.m,1,rescale, to=c(0,1)))
scale.expr.m <- round( scale.expr.m, digit = 3)
#--------------------------------------#
df2 <- cbind( melt(expr.m),  melt(scale.expr.m)[,3,drop=FALSE], melt(pct.m)[,3,drop=FALSE] )
colnames(df2) <- c("gene","cluster","expr","scale.expr","pct")
#--------------------------------------#
df2$cluster <- factor(df2$cluster, levels=colnames(scale.expr.m)[hclust(dist(t(scale.expr.m)))$order])
df2$gene    <- factor(df2$gene,    levels=rownames(scale.expr.m)[hclust(dist(scale.expr.m))$order])
#------------------------------------------#

#------------------------------------------#
i  <- 0
gs1 <- list()

gs1[[inc(i)]] <- AugmentPlot(scatter.plot(df, "seurat_clusters", "umap1","umap2",label.coord.data=clust.coord, title="clusters", with.legend="yes"),dpi=200)
gs1[[inc(i)]] <- AugmentPlot(scatter.plot(df, "orig.ident"     , "umap1","umap2",label.coord.data=clust.coord, title="samples", with.legend="yes"),dpi=200)

for (category in c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo")) {
    gs1[[inc(i)]] <- AugmentPlot(scatter.plot(df, category, "umap1","umap2",title=category, gradient.colors=c("lightgrey","red","darkred"),with.legend="no"),dpi=200)
}
for (category in c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo")) {
    sub.df <- df[,c("seurat_clusters",category)]
    colnames(sub.df) <- c("cluster","metadata")
    gs1[[inc(i)]] <- ggplot(sub.df, aes(x=cluster, y=metadata , fill=cluster)) + geom_boxplot(outlier.shape=NA) + coord_flip() + theme_bw() + scale_x_discrete(limits=rev) + theme( legend.position = 'none' , axis.title.x = element_blank(),axis.title.y = element_blank())
}
for (category in c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo")) {
    sub.df <- df[,c("orig.ident",category)]
    colnames(sub.df) <- c("cluster","metadata")
    gs1[[inc(i)]] <- ggplot(sub.df, aes(x=cluster, y=metadata , fill=cluster)) + geom_boxplot(outlier.shape=NA) + coord_flip() + theme_bw() + scale_x_discrete(limits=rev) + theme( legend.position = 'none' , axis.title.x = element_blank(),axis.title.y = element_blank())
}

for (category in c("scateroutlier","DoubletFinder","outlier")) {
    gs1[[inc(i)]] <- AugmentPlot(scatter.plot(df, category, "umap1","umap2",title=category, colors=c("green","red"),with.legend="no"),dpi=200)
    
    if (category == "scateroutlier") {
        gs1[[inc(i)]] <- AugmentPlot(scatter.plot(df, category, "scateroutlierPC1","scateroutlierPC2",title=category, colors=c("green","red"),with.legend="no"),dpi=200)
    }
}
for (category in c("scateroutlier","DoubletFinder","outlier")) {
    gs1[[inc(i)]] <- ggplot(cluster.contingency.table(df,category,"seurat_clusters"), aes(x=cluster,y=n,fill=metadata)) + geom_bar(position="fill", stat="identity") +  scale_fill_manual(values=colorRampPalette(c("#DDDDDD", "#222222"))(length(levels(df[[category]])))) + theme_bw() + coord_flip() + theme( plot.title = element_text(size = 10, hjust = 0.5), legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank())
}
for (category in c("scateroutlier","DoubletFinder","outlier")) {
    gs1[[inc(i)]] <- ggplot(cluster.contingency.table(df,category,"orig.ident"), aes(x=cluster,y=n,fill=metadata)) + geom_bar(position="fill", stat="identity") +   scale_fill_manual(values=colorRampPalette(c("#DDDDDD", "#222222"))(length(levels(df[[category]])))) + theme_bw() + coord_flip() + theme( plot.title = element_text(size = 10, hjust = 0.5), legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank())
}
#------------------------------------------#

#------------------------------------------#
i   <- 0
gs2 <- list()
gs2[[inc(i)]] <- AugmentPlot(scatter.plot(df, "seurat_clusters", "umap1","umap2",label.coord.data=clust.coord, title="clusters", with.legend="yes"),dpi=200)
gs2[[inc(i)]] <- AugmentPlot(scatter.plot(df, "orig.ident"     , "umap1","umap2",label.coord.data=clust.coord, title="samples", with.legend="yes"),dpi=200)

for (category in c("S.Score","G2M.Score","Phase")) {
    if (category != "Phase") {
        gs2[[inc(i)]] <- AugmentPlot(scatter.plot(df, category, "umap1","umap2",title=category, gradient.colors=c("lightgrey","red","darkred"),with.legend="no"),dpi=200)
    } else {
        gs2[[inc(i)]] <- AugmentPlot(scatter.plot(df, category, "umap1","umap2",title=category, colors=c("lightgrey","yellow","red"),with.legend="no"),dpi=200)
    }
}

gs2[[inc(i)]] <- ggplot(cluster.contingency.table(df,"seurat_clusters","orig.ident"), aes(x=cluster,y=n,fill=metadata)) + geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + theme( plot.title = element_text(size = 10, hjust = 0.5), legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank())

for (category in c("S.Score","G2M.Score","Phase")) {
    if (category != "Phase") {
        sub.df <- df[,c("seurat_clusters",category)]
        colnames(sub.df) <- c("cluster","metadata")

        gs2[[inc(i)]] <- ggplot(sub.df, aes(x=cluster, y=metadata , fill=cluster)) + geom_boxplot(outlier.shape=NA) + coord_flip() + theme_bw() + scale_x_discrete(limits=rev) + theme( legend.position = 'none' , axis.title.x = element_blank(),axis.title.y = element_blank())
    } else {
        gs2[[inc(i)]] <- ggplot(cluster.contingency.table(df,category,"seurat_clusters"), aes(x=cluster,y=n,fill=metadata)) + geom_bar(position="fill", stat="identity") +  scale_fill_manual(values=colorRampPalette(c("#DDDDDD", "#222222"))(length(levels(df[[category]])))) + theme_bw() + coord_flip() + theme( plot.title = element_text(size = 10, hjust = 0.5), legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank())
    }
}

gs2[[inc(i)]] <- ggplot(cluster.contingency.table(df,"orig.ident","seurat_clusters"), aes(x=cluster,y=n,fill=metadata)) + geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + theme( plot.title = element_text(size = 10, hjust = 0.5), legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank())

for (category in c("S.Score","G2M.Score","Phase")) {
    if (category != "Phase") {
        sub.df <- df[,c("orig.ident",category)]
        colnames(sub.df) <- c("cluster","metadata")

        gs2[[inc(i)]] <- ggplot(sub.df, aes(x=cluster, y=metadata , fill=cluster)) + geom_boxplot(outlier.shape=NA) + coord_flip() + theme_bw() + scale_x_discrete(limits=rev) + theme( legend.position = 'none' , axis.title.x = element_blank(),axis.title.y = element_blank())
    } else {
        gs2[[inc(i)]] <- ggplot(cluster.contingency.table(df,category,"orig.ident"), aes(x=cluster,y=n,fill=metadata)) + geom_bar(position="fill", stat="identity") +   scale_fill_manual(values=colorRampPalette(c("#DDDDDD", "#222222"))(length(levels(df[[category]])))) + theme_bw() + coord_flip() + theme( plot.title = element_text(size = 10, hjust = 0.5), legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank())
    }
}

gs2[[inc(i)]] <- spot.plot( df2 , title = "Find markers"        ) + theme(  plot.title = element_text(size = 10, hjust = 0.5), legend.position = 'none')
#------------------------------------------#

#------------------------------------------#
categories <- c("orig.ident","nCount_RNA","nFeature_RNA","percent.mito","percent.ribo","scateroutlier","scateroutlierPC1","scateroutlierPC2","DoubletFinder","nCount_SCT","nFeature_SCT","S.Score","G2M.Score","Phase","old.ident","seurat_clusters", colnames(seurat.obj@meta.data)[grep("^SCT_snn_res",colnames(seurat.obj@meta.data))])
other.categories <- setdiff( colnames(seurat.obj@meta.data), categories )

i   <- 0
gs3 <- list()
for (category in other.categories) {

    if ( is.numeric(seurat.obj@meta.data[[category]]) == TRUE ) {
        gs3[[inc(i)]] <- AugmentPlot(scatter.plot(df, category, "umap1","umap2",title=category, gradient.colors=c("lightgrey","red","darkred"),with.legend="no"),dpi=200)

        sub.df <- df[,c("seurat_clusters",category)]
        colnames(sub.df) <- c("cluster","metadata")
        gs3[[inc(i)]] <- ggplot(sub.df, aes(x=cluster, y=metadata , fill=cluster)) + geom_boxplot(outlier.shape=NA) + coord_flip() + theme_bw() + scale_x_discrete(limits=rev) + theme( legend.position = 'none' , axis.title.x = element_blank(),axis.title.y = element_blank())

        sub.df <- df[,c("orig.ident",category)]
        colnames(sub.df) <- c("cluster","metadata")
        gs3[[inc(i)]] <- ggplot(sub.df, aes(x=cluster, y=metadata , fill=cluster)) + geom_boxplot(outlier.shape=NA) + coord_flip() + theme_bw() + scale_x_discrete(limits=rev) + theme( legend.position = 'none' , axis.title.x = element_blank(),axis.title.y = element_blank())

    } else {
        gs3[[inc(i)]] <- AugmentPlot(scatter.plot(df, category, "umap1","umap2",title=category, with.legend="yes"),dpi=200)
        gs3[[inc(i)]] <- ggplot(cluster.contingency.table(df,category,"seurat_clusters"), aes(x=cluster,y=n,fill=metadata)) + geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + theme( plot.title = element_text(size = 10, hjust = 0.5), legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank())
        gs3[[inc(i)]] <- ggplot(cluster.contingency.table(df,category,"orig.ident")     , aes(x=cluster,y=n,fill=metadata)) + geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + theme( plot.title = element_text(size = 10, hjust = 0.5), legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank())
    
    }
}
#------------------------------------------#

#------------------------------------------#
pdf(save_pdf_path,height=28,width=14)
grid.arrange(grobs=gs1,ncol=4,nrow=8,layout_matrix = rbind(c(1, 1, 2, 2 ),
                                                   c(1, 1, 2, 2 ),
                                                   c(3, 4, 5, 6 ),
                                                   c(7, 8, 9, 10),
                                                   c(11,12,13,14),
                                                   c(15,16,17,18),
                                                   c(19,NA,20,21),
                                                   c(22,NA,23,24)))

grid.arrange(grobs=gs2,ncol=4,nrow=8,layout_matrix = rbind(c(1, 1, 2, 2 ),
                                                   c(1, 1, 2, 2 ),
                                                   c(NA,3, 4, 5 ),
                                                   c(6, 7, 8, 9),
                                                   c(10,11,12,13),
                                                   c(14,14,14,14),
                                                   c(14,14,14,14),
                                                   c(14,14,14,14)))
                                                   
if (length(gs3) > 0) {
    grid.arrange(grobs=gs3,ncol=3)
}                                                  
   
dev.off()
#------------------------------------------#

quit(save="no")
