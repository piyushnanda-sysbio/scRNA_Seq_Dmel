#Single cell RNA Sequencing for Ahmad's analysis

#install.packages('Seurat')
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)

setwd('/Users/piyushnanda/Documents/scRNASeq/')

data_normal='/Users/piyushnanda/Documents/scRNASeq/Normal/Sample'
data_starved='/Users/piyushnanda/Documents/scRNASeq/Starved/Sample'

#Normal Cells

expression_matrix <- Read10X(data.dir = data_normal)
normal = CreateSeuratObject(counts = expression_matrix)

#normal[["percent.mt"]] <- PercentageFeatureSet(normal, pattern = "^MT-")

VlnPlot(normal, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

normal <- NormalizeData(normal, normalization.method = "LogNormalize", scale.factor = 10000)

normal <- FindVariableFeatures(normal, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(normal), 50)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(normal)
plot_norm <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot_norm

#Starved Cells

expression_matrix <- Read10X(data.dir = data_starved)
starved = CreateSeuratObject(counts = expression_matrix)

#normal[["percent.mt"]] <- PercentageFeatureSet(normal, pattern = "^MT-")

VlnPlot(starved, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

starved <- NormalizeData(starved, normalization.method = "LogNormalize", scale.factor = 10000)

starved <- FindVariableFeatures(starved, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(starved), 50)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(starved)
plot_starved <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot_starved

plot_norm+plot_starved

all.genes <- rownames(normal)
normal <- ScaleData(normal, features = all.genes)

all.genes <- rownames(starved)
starved <- ScaleData(starved, features = all.genes)

#PCA

normal <- RunPCA(normal, features = VariableFeatures(object = normal))
print(normal[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(normal, dims = 1:2, reduction = "pca")

DimPlot(normal, reduction = "pca")

ElbowPlot(normal)

normal <- FindNeighbors(normal, dims = 1:10)
normal <- FindClusters(normal, resolution = 0.5)

head(Idents(normal), 100)

normal <- RunUMAP(normal, dims = 1:10)
DimPlot(normal, reduction = "umap")

cluster2.markers <- FindMarkers(normal, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 20)

saveRDS(normal, file = "normal.rds")

FeaturePlot(normal, features = c("Ilp2","Akh",'mnd','Ilp5','Ilp3','slif','CG9413','Ilp6'))

normal.markers <- FindAllMarkers(normal, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

normal.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(normal, features = top10$gene) + NoLegend()

##For Starved

starved <- RunPCA(starved, features = VariableFeatures(object = starved))
print(starved[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(starved, dims = 1:2, reduction = "pca")

DimPlot(starved, reduction = "pca")

ElbowPlot(starved)

starved <- FindNeighbors(starved, dims = 1:15)
starved <- FindClusters(starved, resolution = 0.5)

head(Idents(starved), 5)

starved <- RunUMAP(starved, dims = 1:15)
DimPlot(starved, reduction = "umap")

cluster2.markers <- FindMarkers(starved, ident.1 = 6, min.pct = 0.25)
head(cluster2.markers, n = 5)

saveRDS(starved, file = "starved.rds")

FeaturePlot(starved, features = c("Ilp2","Akh",'mnd','Ilp5','Ilp3','slif','CG9413','Ilp6'))

VlnPlot(starved,features = c('Ilp2'))


starved.markers <- FindAllMarkers(starved, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

starved.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(starved, features = top10$gene) + NoLegend()

clust0_d<-subset(x=starved,idents=0)
clust0<-subset(x=starved,idents=0,subset=Akh>1)
expr<-AverageExpression(clust0)

cells_wakh<-subset(x=starved,subset=Akh>1)
cells_wilp<-subset(x=starved,subset=Ilp2>5)

