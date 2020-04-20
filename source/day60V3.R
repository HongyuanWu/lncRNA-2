# Title     : 30 Days clustering pipeline 3.0
# Objective : Extract lncRNA sequences from the different populations at 3 different times and pool them together by population and timepoint
# Created by: DAVID
# Created on: 01/04/2020

library(tidyverse)
library(Seurat)
library(patchwork)

#Read the data
#-------------
day60.data <-Read10X(data.dir = "data/V3/FGF8plus_day60_V3/filtered_feature_bc_matrix")
day60 <- CreateSeuratObject(counts = day60.data, project="lncRNA", min.cells = 5, min.features = 100)
#min cells takes features present in at leas that many cells.
#min features includes cells with at least that many features

#Pre-Processing
#--------------
day60[["percent.mt"]]<-PercentageFeatureSet(day60, pattern = "^MT-")#calculate mt% to clean it later

#visualize features
VlnPlot(day60, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#check nFeature plot, pretty wide


plot1 <- FeatureScatter(day60, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(day60, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2



day60 <- subset(day60, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 3)
#subset cells with at least 1000 genes and max 4000 and less than 3% mt contamination


#Normalize the data
#------------------

day60 <- NormalizeData(day60,normalization.method = "LogNormalize", scale.factor = 10000)
#Is logarithmic normalization good enough(?)


#Feature selection
#-----------------
day60 <- FindVariableFeatures(day60, selection.method = "vst", nfeatures = 2000)
#taking 3000 most variable features

top20 <- head(VariableFeatures(day60),20)
plot1 <- VariableFeaturePlot(day60)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2


#Scaling the data
#----------------
all.genes <- rownames(day60)
day60 <- ScaleData(day60, features = all.genes)


# Linear dimensional reduction
#-----------------------------
day60 <- RunPCA(day60, features = VariableFeatures(object = day60))

VizDimLoadings(day60, dims = 1:4, reduction = "pca")

DimPlot(day60, reduction = "pca")
DimHeatmap(day60, dims = 1, cells = 1000, balanced = TRUE)
DimHeatmap(day60, dims = 1:21, cells = 1000, balanced = TRUE)


#Dimensionality
#--------------

day60 <- JackStraw(day60,num.replicate = 100)
day60 <- ScoreJackStraw(day60, dims = 1:20)
JackStrawPlot(day60, dims = 1:20)
ElbowPlot(day60) #I'll take 4

#Clustering
#----------
day60 <- FindNeighbors(day60,dims = 1:4)
day60 <- FindClusters(day60, resolution = 0.25)

#Non-linear reduction
#--------------------
day60 <- RunUMAP(day60, dims = 1:4)
day60 <-RunTSNE(day60, dims = 1:4)
DimPlot(day60, reduction = "umap")#With 0.1 resolution looks better
DimPlot(day60, reduction = "tsne") #With 0.2 resolution looks better
#saveRDS(day60, file = "output/day60.rds")


#Find clusters
#-------------
day60.markers <- FindAllMarkers(day60, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
#cluster
day60.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#Sort the clusters

#-------------
head(day60.markers,20)
VlnPlot(day60, features = c("MALAT1"))

top3<-day60.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

FeaturePlot(day60, features = top3$gene, reduction = "tsne")
FeaturePlot(day60, features = top3$gene, reduction = "umap")
#--------------

#HEatmap
top15 <- day60.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
DoHeatmap(day60, features = top15$gene) + NoLegend()
