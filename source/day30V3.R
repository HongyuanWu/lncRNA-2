# Title     : 30 Days clustering pipeline 3.0
# Objective : Extract lncRNA sequences from the different populations at 3 different times and pool them together by population and timepoint
# Created by: DAVID
# Created on: 01/04/2020

library(tidyverse)
library(Seurat)
library(patchwork)

#Read the data
#-------------
day30.data <-Read10X(data.dir = "data/V3/FGF8plus_day30_V3/filtered_feature_bc_matrix")
day30 <- CreateSeuratObject(counts = day30.data, project="lncRNA", min.cells = 5, min.features = 100)
#min cells takes features present in at leas that many cells.
#min features includes cells with at least that many features

#Pre-Processing
#--------------
day30[["percent.mt"]]<-PercentageFeatureSet(day30, pattern = "^MT-")#calculate mt% to clean it later

#visualize features
VlnPlot(day30, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#check nFeature plot, pretty wide


plot1 <- FeatureScatter(day30, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(day30, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2



day30 <- subset(day30, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 3)
#subset cells with at least 1000 genes and max 4000 and less than 3% mt contamination


#Normalize the data
#------------------

day30 <- NormalizeData(day30,normalization.method = "LogNormalize", scale.factor = 10000)
#Is logarithmic normalization good enough(?)


#Feature selection
#-----------------
day30 <- FindVariableFeatures(day30, selection.method = "vst", nfeatures = 500)
#taking 3000 most variable features

top20 <- head(VariableFeatures(day30),20)
plot1 <- VariableFeaturePlot(day30)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2


#Scaling the data
#----------------
all.genes <- rownames(day30)
day30 <- ScaleData(day30, features = all.genes)


# Linear dimensional reduction
#-----------------------------
day30 <- RunPCA(day30, features = VariableFeatures(object = day30))

VizDimLoadings(day30, dims = 1:4, reduction = "pca")

DimPlot(day30, reduction = "pca")
DimHeatmap(day30, dims = 1, cells = 1000, balanced = TRUE)
DimHeatmap(day30, dims = 1:21, cells = 1000, balanced = TRUE)


#Dimensionality
#--------------

day30 <- JackStraw(day30,num.replicate = 100)
day30 <- ScoreJackStraw(day30, dims = 1:20)
JackStrawPlot(day30, dims = 1:20)
ElbowPlot(day30) #I'll take 3

#Clustering
#----------
day30 <- FindNeighbors(day30,dims = 1:3)
day30 <- FindClusters(day30, resolution = 0.05)

#Non-linear reduction
#--------------------
day30 <- RunUMAP(day30, dims = 1:3)
day30 <-RunTSNE(day30, dims = 1:3)
DimPlot(day30, reduction = "umap")#With 0.1 resolution looks better
DimPlot(day30, reduction = "tsne") #With 0.2 resolution looks better
#saveRDS(day30, file = "output/day30.rds")


#Find clusters
#-------------
day30.markers <- FindAllMarkers(day30, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
#cluster
day30.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#Sort the clusters

#-------------
head(day30.markers,20)
VlnPlot(day30, features = c("TH"))

top3<-day30.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

FeaturePlot(day30, features = top3$gene, reduction = "tsne")
FeaturePlot(day30, features = top3$gene, reduction = "umap")
#--------------

#HEatmap
top15 <- day30.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
DoHeatmap(day30, features = top15$gene) + NoLegend()
