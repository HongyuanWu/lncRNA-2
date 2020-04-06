# Title     : 16 Days clustering
# Objective : Extract lncRNA sequences from the different populations at 3 different times and pool them together by population and timepoint
# Created by: DAVID
# Created on: 01/04/2020

library(tidyverse)
library(Seurat)
library(patchwork)

#Read the data
#-------------
day16.data <-Read10X(data.dir = "data/FGF8plus_day16/filtered_gene_bc_matrices/GRCh38")
day16 <- CreateSeuratObject(counts = day16.data, project="lncRNA", min.cells = 5, min.features = 100)
#min cells takes features present in at leas that many cells.
#min features includes cells with at least that many features

#Pre-Processing
#--------------
day16[["percent.mt"]]<-PercentageFeatureSet(day16, pattern = "^MT-")#calculate mt% to clean it later

#visualize features
VlnPlot(day16, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#check nFeature plot, pretty wide


plot1 <- FeatureScatter(day16, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(day16, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


day16 <- subset(day16, subset = nFeature_RNA > 900 & nFeature_RNA < 4000 & percent.mt < 3)
#subset cells with at least 500 genes and max 4000 and less than 3% mt contamination

#Normalize the data
#------------------

day16 <- NormalizeData(day16,normalization.method = "LogNormalize", scale.factor = 10000)
#Is logarithmic normalization good enough(?)

#Feature selection
#-----------------

sag