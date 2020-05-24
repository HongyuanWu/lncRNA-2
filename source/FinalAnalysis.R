
source("source/AnalysisFunctions.R") # load functions and packages

cells <- GetGeneList("data/GeneLists/") #cell markers
featCounts <- c(1000) #Feature counts
Res <-c(0.1)
out <- "output/Results_Pval_0.000001/" #Root for results

dir16 <-"data/V3/FGF8plus_day16_V3/filtered_feature_bc_matrix"
project <- strsplit(strsplit(dir16, split = "/")[[1]][3],split = "_")[[1]][2]
out1 <-paste0(out,project,"/")
day16 <- OpenData(dir=dir16, project.name = paste0("lncRNA_",project),outdir = out1 )
write(c("BeforeQC ", ncol(day16)), file = paste0(out1,"CelllogQC.txt"))
for (features in featCounts){
  out2<-paste0(out1, features,"/")
  day16 <- subset(day16, subset = nFeature_RNA > 1100 & nFeature_RNA < 4000 & percent.mt < 3)
  day16 <- NormalizeAndScale(data.object = day16, nFeatures = features, outdir = out2)
  write(c("AfterQC ", ncol(day16)), file = paste0(out1,"celllogQC.txt"), append = TRUE)
  day16 <- LinearAnalysis(day16,dims = 40, saveImg = TRUE,outdir = out2)
  dims <- DimThreshold(day16,pval = 0.000001)
  out3 <- paste0(out2, length(dims),"_dims/")
  for(res in Res){
    out4 <- paste0(out3,res,"_Res/")
    day16 <- Cluster(data.object = day16, dims = dims, UMAPres = res, tSNEREs = res, outdir = out4)
    out5 <- paste0(out4,"cells/")
    day16.markers<- ExtractMarkers(day16, min.pct = 0.2, logfc.threshold = 0.4)
    DotPlotGenes(data.object = day16, geneList = cells, outdir = out5)
    VlnPlotGenes(data.object = day16, ArrayOf = 9, geneList = cells, outdir = out5)
  }
}

#After looking at the Vplot and DotPlot name the clusters
new.cluster.ids16 <- c("InDevelopment","Neurons","Oligodendrocytes+PanNeurons+Extraneuronal","FBLC")
names(new.cluster.ids16) <- levels (day16)
day16<-RenameIdents(day16, new.cluster.ids16)#update the object


png(filename=paste0(out1,"ClusterDimPlot.png"))
print(DimPlot(day16,reduction="umap", label=TRUE))
dev.off()

saveRDS(day16, file=paste0(out1,"day16_1000nFeat_0.1Res.rds")) #Save the object, you dotn want to run this each time

write.table(table(Idents(day16)),file=paste0(out1,"ResultLogDay16.txt"), row.names = FALSE)
write.table(prop.table(table(Idents(day16))),file=paste0(out1,"ResultLogDay16.txt"), row.names = FALSE, append = TRUE) #Some population stats

write.table(WhichCells(day16, idents = "Neurons"),file=paste0(out1,"NeuronBarcodesDay16.txt"),row.names = FALSE, col.names = FALSE) # print the barcodes to run lncRNA enrichments


#----------------------------------------


dir30 <- "data/V3/FGF8plus_day30_V3/filtered_feature_bc_matrix"
project <- strsplit(strsplit(dir30, split = "/")[[1]][3],split = "_")[[1]][2]
out1 <-paste0(out,project,"/")
day30 <- OpenData(dir=dir30, project.name = paste0("lncRNA_",project),outdir = out1 ) 
write(c("BeforeQC ", ncol(day30)), file = paste0(out1,"CelllogQC.txt"))
for (features in featCounts){
  out2<-paste0(out1, features,"/")
  day30 <- subset(day30, subset = nFeature_RNA > 1100 & nFeature_RNA < 4000 & percent.mt < 3)
  day30 <- NormalizeAndScale(data.object = day30, nFeatures = features, outdir = out2)
  write(c("AfterQC ", ncol(day30)), file = paste0(out1,"celllogQC.txt"), append = TRUE)
  day30 <- LinearAnalysis(day30,dims = 40, saveImg = TRUE,outdir = out2)
  dims <- DimThreshold(day30,pval = 0.000001)
  out3 <- paste0(out2, length(dims),"_dims/")
  for(res in Res){
    out4 <- paste0(out3,res,"_Res/")
    day30 <- Cluster(data.object = day30, dims = dims, UMAPres = res, tSNEREs = res, outdir = out4)
    out5 <- paste0(out4,"cells/")
    day30.markers<- ExtractMarkers(day30, min.pct = 0.2, logfc.threshold = 0.4)
    DotPlotGenes(data.object = day30, geneList = cells, outdir = out5)
    VlnPlotGenes(data.object = day30, ArrayOf = 9, geneList = cells, outdir = out5)
  }
}

#After looking at the Vplot and DotPlot name the clusters
new.cluster.ids30 <- c("Neurons","Oligodendrocytes+Progenitors","InDevelopment","FBLC")
names(new.cluster.ids30) <- levels (day30)
day30<-RenameIdents(day30, new.cluster.ids30)

png(filename=paste0(out1,"ClusterDimPlot.png"))
print(DimPlot(day30,reduction="umap", label=TRUE))
dev.off()

saveRDS(day30, file=paste0(out1,"day30_1000nFeat_0.1Res.rds"))

write.table(table(Idents(day30)),file=paste0(out1,"ResultLogDay30.txt"), row.names = FALSE)
write.table(prop.table(table(Idents(day30))),file=paste0(out1,"ResultLogDay30.txt"), row.names = FALSE, append = TRUE)

write.table(WhichCells(day30, idents = "Neurons"),file=paste0(out1,"NeuronBarcodesDay30.txt"),row.names = FALSE, col.names = FALSE)


#----------------------------------------


dir60 <- "data/V3/FGF8plus_day60_V3/filtered_feature_bc_matrix"
project <- strsplit(strsplit(dir60, split = "/")[[1]][3],split = "_")[[1]][2]
out1 <-paste0(out,project,"/")
day60 <- OpenData(dir=dir60, project.name = paste0("lncRNA_",project),outdir = out1 )
write(c("BeforeQC ", ncol(day60)), file = paste0(out1,"CelllogQC.txt"))
for (features in featCounts){
  out2<-paste0(out1, features,"/")
  day60 <- subset(day60, subset = nFeature_RNA > 1100 & nFeature_RNA < 4000 & percent.mt < 3)
  day60 <- NormalizeAndScale(data.object = day60, nFeatures = features, outdir = out2)
  write(c("AfterQC ", ncol(day60)), file = paste0(out1,"celllogQC.txt"), append = TRUE)
  day60 <- LinearAnalysis(day60,dims = 40, saveImg = TRUE,outdir = out2)
  dims <- DimThreshold(day60,pval = 0.000001)
  out3 <- paste0(out2, length(dims),"_dims/")
  for(res in Res){
    out4 <- paste0(out3,res,"_Res/")
    day60 <- Cluster(data.object = day60, dims = dims, UMAPres = res, tSNEREs = res, outdir = out4)
    out5 <- paste0(out4,"cells/")
    day60.markers<- ExtractMarkers(day60, min.pct = 0.2, logfc.threshold = 0.4)
    DotPlotGenes(data.object = day60, geneList = cells, outdir = out5)
    VlnPlotGenes(data.object = day60, ArrayOf = 9, geneList = cells, outdir = out5)
  }
}

#After looking at the Vplot and DotPlot name the clusters
new.cluster.ids60 <- c("Neurons","InDevelopment","Oligodendrocytes","FBLC")
names(new.cluster.ids60) <- levels (day60)
day60<-RenameIdents(day60, new.cluster.ids60)

png(filename=paste0(out1,"ClusterDimPlot.png"))
print(DimPlot(day60,reduction="umap", label=TRUE))
dev.off()

saveRDS(day60, file=paste0(out1,"day60_1000nFeat_0.1Res.rds"))

write.table(table(Idents(day60)),file=paste0(out1,"ResultLogDay60.txt"), row.names = FALSE)
write.table(prop.table(table(Idents(day60))),file=paste0(out1,"ResultLogDay60.txt"), row.names = FALSE, append = TRUE)

write.table(WhichCells(day60, idents = "Neurons"),file=paste0(out1,"NeuronBarcodesDay60.txt"),row.names = FALSE, col.names = FALSE)

#######################

populations <- read.csv("output/Results_Pval_0.000001/cell_cluster2.csv", sep = ";") # Read the csv file created from the ident proportions from each dataset

populations$Cluster<-fct_relevel(populations$Cluster,"Neurons", "In Development", "FBLC","Oligodendrocytes") #Reorder so it looks nicer

#Plot em by day
ggplot(populations, aes(Timepoint, Population, fill=Cluster, group=Cluster)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_brewer(palette = "Set1") 

#Plot em by cluster
ggplot(populations, aes(Cluster, Population, fill=Timepoint, group=Timepoint)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_brewer(palette = "Set1") 



