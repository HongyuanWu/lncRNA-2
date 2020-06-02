library("DESeq2")
library("pheatmap")
library("tidyverse")
library("RColorBrewer")


# Input matrixes creation
data <- read.table(file = "data/SCells_LncRNA_s1_CombinedSamples.txt", header = TRUE)
timepoints <- c("Day_16_D","Day_30_D","Day_60_D","Day_16_N","Day_30_N","Day_60_N")
colnames(data)[7:12] <- timepoints
rownames(data)<-data[,1]
count_data <- data[,c(7:12)]

coldata <- data.frame(state=c("InDev","InDev","InDev","Neurons","Neurons","Neurons"),
                      Timepoint=c("Day_16","Day_30","Day_60"), 
                      row.names = timepoints)

#Deseq object creation
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~Timepoint)
levels(dds$Timepoint)
dds$Timepoint <- relevel(dds$Timepoint, ref = "Day_16") #Stablish day 16 as reference

dds <- DESeq(dds)

res16_30 <- results(dds, name = "Timepoint_Day_30_vs_Day_16")
res16_60 <- results(dds, name = "Timepoint_Day_60_vs_Day_16")


## VISZUALIZE
outdir <- "output/DE/"
res <- results(dds, alpha = 0.05)
res<-res[order(res$padj),]
head(res,10)

sink(file = (paste0(outdir,"Summary.txt")))
summary(res)
sink()

# Expression versus significance
png(paste0(outdir, "MAPlot_Pval_0.05.png"))
plotMA(dds, ylim=c(-2,2),  main="MA Plot")
dev.off()


# Data transformations for visualizations
rld <- rlogTransformation(dds,blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)


# heatmap countmatrix
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:30]
df<- as.data.frame(colData(dds)[,c("Timepoint","state")])
top50 <- subset(res, pvalue<0.05)
countTop50 <- subset(counts(dds), rownames(counts(dds)) %in% rownames(top50))[1:25]



top25pval <- res[order(res$padj),] [1:25,]
png(paste0(outdir,"Heatmap_top25byPvalVSD2.png"))
pheatmap(assay(vsd)[rownames(top25pval),], cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
dev.off()


pheatmap(countTop25)

png(paste0(outdir,"Heatmap_top25byPvalVSD.png"))
pheatmap(assay(vsd)[countTop50,], cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
dev.off()

png(paste0(outdir,"Heatmap_countsVSD.png"))
pheatmap(assay(vsd)[select,], cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
dev.off()

png(paste0(outdir,"Heatmap_countsRLD.png"))
pheatmap(assay(rld)[select,], cluster_rows=FALSE,show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
dev.off()

# Heatmap of similarity between replicates
distVSD <- dist(t(assay(vsd)))
matrix <- as.matrix(distVSD)
rownames(matrix) <- paste(vsd$Timepoint,vsd$state, sep = "-")
colnames(matrix) <- paste(vsd$Timepoint,vsd$state, sep = "-")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

png(paste0(outdir,"Heatmap_DistancesNT.png"))
pheatmap(matrix,clustering_distance_rows = distVSD, 
         clustering_distance_cols = distVSD,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = TRUE,show_colnames = TRUE, fontsize = 15,
         color = hmcol, main = "Distance Matrix")
dev.off()

png(paste0(outdir,"Heatmap_Distances.png"))
pheatmap(matrix,clustering_distance_rows = distVSD, 
         clustering_distance_cols = distVSD,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE,show_colnames = TRUE, 
         color = hmcol, main = "Distance Matrix")
dev.off()

# PCA plot
png(paste0(outdir,"PCA.png"))
print(plotPCA(rld, intgroup=c("Timepoint")))
dev.off()
