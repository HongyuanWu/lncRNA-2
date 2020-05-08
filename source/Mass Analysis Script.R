

source("source/AnalysisFunctions.R")

datadirs <-  list("data/V3/FGF8plus_day16_V3/filtered_feature_bc_matrix", "data/V3/FGF8plus_day30_V3/filtered_feature_bc_matrix","data/V3/FGF8plus_day60_V3/filtered_feature_bc_matrix")

featCounts <- c(500)

Res <-c(0.1)

cells <- GetGeneList("data/GeneLists/")

for (datadir in datadirs){
  project <- strsplit(strsplit(datadir, split = "/")[[1]][3],split = "_")[[1]][2]
  out1 <-paste0("output/archive/",project,"/")
  dataset <- OpenData(dir=datadir, project.name = paste0("lncRNA_",project),outdir = out1 ) 
  for (features in featCounts){
    out2<-paste0(out1, features,"/")
    dataset <- NormalizeAndScale(data.object = dataset,nfeatures = features, outdir = out2)
    dataset <- LinearAnalysis(data.object = dataset, outdir = out2)
    for (dims in c(10)){#Add Pval threshoild
      out3 <- paste0(out2, dims,"_dims/")
      for(res in Res){
        out4 <- paste0(out3,res,"_Res/")
      dataset <- Cluster(data.object = dataset, dims = dims, UMAPres = res, tSNEREs = res, outdir = out4)
        out5 <- paste0(out4,"cells/")
        DotPlotGenes(data.object = dataset, geneList = cells, outdir = out5)
        VlnPlotGenes(data.object = dataset, ArrayOf = 9, geneList = cells, outdir = out5)
      }
    }
  }
  rm(dataset)
}




VlnPlotGenes(data.object = dataset, ArrayOf = 9, geneList = cells, outdir = out5)
