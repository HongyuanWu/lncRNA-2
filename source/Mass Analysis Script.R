

source("source/AnalysisFunctions.R")

datadirs <-  list("data/V3/FGF8plus_day16_V3/filtered_feature_bc_matrix", "data/V3/FGF8plus_day30_V3/filtered_feature_bc_matrix","data/V3/FGF8plus_day60_V3/filtered_feature_bc_matrix")

featCounts <- c(1000)

Res <-c(0.3)

cells <- GetGeneList("data/GeneLists/")

for (datadir in datadirs){
  project <- strsplit(strsplit(datadir, split = "/")[[1]][3],split = "_")[[1]][2]
  out1 <-paste0("output/archive_PvalThreshold 0.000001/",project,"/")
  dataset <- OpenData(dir=datadir, project.name = paste0("lncRNA_",project),outdir = out1 ) 
  for (features in featCounts){
    out2<-paste0(out1, features,"/")
    dataset <- NormalizeAndScale(data.object = dataset,nfeatures = features, outdir = out2)
    dataset <- LinearAnalysis(dataset,dims = 40, saveImg = TRUE,outdir = out2)
    dims <- DimThreshold(dataset,pval = 0.000001)
    out3 <- paste0(out2, length(dims),"_dims/")
    for(res in Res){
      out4 <- paste0(out3,res,"_Res/")
      dataset <- Cluster(data.object = dataset, dims = dims, UMAPres = res, tSNEREs = res, outdir = out4)
      out5 <- paste0(out4,"cells/")
      DotPlotGenes(data.object = dataset, geneList = cells, outdir = out5)
      VlnPlotGenes(data.object = dataset, ArrayOf = 9, geneList = cells, outdir = out5)
    }
    }
  rm(dataset)
  }
  






