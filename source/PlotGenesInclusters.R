

source("source/AnalysisFunctions.R")


datadirs <-  list("data/V3/FGF8plus_day16_V3/filtered_feature_bc_matrix")
  
"data/V3/FGF8plus_day30_V3/filtered_feature_bc_matrix","data/V3/FGF8plus_day60_V3/filtered_feature_bc_matrix")

cells <- GetGeneList("data/GeneLists/")

for (datadir in datadirs){
  project <- strsplit(strsplit(datadir, split = "/")[[1]][3],split = "_")[[1]][2]
  results <- paste0("output/NameClusters1/",project,"/")
  dataset <- OpenData(dir = datadir, project.name = project, saveImg = TRUE, outdir = results)
  dataset <- NormalizeAndScale(dataset, nfeatures = 500, saveImg = TRUE,outdir = results)
  dataset <- LinearAnalysis(dataset,dims = 40, saveImg = TRUE,outdir = results)
  dims <- DimThreshold(dataset,pval = 0.000001)
  results <- paste0("output/NameClusters1/",project,"/")
  dataset <- Cluster(dataset,dims = dims, UMAPres = 0.1, tSNEREs = 0.1, saveImg = TRUE,outdir = results)
  DotPlotGenes(geneList = cells, data.object = dataset, outdir = results )
  
}




DimThreshold <- function(data.object, pval){
  scores <- as.data.frame((dataset[["pca"]]@"jackstraw")$overall.p.values) 
  ValidScores<-subset(scores, Score<pval)
  result<-list(many=nrow(ValidScores), Valid=as.vector(ValidScores$PC))
  return(result)
}


DimThreshold(dataset,0.000001)

