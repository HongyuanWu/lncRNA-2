source("source/AnalysisFunctions.R") #Get the fucntions and load the packages

datadirs <-  list("data/V3/FGF8plus_day16_V3/filtered_feature_bc_matrix", "data/V3/FGF8plus_day30_V3/filtered_feature_bc_matrix","data/V3/FGF8plus_day60_V3/filtered_feature_bc_matrix") #Where are my datafiles

featCounts <- c(500,1000,1500, 2000) #feature counts to work with
Res <-c(0.1) # Resolutions to work with
out <- "output/archive_3K_0.000001/" # Root output dir
cells <- GetGeneList("data/GeneLists/") #cell markers

for (datadir in datadirs){ #For each dataset
  project <- strsplit(strsplit(datadir, split = "/")[[1]][3],split = "_")[[1]][2] # this just cuts dayxx                                                                                      from the path
  out1 <-paste0(out,project,"/") #Subpath for the dataset
  dataset <- OpenData(dir = datadir, project.name = paste0("lncRNA_",project),outdir = out1 ) #It will create the folder
  for (features in featCounts){
    out2<-paste0(out1, features,"/")# subpath for the featureCounts
    dataset <- subset(dataset, subset = nFeature_RNA > 1000 & nFeature_RNA < 3000 & percent.mt < 3)
    dataset <- NormalizeAndScale(data.object = dataset,nFeatures = features, outdir = out2) #It will create the directory
    dataset <- LinearAnalysis(dataset,dims = 40, saveImg = TRUE,outdir = out2)
    dims <- DimThreshold(dataset,pval = 0.000001)
    out3 <- paste0(out2, length(dims),"_dims/") #Path for the number of dimensions selected
    for(res in Res){
      out4 <- paste0(out3,res,"_Res/") #Path for the resoution selected
      dataset <- Cluster(data.object = dataset, dims = dims, UMAPres = res, tSNEREs = res, outdir = out4)
      out5 <- paste0(out4,"cells/")
      DotPlotGenes(data.object = dataset, geneList = cells, outdir = out5) # Dot plto the features
      VlnPlotGenes(data.object = dataset, ArrayOf = 9, geneList = cells, outdir = out5)
    }
    }
  rm(dataset) #Clean the seurat object and repeat
  }
  






