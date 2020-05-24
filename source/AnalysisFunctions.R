
library(tidyverse)
library(Seurat)
library(patchwork)

OpenData <- function(dir="data/",
                     project.name="project",
                     min.cels=5,
                     min.feat=100,
                     outdir="output/",
                     saveImg=TRUE){

  data <- Read10X(data.dir=dir)

  data.object<-CreateSeuratObject(counts=data,
                                  project=project.name,
                                  min.cells=min.cels,
                                  min.features=min.feat)
  data.object[["percent.mt"]]<-PercentageFeatureSet(data.object, pattern = "^MT-") #% od mitochondrial contamination

  Vplot<-VlnPlot(data.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
  #create the Vplot with dots and without
  VplotND<-VlnPlot(data.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, pt.size = 0)
  

  plot1 <- FeatureScatter(data.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(data.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- plot1 + plot2
  
  if (saveImg==TRUE){
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE) #Create the directories and save the files
    png(filename=paste0(outdir,"FeatViolinPlot.png"))
    print(Vplot)
    dev.off()
    
    png(filename=paste0(outdir,"FeatViolinPlotND.png"))
    print(VplotND)
    dev.off()
    
    png(filename=paste0(outdir,"FeatScatterPlot.png"))
    print(plot3)
    dev.off()
    }

  return(data.object)
}


NormalizeAndScale <- function(data.object,
                              nFeatures = 500,
                              outdir = "output/",
                              saveImg = TRUE){

  data.object <- NormalizeData(data.object,
                               normalization.method = "LogNormalize",
                               scale.factor = 10000) #Normalize the expression data

  data.object <- FindVariableFeatures(data.object,
                                      selection.method = "vst",
                                      nfeatures = nFeatures) #find the top variable features

  all.genes <- rownames(data.object)
  data.object <- ScaleData(data.object, features = all.genes)#Scale the data to a mean expression of 0 and variance 1.


  top20<-head(VariableFeatures(data.object),20)
  plot1 <- VariableFeaturePlot(data.object)
  plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
  
  if(saveImg==TRUE){
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE) #Create the output dir and print the image
    png(filename=paste0( outdir , "VarFeatPlot.png"))
    print(plot2)
    dev.off()
    }
  return(data.object)
}



LinearAnalysis <- function (data.object,
                            dims=35,
                            cells=1000,
                            balance=TRUE,
                            JackStrawReplicates=100,
                            outdir="output/",
                            saveImg=TRUE){


  data.object <- RunPCA(data.object, features = VariableFeatures(object = data.object)) #PCA analysis
  data.object <- JackStraw(data.object, dims= dims, num.replicate = JackStrawReplicates) #JAckstraw
  data.object <- ScoreJackStraw(data.object, dims = 1:dims)

  # Create the plots
  heatm <- DimHeatmap(data.object, dims = 1:(dims/4), cells = 1000, balanced = TRUE)
  dimplot <- DimPlot(data.object, reduction = "pca")
  JSPlot <- JackStrawPlot(data.object, dims = 1:dims)
  elbow <- ElbowPlot(data.object)

  if(saveImg==TRUE){ #print the plots
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    
    png(filename=paste0(outdir,"HeatMap.png"))
    print(heatm)
    dev.off()
    
  
    png(filename=paste0(outdir,"DimplotPCA.png"))
    print(dimplot)
    dev.off()
  
    png(filename=paste0(outdir,"JackStrawPlot",dims,"dims.png"))
    print(JSPlot)
    dev.off()
  
    png(filename=paste0(outdir,"ElbowPlotPCA.png"))
    print(elbow)
    dev.off()
    }

  return(data.object)
}


DimThreshold <- function(data.object, pval){
  scores <- as.data.frame((data.object[["pca"]]@"jackstraw")$overall.p.values) # Access the seurat object and look for the scores
  ValidScores<-subset(scores, Score<pval) # Select which pass the threshold
  result<-as.vector(ValidScores$PC) #Return the PC number in a vector this way it can be 1,2,3,5,8...
  return(as.vector(result))
}

Cluster <-  function (data.object,
                      dims=10,
                      UMAPres=0.1,
                      tSNEREs=0.1,
                      outdir="output/",
                      saveImg=TRUE){
  
  data.object <- FindNeighbors(data.object,dims = dims)

  if(UMAPres==tSNEREs){ #just to save processing time
    data.object <- FindClusters(data.object, resolution = UMAPres) #Find clusters

    data.object <- RunUMAP(data.object, dims = dims) #RunUMAP
    UMAPplot <- DimPlot(data.object, reduction = "umap") #create the plot

    data.object <- RunTSNE(data.object, dims = dims) #RuntSNE
    tSNEPlot <- DimPlot(data.object, reduction = "tsne") #create the plot
  }
  else{ #If the resolutions are different run each one separately
    data.object <- FindClusters(data.object, resolution = UMAPres) 
    data.object <- RunUMAP(data.object, dims = dims)
    UMAPplot <- DimPlot(data.object, reduction = "umap")

    data.object <- FindClusters(data.object, resolution = tSNEREs)
    data.object <- RunTSNE(data.object, dims = dims)
    tSNEPlot <- DimPlot(data.object, reduction = "tsne")

  }
  
  if(saveImg==TRUE){ #save the images
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    png(filename=paste0(outdir,"UMAPPlot_",length(dims),"dims_",UMAPres,"Res.png"))
    print(UMAPplot)
    dev.off()
  
    png(filename=paste0(outdir,"tSNEplot_",length(dims),"dims_",tSNEREs,"Res.png"))
    print(tSNEPlot)
    dev.off()
    }

  return(data.object)
}

ExtractMarkers <- function (data.object, min.pct=0.1, logfc.threshold=0.2){ #extract top 2 markers/ cluster
  markers <- FindAllMarkers(data.object, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold)
  markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
  return(markers)
}




GetGeneList <- function(directory){
  FileList <- list.files(directory, pattern = ".txt") # ls the directory
  ResList <- list() #placeholder
  for (file in FileList){
    path <- paste0(directory,file)
    Filedata <- as.list(read.delim(path, sep = "\n")) #Extract the info
    ResList <- append(ResList, Filedata) #Add the info to the result object
    
    
  }
  
  return(ResList)
}

DotPlotGenes <- function ( geneList, data.object, 
                           outdir = "output/",
                           saveImg=TRUE){
  for (cell in names(geneList)){ #Iterate for each element (cell type)
    GenesToPlot <- cells[cell][[1]] #Extract the genes
    
    img<- DotPlot(data.object, features = GenesToPlot) + coord_flip()
    
    if(saveImg==TRUE){
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      png(filename=paste0(outdir,cell,".png"))
      print(img)
      dev.off()
    }
  }
}

VlnPlotGenes <- function ( geneList, data.object,ArrayOf=6, 
                           outdir="output/", 
                           saveImg=TRUE){
  for (cell in names(geneList)){ #Iterate for each element (cell type)
    GenesToPlot <- cells[cell][[1]] #Extract the genes
    #Sliding window to plot the features
    lower=1
    upper=ArrayOf
    Identifier <- 0 #Identifier to orther the result files
    if(length(GenesToPlot)>ArrayOf){ #If there are more markers than we want in 1 file
      while(upper<=length(GenesToPlot)){
        plotgenes <- GenesToPlot[lower:upper] #Select the ones we want
        img<- VlnPlot(data.object, features = plotgenes)
        
        lower <- lower + ArrayOf #move the window
        upper <- upper + ArrayOf
        Identifier <- Identifier + 1 
        
        if(saveImg==TRUE){#Save the images
          dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
          png(filename=paste0(outdir,cell,"_VPlot_",Identifier,".png"))
          print(img)
          dev.off()
        }
      }
      
      plotgenes <- GenesToPlot[lower:length(GenesToPlot)] # add the last iteration image
      img<- VlnPlot(data.object, features = plotgenes)
      if(saveImg==TRUE){
        dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
        png(filename=paste0(outdir,cell,"_VPlot_",Identifier,".png"))
        print(img)
        dev.off()
      }
    }
      
    else{ # In case we have less markers thatn desired plots per file
    img<- VlnPlot(data.object, features = GenesToPlot)

    if(saveImg==TRUE){
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      png(filename=paste0(outdir,cell,"_VPlot.png"))
      print(img)
      dev.off()
    }
    }

    
  }
}

  
  
  

