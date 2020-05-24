
# Analysis of long non-coding RNA expression from single cell datasets

# Objective
The objective of this proejct is to design a pipeline/tools to analyze single cell RNA expression data using R and packages such as [Seurat](https://satijalab.org/seurat/) and [DeSeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) for especific diferential expression analysis and [Tidyverse](https://www.tidyverse.org/) as a general toolbox for different plots. (Note that ggplot2 is necesary for Seurat and its contained in Tidyverse).
## Procedure
After runing the single cell sequencing protocol on our cell pupulations per manufacturer's protocol using 10x Genomics Single Cell 3' chips in a nextseq machine. The raw basecalls were translated into fastq format using cellranger mkfastq that uses bcl2fastq provided by illumina then the resulting fasta files were processed using the cellranger count V 3.0 provided by 10x.

The different fucntions created for downstream analysis can be found on "source/AnalysisFunctions.R". Detailed explanation in the `Source section`

Then the data was analyzed using different combinations of parameters using the "Mass analysis Script.R" found in Source, this script allows the user to enter different parameters to run the analysis to determine which parameters produce results that are biologically coherent.

After inspection of the results the specific parameters were determined ( using 1000 top variable features on which perform the statistical analysis and 0.1 resolution to run the clustering). These parameters were used to run the "source/FinalAnalysis.R", this script will also plot the different gene lists over the data to annotate the cell populations present on each cluster.

# Source
## Final Analysis
This script will reproduce my clustering results on my datasets. Please note that since these results are not published yet they cannot be shared.
## Analysis functions
This file calls the different R packages necesary for the analysis and contains the fucntions necesary to find the different clusters in this type of data.
The functions that contains an `saveIMG` and `outdir` arguments will save tyhe plot generated in the specified outdir after trying to create it, if `saveIMG` is set to `FALSE` (`TRUE` by default) there is no need to specify and `outdir`.
- **OpenData**: returns the seurat object from a dorectory where all 3 matrix files are located (barcodes, features, matrix). This function will aslo print violin plots for number of features, feature counts and mitochondrial gene presence in our data and FeatureScatter plots 
 Usage: `DataObject <- OpenData(dir=PathToData, project.name = Name,outdir = PathToOutput )`
 Note that it is necesary to subset the data after looking at the previous plots.
 `DataObject <- subset(DataObject , subset = nFeature_RNA > 1100 & nFeature_RNA < 4000 & percent.mt < 3)`
- **NormalizeAndScale**: This function takes a seurat object and a number of variable features, and runs the data normalization using seurat default values ("LogNormalize" as method and a scale factor of 10000). It will return a seurat object with 
Usage: `DataObject <- NormalizeAndScale(data.object = DataObject, nFeatures = 500, outdir = PathToOutput )`
Note: If you want to use different parameters for the normalization use the standard seurat functions, an example can be found [here](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html).
- **LinearAnalysis**: This fucniton will take a seurat object as input and will return a seurat object with the added results. It will run Principal Component Analysis (PCA) and a JackStraw test using the specified number of dimensions `dims` and replciates `JackStrawReplicates`. It will also score the JackStraw test. The arguments `cells` and `balance` affect the heatmap printed after the analysis. This function will also print a Jack Straw Plot with the evaluated components and their pvalues, an elbow plot and a scatter plot of the first 2 components.
Usage: `DataObject <- LinearAnalysis(DataObject ,dims = 40, saveImg = TRUE,outdir = PathToOutput )`
- **DimThreshold**: This function will take a seurat object and a pvalue and will return a list with all the comonents with a JackStrawScore < pvalue.
Usage: `dims <- DimThreshold(DataObject ,pval = 0.000001)`
Note that the seurat object needs to contain the results created in the previopus step.
- **Cluster**: This function takes a seurat object, a number of dimensions (user determined or the resulting vector from the `DimThreshold` function and resolutions for UMAP and tSNE. It will find the different neighbors in our data and it will also run UMAP and tSNE clustering with the specified resolutions and dimensions. It will return the data object with the added results and it will also print the clustering graphs in the specified outdir.
Usage: ` DataObject <- Cluster(data.object = DataObject , dims = dims, UMAPres = 0.1, tSNEREs = 0.15, outdir = PathToOutput)`
- **ExtractMarkers **: This fuction takes a seurat object and 2 parameter thresholds and will return a tible with the top 2 most variable markers per clusters.
Usage: `Object.markers<- ExtractMarkers(data.object = DataObject, min.pct = 0.2, logfc.threshold = 0.4)`
- **GetGeneList**: This fucntion will take the path to a directory in which there are .txt files with marker names (1 per line) and will return a list with one element per file. The element name will be the original file name and the file content will be the list content.
Usage: `cells <- GetGeneList("data/GeneLists/")`
- **DotPlotGenes**: This function will take a seurat object and a list with features to plot on our data. (Note that this is meant to worl with a list with the same structure as the one created by the `GetGeneList` function). It will output a DotPlot for each element of the list plotting the features inside the list.
Usage:`DotPlotGenes(data.object = DataObject, geneList = cells, outdir = PathToOutput)`
- **VlnPlotGenes**: This function will take a seurat object and a list with features to plot on our data, the argument `ArrayOf` determines how many Violin plots are present per image. (Note that this is meant to worl with a list with the same structure as the one created by the `GetGeneList` function). It will output a ViolinPlot for each element of the list plotting the features inside the list.
Usage: `VlnPlotGenes(data.object = DataObject, ArrayOf = 9, geneList = cells, outdir = PathToOutput)`

## Mass analysis
This script is meant to be used a brute force analysis tool to test different parameters for our analysis.
The user has to specify a list of paths to the diferent datasets, a vector with the desired number of most variable counts (`featCounts`), a vector with the desired resolutions (`Res`), a diretory from which write the results and a cell feature list.
It will create a file tree with branches at each different level (Dataset, Variable Features and resolution) with the branched folder named as the varying parameter.




# Data

In this directory you can find a subdirectory with the gene markers used for my results.