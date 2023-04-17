#!/usr/bin/env Rscript
#Libraries
library(Seurat)
library(argparse)
library(ggplot2)

# Create a parser
parser <- ArgumentParser("Parameters for HTODemux Visualisation")
parser$add_argument("--hashtagPath",help="folder where rds object was saved from the first part of HTODemux")
parser$add_argument("--assay",help="Name of the Hashtag assay HTO by default", default = "HTO")
#Output graphs - Ridge Plot
parser$add_argument("--ridgePlot", help = "Generates a ridge plot from the results, True to generate", default = "TRUE")
parser$add_argument("--ridgeNCol", help = "Number of columns for ridgePlot", default = 3, type = "integer")

#Output graphs - Scatter Feature
parser$add_argument("--featureScatter",help = "Generates a ridge plot from the results, True to generate", default = "TRUE")
parser$add_argument("--scatterFeat1", help = "Feature 1 for Feature Scatter Plot", default = "hto_HTO-A")
parser$add_argument("--scatterFeat2", help = "Feature 2 for Feature Scatter Plot", default = "hto_HTO-B")

#Output graphs - Violin Plot
parser$add_argument("--vlnPlot", help = "Generates a violin plot from the results, True to generate", default = "TRUE")
parser$add_argument("--vlnFeatures", help = "Features to plot (gene expression, metrics, PC scores, anything that can be retreived by FetchData)", default = "nCount_RNA")
parser$add_argument("--vlnLog", help = "plot the feature axis on log scale", action = "store_true")

#Output graphs - tSNE
parser$add_argument("--tSNE", help = "Generate a two dimensional tSNE embedding for HTOs", default = "TRUE")
parser$add_argument("--tSNEIdents", help = "What should we remove from the object (we have Singlet,Doublet and Negative)", default = "Negative")
parser$add_argument("--tSNEInvert", action = "store_true") # TRUE
parser$add_argument("--tSNEVerbose", action = "store_true") # FALSE
parser$add_argument("--tSNEApprox", action = "store_true") # FALSE
parser$add_argument("--tSNEDimMax", help = "max number of donors ",type = "integer", default = 1)
parser$add_argument("--tSNEPerplexity", help = "value for perplexity", type = "integer",  default = 100)

#Output graphs - Heatmap
parser$add_argument("--heatMap", help = "Generate a Heatmap", default = "FALSE")
parser$add_argument("--heatMapNcells", help ="value for number of cells", type = "integer",  default = 500)
parser$add_argument("--outputdir", help='Output directory')

args <- parser$parse_args()

Argument <- c("hashtagPath", "assay", "ridgePlot", "ridgeNCol", "featureScatter", "scatterFeat1", "scatterFeat2", "vlnPlot", "vlnFeatures", "vlnLog", "tSNE", "tSNEIdents", "tSNEInvert", "tSNEVerbose", "tSNEApprox", "tSNEDimMax", "tSNEPerplexity")
Value <- c(args$hashtagPath, args$assay, args$ridgePlot, args$ridgeNCol, args$featureScatter, args$scatterFeat1, args$scatterFeat2, args$vlnPlot, args$vlnFeatures, args$vlnLog, args$tSNE, args$tSNEIdents, args$tSNEInvert, args$tSNEVerbose, args$tSNEApprox, args$tSNEDimMax, args$tSNEPerplexity)

params <- data.frame(Argument, Value)


print(paste0("Loading RDS object in ", args$hashtagPath))
if (!endsWith(args$hashtagPath, ".rds")){
    hash_file <- list.files(args$hashtagPath, pattern = "\\.rds$", full.names = TRUE)[1]
}else{
    hash_file <- args$hashtagPath
}
print(hash_file)
hashtag <-readRDS(hash_file)


# Ridge Plot
# Group cells based on the max HTO signal
if(args$ridgePlot == "TRUE"){
  Idents(hashtag) <- paste0(args$assay, "_maxID")
  RidgePlot(hashtag, assay = args$assay, features = rownames(hashtag[[args$assay]]), ncol = args$ridgeNCol)
  ggsave(paste0(args$outputdir, '/ridge.jpeg'), device = 'jpeg', dpi = 500, height = 10, width = 10)
}

if(args$featureScatter == "TRUE"){
  FeatureScatter(hashtag, feature1 = args$scatterFeat1, feature2 = args$scatterFeat2)
  ggsave(paste0(args$outputdir, '/featureScatter.jpeg'), device = 'jpeg',dpi = 500)
}

if(args$vlnPlot == "TRUE"){
  Idents(hashtag) <- paste0(args$assay, "_classification.global")
  VlnPlot(hashtag, features = args$vlnFeatures, pt.size = 0.1, log = args$vlnLog)
  ggsave(paste0(args$outputdir, '/violinPlot.jpeg'), device = 'jpeg',dpi = 500)
}

if(args$tSNE == "TRUE"){
  hashtag.subset <- subset(hashtag, idents = args$tSNEIdents, invert = args$tSNEInvert)
  DefaultAssay(hashtag.subset) <- args$assay
  hashtag.subset <- ScaleData(hashtag.subset, features = rownames(hashtag.subset),
                                   verbose = args$tSNEVerbose)
  hashtag.subset <- RunPCA(hashtag.subset, features = rownames(hashtag.subset), approx = args$tSNEApprox)
  hashtag.subset <- RunTSNE(hashtag.subset, dims = 1:args$tSNEDimMax, perplexity = args$tSNEPerplexity)
  DimPlot(hashtag.subset)
  ggsave(paste0(args$outputdir, '/tSNE.jpeg'), device = 'jpeg',dpi = 500)
}

if(args$heatMap == "TRUE"){
  HTOHeatmap(hashtag, assay = args$assay, ncells = args$heatMapNcells)
  ggsave(paste0(args$outputdir, '/heatMap.jpeg'), device = 'jpeg',dpi = 500)
}

write.csv(params, paste0(args$outputdir, "/visual_params.csv"))





