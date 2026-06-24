#!/usr/bin/env Rscript

################################################
################################################
## Fucntions                                  ##
################################################
################################################

# Helper function for NULL condition
string_to_null <- function(x, val = "null") if (x == val) NULL else x
null_to_string <- function(x, val = "NULL") if (is.null(x)) val else x

################################################
################################################
## USE PARAMETERS FROM NEXTFLOW               ##
################################################
################################################

# cast parameters from nextflow
seurat_object <- '$seurat_object'
assay <- '$assay'
ridgePlot <- as.logical('$ridgePlot')
ridgeNCol <- as.numeric('$ridgeNCol')
featureScatter <- as.logical('$featureScatter')
scatterFeat1 <- string_to_null('$scatterFeat1')
scatterFeat2 <- string_to_null('$scatterFeat2')
vlnPlot <- as.logical('$vlnPlot')
vlnFeatures <- '$vlnFeatures'
vlnLog <- as.logical('$vlnLog')
tSNE <- as.logical('$tSNE')
tSNEIdents <- '$tSNEIdents'
tSNEInvert <- as.logical('$tSNEInvert')
tSNEVerbose <- as.logical('$tSNEVerbose')
tSNEApprox <- as.logical('$tSNEApprox')
tSNEDimMax <- as.numeric('$tSNEDimMax')
tSNEPerplexity <- as.numeric('$tSNEPerplexity')
heatMap <- as.logical('$heatMap')
heatMapNcells <- as.numeric('$heatMapNcells')
prefix <- '$prefix'

# check if the file exists
if (! file.exists(seurat_object)){
    stop(paste0(seurat_object, ' is not a valid file'))
}

################################################
################################################
## Load libraries                             ##
################################################
################################################

library(Seurat)
library(ggplot2)

################################################
################################################
## Main Process                               ##
################################################
################################################

hashtag <- readRDS(seurat_object)

# Ridge Plot
if (ridgePlot) {
  print("Generating ridge plot...")
  Idents(hashtag) <- paste0(assay, "_maxID")
  RidgePlot(hashtag, assay = assay, features = rownames(hashtag[[assay]]), ncol = ridgeNCol)
  ggsave(paste0(prefix, "_ridge_htodemux.jpeg"), device = "jpeg", dpi = 500)
}

# Feature Scatter Plot
if (featureScatter) {
  if (is.null(scatterFeat1) || is.null(scatterFeat2)) {
    available_features <- rownames(hashtag[[assay]])
    if (length(available_features) >= 2) {
      scatterFeat1 <- available_features[1]
      scatterFeat2 <- available_features[2]
    } else {
      stop("Error: Not enough features available for scatter plot")
    }
  }
  FeatureScatter(hashtag, feature1 = scatterFeat1, feature2 = scatterFeat2)
  ggsave(paste0(prefix, "_featureScatter_htodemux.jpeg"), device = "jpeg", dpi = 500)
}

# Violin Plot
if (vlnPlot) {
  print("Generating violin plot...")
  Idents(hashtag) <- paste0(assay, "_classification.global")
  VlnPlot(hashtag, features = vlnFeatures, pt.size = 0.1, log = vlnLog)
  ggsave(paste0(prefix, "_violinPlot_htodemux.jpeg"), device = "jpeg", dpi = 500)
}

# tSNE Plot
if (tSNE) {
  print("Generating tSNE plot...")
  if (tSNEIdents %in% levels(Idents(hashtag))) {
    hashtag.subset <- subset(hashtag, idents = tSNEIdents, invert = tSNEInvert)
  } else {
    hashtag.subset <- hashtag
  }
  DefaultAssay(hashtag.subset) <- assay
  hashtag.subset <- ScaleData(hashtag.subset,
    features = rownames(hashtag.subset),
    verbose = tSNEVerbose
  )
  hashtag.subset <- RunPCA(hashtag.subset, features = rownames(hashtag.subset), approx = tSNEApprox)
  hashtag.subset <- RunTSNE(hashtag.subset, dims = 1:tSNEDimMax, perplexity = tSNEPerplexity, check_duplicates = FALSE)
  DimPlot(hashtag.subset)
  ggsave(paste0(prefix, "_tSNE_htodemux.jpeg"), device = "jpeg", dpi = 500)
}

# Heatmap
if (heatMap) {
  print("Generating heatmap...")
  HTOHeatmap(hashtag, assay = assay, ncells = heatMapNcells)
  ggsave(paste0(prefix, "_heatMap_htodemux.jpeg"), device = "jpeg", dpi = 500)
}

################################################
################################################
## SAVING RESULTS                             ##
################################################
################################################

# Save parameters
Argument <- c(
  "seurat_object",
  "assay",
  "ridgePlot",
  "ridgeNCol",
  "featureScatter",
  "scatterFeat1",
  "scatterFeat2",
  "vlnPlot",
  "vlnFeatures",
  "vlnLog",
  "tSNE",
  "tSNEIdents",
  "tSNEInvert",
  "tSNEVerbose",
  "tSNEApprox",
  "tSNEDimMax",
  "tSNEPerplexity",
  "heatMap",
  "heatMapNcells"
)

Value <- c(
  seurat_object,
  assay,
  ridgePlot,
  ridgeNCol,
  featureScatter,
  null_to_string(scatterFeat1),
  null_to_string(scatterFeat2),
  vlnPlot,
  vlnFeatures,
  vlnLog,
  tSNE,
  tSNEIdents,
  tSNEInvert,
  tSNEVerbose,
  tSNEApprox,
  tSNEDimMax,
  tSNEPerplexity,
  heatMap,
  heatMapNcells
)

params <- data.frame(Argument, Value)
write.csv(params, paste0(prefix, "_visual_params_htodemux.csv"))

################################################
################################################
## SAVE VERSIONS                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
seurat.version <- as.character(packageVersion('Seurat'))
ggplot2.version <- as.character(packageVersion('ggplot2'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    r-seurat:', seurat.version),
        paste('    r-ggplot2:', ggplot2.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
