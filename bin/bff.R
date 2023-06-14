#!/usr/bin/env Rscript

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager",repos = "http://cran.us.r-project.org")
#is_demuxmix<-require("demuxmix")
#if (!require("demuxmix", quietly = TRUE))
#    BiocManager::install("demuxmix")
#  BiocManager::install("demuxmix")
#if (!require("DropletUtils", quietly = TRUE))
#  BiocManager::install("DropletUtils")
#if (!require("preprocessCore", quietly = TRUE))
#  BiocManager::install("preprocessCore")
#if (!require("cellhashR", quietly = TRUE))
#  devtools::install_github(repo = 'bimberlab/cellhashR', ref = 'master', dependencies = TRUE, upgrade = 'always')

library(DropletUtils)
library(Seurat)
library(ggplot2)
library(cowplot)
library(cellhashR)
library(here)
library(dplyr)
library(argparse)
library(tidyverse)

# Create a parser
parser <- ArgumentParser("Parameters for BFF")
parser$add_argument("--fileHto", "Path to file HTO count matrix.")
parser$add_argument("--methods", help='A vector of one or more calling methods to use.', default="bff_raw,bff_cluster")
parser$add_argument("--methodsForConsensus", help='By default, a consensus call will be generated using all methods', default=NULL)
parser$add_argument("--cellbarcodeWhitelist", help='A vector of expected cell barcodes. This allows reporting on the total set of expected barcodes, not just those in the filtered count matrix')
parser$add_argument("--metricsFile", help='If provided, summary metrics will be written to this file.', default="metrics_cell_hash_r.csv")
parser$add_argument("--doTSNE", help='If true, tSNE will be run on the resulting hashing calls after each caller.', default=TRUE)
parser$add_argument("--doHeatmap", help='f true, Seurat::HTOHeatmap will be run on the results of each caller.', default=TRUE)
parser$add_argument("--perCellSaturation", help='An optional dataframe with the columns cellbarcode and saturation',default=NULL)
parser$add_argument("--majorityConsensusThreshold", help='This applies to calculating a consensus call when multiple algorithms are used',default=NULL)
parser$add_argument("--chemistry", help='This string is passed to EstimateMultipletRate. Should be either 10xV2 or 10xV3.', default="10xV3")
parser$add_argument("--callerDisagreementThreshold", help='If provided, the agreement rate will be calculated between each caller and the simple majority call, ignoring discordant and no-call cells.',default=NULL)
parser$add_argument("--rawFeatureMatrixH5", help="filepath to the 10x h5 gene expression counts file")
parser$add_argument("--assignmentOutBff", help="Prefix name for the file containing the output of BFF assignment", type = "character", default = "bff")

parser$add_argument("--outputdir", help='Output directory')
args <- parser$parse_args()

#Open Seurat object containing HTO asssay
args <- parser$parse_args()
if (!endsWith(args$seuratObject, ".rds")){
  seuratObj <- list.files(args$seuratObject, pattern = "\\.rds$", full.names = TRUE)[1]
}else{
  seuratObj <- args$seuratObject
}
#Parameters originally Null
methodsForConsensus <- args$methodsForConsensus
if(is.null(methodsForConsensus)){
  methodsForConsensus <- "NULL"
}
cellbarcodeWhitelist <- args$cellbarcodeWhitelist
if(is.null(cellbarcodeWhitelist)){
  cellbarcodeWhitelist <- "NULL"
}
perCellSaturation <- args$perCellSaturation
if(is.null(perCellSaturation)){
  perCellSaturation <- "NULL"
}
majorityConsensusThreshold <- args$majorityConsensusThreshold
if(is.null(majorityConsensusThreshold)){
  majorityConsensusThreshold <- "NULL"
}
callerDisagreementThreshold <- args$callerDisagreementThreshold
if(is.null(callerDisagreementThreshold)){
  callerDisagreementThreshold <- "NULL"
}

#saving parameters in a dataframe
Argument <- c("HTO-File", "methods", "methodsForConsensus", "cellbarcodeWhitelist", "metricsFile", "perCellSaturation","majorityConsensusThreshold","callerDisagreementThreshold", "doTSNE","doHeatmap","chemistry")
Value <- c(args$fileHto, args$methods, methodsForConsensus, cellbarcodeWhitelist, args$metricsFile, perCellSaturation, majorityConsensusThreshold, allerDisagreementThreshold, args$doTSNE, args$doHeatmap,args$chemistry)
params <- data.frame(Argument, Value)

# Loading Seurat object
counts <- Read10X(data.dir = args$fileHto)

#get methods used for demultiplexing
methods_cell_hash_r <- str_split_1(args$methods, ",")
methods_cell_hash_r <- c(methods_cell_hash_r)

#loading h5 data if available
if(isTRUE('demuxmix' %in% methods_cell_hash_r | 'demuxem' %in% methods_cell_hash_r )){
  cell_hash_R_res <- GenerateCellHashingCalls(barcodeMatrix = counts, methods = methods_cell_hash_r, , doTSNE = args$doTSNE, doHeatmap = args$doHeatmap,methodsForConsensus = args$methodsForConsensus,cellbarcodeWhitelist = args$cellbarcodeWhitelist, metricsFile= args$metricsFile, perCellSaturation= args$perCellSaturation, majorityConsensusThreshold = args$majorityConsensusThreshold, chemistry = args$chemistry, callerDisagreementThreshold = args$callerDisagreementThreshold, rawFeatureMatrixH5 = args$rawFeatureMatrixH5 )
}else{
  cell_hash_R_res <- GenerateCellHashingCalls(barcodeMatrix = counts, methods = methods_cell_hash_r, , doTSNE = args$doTSNE, doHeatmap = args$doHeatmap,methodsForConsensus = args$methodsForConsensus,cellbarcodeWhitelist = args$cellbarcodeWhitelist, metricsFile= args$metricsFile, perCellSaturation= args$perCellSaturation, majorityConsensusThreshold = args$majorityConsensusThreshold, chemistry = args$chemistry, callerDisagreementThreshold = args$callerDisagreementThreshold )
}




write.csv(params, paste0(args$outputdir, "/params.csv"))
write.csv(cell_hash_R_res, paste0(args$outputdir, "/", args$assignmentOutBff, "_assignment_bff.csv"), row.names=FALSE)

