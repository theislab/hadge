#!/usr/bin/env Rscript

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org")
is_demuxmix<-require("demuxmix")
if (!require("demuxmix", quietly = TRUE))
    BiocManager::install("demuxmix")
  BiocManager::install("demuxmix")
if (!require("DropletUtils", quietly = TRUE))
  BiocManager::install("DropletUtils")
if (!require("preprocessCore", quietly = TRUE))
  BiocManager::install("preprocessCore")
if (!require("cellhashR", quietly = TRUE))
  devtools::install_github(repo = 'bimberlab/cellhashR', ref = 'master', dependencies = TRUE, upgrade = 'always')

library(DropletUtils)
library(Seurat)
library(ggplot2)
library(cowplot)
library(cellhashR)
library(here)
library(dplyr)
library(argparse)

# Create a parser
parser <- ArgumentParser("Parameters for BFF")
parser$add_argument("--seuratObject", help = "Seurat object. Assumes that the hash tag oligo (HTO) data has been added and normalized.")
parser$add_argument("--assay", help='Assay name')
parser$add_argument("--methods", help='A vector of one or more calling methods to use.')
parser$add_argument("--methodsForConsensus", help='By default, a consensus call will be generated using all methods')
parser$add_argument("--cellbarcodeWhitelist", help='A vector of expected cell barcodes. This allows reporting on the total set of expected barcodes, not just those in the filtered count matrix')
parser$add_argument("--metricsFile", help='If provided, summary metrics will be written to this file.')
parser$add_argument("--doTSNE", help='If true, tSNE will be run on the resulting hashing calls after each caller.')
parser$add_argument("--doHeatmap", help='f true, Seurat::HTOHeatmap will be run on the results of each caller.')
parser$add_argument("--perCellSaturation", help='An optional dataframe with the columns cellbarcode and saturation')
parser$add_argument("--majorityConsensusThreshold", help='This applies to calculating a consensus call when multiple algorithms are used')
parser$add_argument("--chemistry", help='This string is passed to EstimateMultipletRate. Should be either 10xV2 or 10xV3.')
parser$add_argument("--callerDisagreementThreshold", help='If provided, the agreement rate will be calculated between each caller and the simple majority call, ignoring discordant and no-call cells.')
parser$add_argument("--assignmentOutBff", help="Prefix name for the file containing the output of BFF assignment", type = "character", default = "bff")
parser$add_argument("--outputdir", help='Output directory')
args <- parser$parse_args()

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
metricsFile <- args$metricsFile
if(is.null(metricsFile)){
  metricsFile <- "NULL"
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
Argument <- c("seuratObject", "assay", "methods", "methodsForConsensus", "cellbarcodeWhitelist", "metricsFile", "perCellSaturation","majorityConsensusThreshold","callerDisagreementThreshold", "doTSNE","doHeatmap","chemistry")
Value <- c(seuratObj, args$assay, args$methods, methodsForConsensus, cellbarcodeWhitelist, metricsFile, perCellSaturation, majorityConsensusThreshold,callerDisagreementThreshold,args$doTSNE,args$doHeatmap,args$chemistry)

params <- data.frame(Argument, Value)
# Loading Seurat object
hashtag <-readRDS(seuratObj)

seurat_hto_counts <- hashtag[["HTO"]]@counts

cell_hash_R_res <- GenerateCellHashingCalls(barcodeMatrix = seurat_hto_counts, methods = c("bff_raw", "bff_cluster"), doTSNE = FALSE, doHeatmap = FALSE)
print("params.csv")

write.csv(params, paste0(args$outputdir, "/params.csv"))
write.csv(cell_hash_R_res, paste0(args$outputdir, "/", args$assignmentOutBff, "_assignment_bff.csv"), row.names=FALSE)

