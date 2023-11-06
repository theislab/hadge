#!/usr/bin/env Rscript

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
parser$add_argument("--fileHto", help="Path to file HTO count matrix.")
parser$add_argument("--methods", help='A vector of one or more calling methods to use.', default="combined_bff")
parser$add_argument("--methodsForConsensus", help='By default, a consensus call will be generated using all methods', default=NULL)
parser$add_argument("--cellbarcodeWhitelist", help='A vector of expected cell barcodes. This allows reporting on the total set of expected barcodes, not just those in the filtered count matrix',default=NULL)
parser$add_argument("--metricsFile", help='If provided, summary metrics will be written to this file.', default="metrics_cell_hash_r.csv")
parser$add_argument("--doTSNE", help='If true, tSNE will be run on the resulting hashing calls after each caller.', default=TRUE)
parser$add_argument("--preprocess", help='If true, the PreProcess function by CellHashR is executed', default=FALSE)
parser$add_argument("--barcodeWhitelist", help='A vector of barcode names to retain, used for preprocessing step', default=NULL)
parser$add_argument("--doHeatmap", help='f true, Seurat::HTOHeatmap will be run on the results of each caller.', default=TRUE)
parser$add_argument("--perCellSaturation", help='An optional dataframe with the columns cellbarcode and saturation',default=NULL)
parser$add_argument("--majorityConsensusThreshold", help='This applies to calculating a consensus call when multiple algorithms are used',default=NULL)
parser$add_argument("--chemistry", help='This string is passed to EstimateMultipletRate. Should be either 10xV2 or 10xV3.', default="10xV3")
parser$add_argument("--callerDisagreementThreshold", help='If provided, the agreement rate will be calculated between each caller and the simple majority call, ignoring discordant and no-call cells.',default=NULL)

parser$add_argument("--assignmentOutBff", help="Prefix name for the file containing the output of BFF assignment", type = "character", default = "bff")

parser$add_argument("--outputdir", help='Output directory')
args <- parser$parse_args()

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
print(args$fileHto)
#Transform logical parameters
do_TSNE <- as.logical(args$doTSNE)
do_Heatmap <- as.logical(args$doHeatmap)

#saving parameters in a dataframe
Argument <- c("HTO-File", "methods", "methodsForConsensus", "cellbarcodeWhitelist", "metricsFile", "perCellSaturation","majorityConsensusThreshold","callerDisagreementThreshold", "doTSNE","doHeatmap","chemistry")
Value <- c(args$fileHto, args$methods, methodsForConsensus, cellbarcodeWhitelist, args$metricsFile, perCellSaturation, majorityConsensusThreshold, callerDisagreementThreshold, args$doTSNE, args$doHeatmap,args$chemistry)
params <- data.frame(Argument, Value)

if(as.logical(args$preprocess)){
  print("Preprocessing activated")
  print(args$preprocess)
  #get barcodes
  string <- args$barcodeWhitelist
  #separate the barcodes by comma
  words <- strsplit(string, ",")[[1]]
  #Remove leading/trailing whitespace from each word
  words <- trimws(words)
  # Step 3: Create a vector from the barcodesl
  vector <- unlist(words)
  print("Preprocessing")
  #counts <- Read10X(args$fileHto) 
  counts <- ProcessCountMatrix(rawCountData = args$fileHto, barcodeBlacklist = vector)
}else{
  print("No preprocessing")
  counts <- Read10X(args$fileHto) 
}

substring_vector <- NULL
if (!is.null(args$methodsForConsensus)) { 
  consensus_methods = args$methodsForConsensus
  substring_vector <- strsplit(consensus_methods, ",")[[1]]
}

perCell_args <- args$perCellSaturation
perCell <- ifelse(perCell_args == "null" || perCell_args == "Null", NULL, perCell_args)

if(args$methodsForConsensus=="bff_raw" || args$methodsForConsensus=="bff_cluster" || args$methodsForConsensus=="bff_raw,bff_cluster" || is.null(args$methodsForConsensus)  )
  #Only Bff in its different variations is available
  if (args$methods == "bff_raw") {
    print("Executing BFF raw")
    cell_hash_R_res <- GenerateCellHashingCalls(barcodeMatrix = counts, methods = c("bff_raw"), doTSNE = do_TSNE, doHeatmap = do_Heatmap, methodsForConsensus = substring_vector,cellbarcodeWhitelist = args$cellbarcodeWhitelist, metricsFile = args$metricsFile, perCellSaturation = args$perCellSaturation, majorityConsensusThreshold = args$majorityConsensusThreshold, chemistry = args$chemistry, callerDisagreementThreshold = args$callerDisagreementThreshold )
  }else if (args$methods == "bff_cluster") {
    print("Executing BFF cluster")
    cell_hash_R_res <- GenerateCellHashingCalls(barcodeMatrix = counts, methods = c("bff_cluster"), doTSNE = do_TSNE, doHeatmap = do_Heatmap,methodsForConsensus = substring_vector,cellbarcodeWhitelist = args$cellbarcodeWhitelist,metricsFile = args$metricsFile, perCellSaturation = args$perCellSaturation, majorityConsensusThreshold = args$majorityConsensusThreshold, chemistry = args$chemistry, callerDisagreementThreshold = args$callerDisagreementThreshold)
  }else if (args$methods == "combined_bff") {
    print("Executing BFF combined")
    cell_hash_R_res <- GenerateCellHashingCalls(barcodeMatrix = counts, methods = c("bff_raw", "bff_cluster") , doTSNE = do_TSNE, doHeatmap = do_Heatmap,methodsForConsensus = substring_vector, cellbarcodeWhitelist = args$cellbarcodeWhitelist ,metricsFile = args$metricsFile, perCellSaturation = NULL, majorityConsensusThreshold = args$majorityConsensusThreshold )
    #cell_hash_R_res <- GenerateCellHashingCalls(barcodeMatrix = counts, methods = c("bff_raw", "bff_cluster") , doTSNE = do_TSNE, doHeatmap = do_Heatmap,methodsForConsensus = substring_vector,cellbarcodeWhitelist = args$cellbarcodeWhitelist, metricsFile = args$metricsFile, perCellSaturation = NULL, majorityConsensusThreshold = args$majorityConsensusThreshold, chemistry = args$chemistry, callerDisagreementThreshold = args$callerDisagreementThreshold )
  
  }else {
    print("Method not available on the pipeline")
}else{
  print("Consensus only available using BFF methods on the pipeline")
}

if(is.null(cell_hash_R_res)){
  print("No results found")
  df <- data.frame()
  write.csv(df, paste0(args$outputdir, "/", args$assignmentOutBff, "_assignment_bff.csv"), row.names=FALSE)
}else{
  write.csv(cell_hash_R_res, paste0(args$outputdir, "/", args$assignmentOutBff, "_assignment_bff.csv"), row.names=FALSE)
}
write.csv(params, paste0(args$outputdir, "/params.csv"))
