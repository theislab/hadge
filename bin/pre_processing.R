#!/usr/bin/env Rscript
#Libraries
library(Seurat)
library(argparse)

# Create a parser
parser <- ArgumentParser("Parameters for Seurat object: Pre-processing")

parser$add_argument("--fileUmi", help = "Path to file UMI count matrix.")
parser$add_argument("--fileHto", help = "Path to file HTO count matrix.")
parser$add_argument("--ndelim", help = "For the initial identity calss for each cell, delimiter for the cell's column name", default = "_")
parser$add_argument("--rdsObject", help = "whether inputs are rds objects", action = "store_true")
parser$add_argument("--selectMethod", help = "Selection method", default = "mean.var.plot")
parser$add_argument("--numberFeatures", help = "Number of features to be used when finding variable features", type = "integer", default = 2000)
parser$add_argument("--assay", help = "Assay name for hashing modality", default = "HTO")
parser$add_argument("--raw_data", help = "Demuxmix and BFF only require the object without pre-processing", default="false")
parser$add_argument("--rna_available", help = "Demuxmix only requires the object without pre-processing", default="true")
parser$add_argument("--normalisationMethod", help = "Normalisation method", default = "CLR")
parser$add_argument("--margin", help = "Margin for normalisation", type="integer", default = 2)

parser$add_argument( "--OutputFile", help="Prefix of output files containing the output of HTODemux hashtag", default = "preprocessed")
parser$add_argument("--outputdir", help='Output directory')
args <- parser$parse_args()



#Check if RNA data is available
if(args$rdsObject){
  if(as.logical(args$rna_available)){
    umi <- readRDS(args$fileUmi)
  }
  counts <- readRDS(args$fileHto)
  
}else{
  if(as.logical(args$rna_available)){ 
    
    umi <- Read10X(data.dir = args$fileUmi)
  }
  
  counts <- Read10X(data.dir = args$fileHto)
}


if(as.logical(args$raw_data)){
  if(as.logical(args$rna_available)){
    Argument <- c("fileUmi", "fileHto")
    Value <- c(args$fileUmi, args$fileHto)
    params <- data.frame(Argument, Value)
    #demuxmix requires a Seurat object created from raw data, without any type of pre-processing (Normalisation, ScaleData, etc)
    # Subset RNA and HTO counts by joint cell barcodes
    #joint.bcs <- intersect(colnames(umi), colnames(counts))
    #umi <- umi[, joint.bcs]
    #counts <- counts[, joint.bcs]
    # Setup Seurat object
    hashtag <- CreateSeuratObject(counts = umi, names.delim = args$ndelim)
    hashtag[[args$assay]] <- CreateAssayObject(counts = counts)
  }else{
    print("create a hto seurat object")
    hashtag <- CreateSeuratObject(counts = counts, names.delim = args$ndelim, assay = args$assay)
    Argument <- c( "fileHto")
    Value <- c(args$fileHto)
    params <- data.frame(Argument, Value)
  }
}else {
  
#If the parameter demuxmix is not passed, then the regular pre-processing is performed
args <- parser$parse_args()
Argument <- c("fileUmi", "fileHto", "ndelim", "rdsObject", "selectMethod", "numberFeatures", "assay", "normalisationMethod", "margin")
Value <- c(args$fileUmi, args$fileHto, args$ndelim, args$rdsObject, args$selectMethod, args$numberFeatures, args$assay, args$normalisationMethod, args$margin)
params <- data.frame(Argument, Value)
# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
joint.bcs <- intersect(colnames(umi), colnames(counts))
# Subset RNA and HTO counts by joint cell barcodes
umi <- umi[, joint.bcs]
counts <- counts[, joint.bcs]
# Setup Seurat object
hashtag <- CreateSeuratObject(counts = umi, names.delim = args$ndelim)
# Normalize RNA data with log normalization
hashtag <- NormalizeData(hashtag)
# Find and scale variable features
hashtag <- FindVariableFeatures(hashtag, selection.method = args$selectMethod)
hashtag <- ScaleData(hashtag, features = VariableFeatures(hashtag))
# Add HTO data as a new assay independent from RNA
hashtag[[args$assay]] <- CreateAssayObject(counts = counts)
# Normalize HTO data
hashtag <- NormalizeData(hashtag, assay = args$assay, normalization.method = args$normalisationMethod, margin = args$margin)
  
}

#Save Results
saveRDS(hashtag, file = paste0(args$outputdir, "/", args$OutputFile, ".rds"))
write.csv(params, paste0(args$outputdir, "/params.csv"))





