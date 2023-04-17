#!/usr/bin/env Rscript
#Libraries
library(Seurat)
library(argparse)
library(ggplot2)

# Create a parser
parser <- ArgumentParser("Parameters for MultiSeqDemux")
parser$add_argument("--seuratObjectPath", help = "Seurat object or the folder containing the seurat object")
parser$add_argument("--quantile", help = "The quantile to use for classification", type = "double", default = 0.7)
parser$add_argument("--autoThresh", help = "Whether to perform automated threshold finding to define the best quantile", action = "store_true")
parser$add_argument("--maxiter", help = "Maximum number of iterations if autoThresh = TRUE ", default = 5, type="integer")
parser$add_argument("--qrangeFrom", help = "A range of possible quantile values to try if autoThresh is TRUE", type = "double", default = 0.1)
parser$add_argument("--qrangeTo", help = "A range of possible quantile values to try if autoThresh is TRUE", type = "double", default = 0.9)
parser$add_argument("--qrangeBy", help = "A range of possible quantile values to try if autoThresh is TRUE", type = "double", default = 0.05)
parser$add_argument("--verbose", help = "Prints the output", action = "store_true")
parser$add_argument("--assay", help = "Name of the multiplexing assay", default="HTO")
parser$add_argument("--assignmentOutMulti",help = "Prefix name for the file containing the output of MULTI-Seq Demux object", default = "mutliseq")
parser$add_argument("--objectOutMulti", help = "Prefix name for the object containing the output of MULTI-Seq Demux object", default = "multiseq")
parser$add_argument("--outputdir", help='Output directory')

args <- parser$parse_args()
if (!endsWith(args$seuratObjectPath, ".rds")){
    seuratObj <- list.files(args$seuratObjectPath, pattern = "\\.rds$", full.names = TRUE)[1]
}else{
    seuratObj <- args$seuratObjectPath
}

Argument <- c("seuratObjectPath", "quantile", "autoThresh", "maxiter", "qrangeFrom", "qrangeTo", "qrangeBy", "verbose", "assay")
Value <- c(seuratObj, args$quantile, args$autoThresh, args$maxiter, args$qrangeFrom, args$qrangeTo, args$qrangeBy, args$verbose, args$assay)

params <- data.frame(Argument, Value)

hashtag <-readRDS(seuratObj)
if (args$autoThresh == TRUE) {
    hashtag <- MULTIseqDemux(hashtag, assay = args$assay,  quantile = args$quantile, autoThresh = TRUE, maxiter = args$maxiter, qrange=seq(from = args$qrangeFrom, to = args$qrangeTo, by = args$qrangeBy), verbose = args$verbose)
}else{
    hashtag <- MULTIseqDemux(hashtag, assay = args$assay, quantile = args$quantile, verbose = args$verbose)
    
}

table(hashtag$MULTI_ID)
print("----------------------------------------------------------------------------")
table(hashtag$MULTI_classification)

hashtag
print("-----------------------------------------------")

dim(x = hashtag)
# head(x = rownames(x = hashtag))
# head(x = colnames(x = hashtag))
print("-----------------------------------------------")
names(x = hashtag)

print("-----------------------------------------------")
hashtag[[args$assay]]

print("-----------------------------------------------")
colnames(x = hashtag[[]])

#Save Results
print("------------------- Following Files are saved ----------------------------")
print(paste0(args$assignmentOutMulti, "_res.csv"))
print(paste0(args$objectOutMulti,".rds"))
print("params.csv")
write.csv(hashtag$MULTI_ID, paste0(args$outputdir, "/", args$assignmentOutMulti, "_res.csv"))
write.csv(params, paste0(args$outputdir, "/params.csv"))
saveRDS(hashtag, paste0(args$outputdir, "/", args$objectOutMulti, ".rds"))
