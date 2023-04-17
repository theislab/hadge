#!/usr/bin/env Rscript
#Libraries
library(Seurat)
library(argparse)
library(ggplot2)

# Create a parser
parser <- ArgumentParser("Parameters for HTODemux")
parser$add_argument("--seuratObject", help = "Seurat object. Assumes that the hash tag oligo (HTO) data has been added and normalized.")
parser$add_argument("--quantile", help = "The quantile of inferred 'negative' distribution for each hashtag, over which the cell is considered positive. ", type = "double", default = 0.99)
parser$add_argument("--kfunc", help = "Clustering function for initial hashtag grouping. Default is clara for fast k-medoids clustering on large applications, also support kmeans for kmeans clustering.", type = "character", default = "clara")
parser$add_argument("--nstarts", help = "nstarts value for k-means clustering.", type = "integer", default = 100)
parser$add_argument("--nsamples", help = "Number of samples to be drawn from the dataset used for clustering, for kfunc = clara", type = "integer", default = 100)
parser$add_argument("--seed", help = "Sets the random seed.", type = "integer", default = 42)
parser$add_argument("--init", help = "Initial number of clusters for hashtags.", default = NULL, type = "integer")
parser$add_argument("--assay", help = "Assay name", default = "HTO")
parser$add_argument("--objectOutHTOdemux", help = "Prefix name for the object containing the output of HTODemux object", type = "character", default = "htodemux")
parser$add_argument("--assignmentOutHTOdemux", help="Prefeix name for the file containing the output of HTODemux assignment", type = "character", default = "htodemux")
parser$add_argument("--outputdir", help='Output directory')


args <- parser$parse_args()
if (!endsWith(args$seuratObject, ".rds")){
    seuratObj <- list.files(args$seuratObject, pattern = "\\.rds$", full.names = TRUE)[1]
}else{
    seuratObj <- args$seuratObject
}

init <- args$init
if(is.null(init)){
    init <- "NULL"
}
Argument <- c("seuratObject", "quantile", "kfunc", "nstarts", "nsamples", "seed", "init", "assay")
Value <- c(seuratObj, args$quantile, args$kfunc, args$nstarts, args$nsamples, args$seed, init, args$assay)

params <- data.frame(Argument, Value)
# Loading Seurat object
hashtag <-readRDS(seuratObj)

# Demultiplex cells based on HTO enrichment
if(args$kfunc == "clara"){
    hashtag <- HTODemux(hashtag, assay = args$assay, positive.quantile = args$quantile, init = args$init, nsamples = args$nsamples, kfunc = "clara", seed = args$seed)
}else{
    hashtag <- HTODemux(hashtag, assay = args$assay, positive.quantile = args$quantile, init = args$init, nstarts = args$nstarts, kfunc = "kmeans", seed = args$seed)
}


# Global classification results
table(hashtag[[paste0(args$assay,"_classification.global")]])

table(hashtag[[paste0(args$assay,"_classification")]])

print("-----------------------------------------------")

hashtag[[args$assay]]

print("-----------------------------------------------")

hashtag

# Saving results

donors <- rownames(hashtag[[args$assay]])
assignment <- hashtag[[paste0(args$assay,"_classification")]]
assignment[[paste0(args$assay,"_classification")]][!assignment[[paste0(args$assay,"_classification")]] %in% c(donors, 'Negative')] <- "Doublet"


print("------------------- Following Files are saved ----------------------------")
print(paste0(args$assignmentOutHTOdemux, "_assignment_htodemux.csv"))
print(paste0(args$assignmentOutHTOdemux, "_classification_htodemux.csv"))
print(paste0(args$objectOutHTOdemux,".rds"))
print("params.csv")
write.csv(params, paste0(args$outputdir, "/params.csv"))
write.csv(assignment, paste0(args$outputdir, "/", args$assignmentOutHTOdemux, "_assignment_htodemux.csv"))
#write.csv(hashtag[[paste0(args$assay,"_classification")]], paste0(args$outputdir, "/", args$assignmentOutHTOdemux, "_assignment_htodemux.csv"))
write.csv(hashtag[[paste0(args$assay,"_classification.global")]], paste0(args$outputdir, "/", args$assignmentOutHTOdemux, "_classification_htodemux.csv"))
saveRDS(hashtag, file=paste0(args$outputdir, "/", args$objectOutHTOdemux,".rds"))


