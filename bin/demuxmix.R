#!/usr/bin/env Rscript
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org")
is_demuxmix<-require("demuxmix")
if (!require("demuxmix", quietly = TRUE))
    BiocManager::install("demuxmix")

library(Seurat)
library(demuxmix)
library(argparse)
library(data.table)


# Create a parser
parser <- ArgumentParser("Parameters for Demuxmix")
parser$add_argument("--seuratObject", help = "Seurat object. Assumes that the hash tag oligo (HTO) data has been added and normalized.")
parser$add_argument("--pAcpt", help='Acceptance probability that must be reached in order to assign a droplet to a hashtag. ')
parser$add_argument("--assay", help='Assay name')
parser$add_argument("--rna_available", help='TRUE if RNA assay is available',default = FALSE)
parser$add_argument("--correctTails", help='If TRUE, droplets meeting the threshold defined by alpha (beta) are classified as "negative" ("positive") even if the mixture model suggests a different classification',default = TRUE)
parser$add_argument("--model", help='A character specifying the type of mixture model to be used. Either "naive", "regpos", "reg" or "auto".', default = 'naive')
parser$add_argument("--alpha_demuxmix", help='Threshold defining the left tail of the mixture distribution where droplets should not be classified as "positive".',type = "double",default=0.9)
parser$add_argument("--beta_demuxmix", help='Threshold for defining the right tail of the mixture distribution where droplets should not be classified as "negative".',type = "double",default=0.9)
parser$add_argument("--tol_demuxmix", help='Convergence criterion for the EM algorithm used to fit the mixture models.', type="double")
parser$add_argument("--maxIter_demuxmix", help='Maximum number of iterations for the EM algorithm',type = "integer",default=100)
parser$add_argument("--k_hto", help='Factor to define outliers in the HTO counts.',type = "double",default=1.5)
parser$add_argument("--k_rna", help='Factor to define outliers in the numbers of detected genes.',type = "double",default=1.5)
parser$add_argument("--assignmentOutDemuxmix", help="Prefix name for the file containing the output of Demuxmix assignment", type = "character", default = "demuxmix")
parser$add_argument("--outputdir", help='Output directory')
args <- parser$parse_args()

#check if seurat object is correct
if (!endsWith(args$seuratObject, ".rds")){
    seuratObj <- list.files(args$seuratObject, pattern = "\\.rds$", full.names = TRUE)[1]
    print(seuratObj)

}else{
    seuratObj <- args$seuratObject
}


Argument <- c("seuratObject", "assay", "model", "alpha_demuxmix", "beta_demuxmix", "tol_demuxmix", "maxIter_demuxmix", "k_hto","k_rna","correct_tails")
Value <- c(seuratObj, args$assay, args$model, args$alpha_demuxmix, args$beta_demuxmix, args$tol_demuxmix, args$maxIter_demuxmix, args$k_hto, args$k_rna,args$correctTails)

params <- data.frame(Argument, Value)


# Loading Seurat object
hashtag <-readRDS(seuratObj)


### Demuxmix
hto_counts <- as.matrix(GetAssayData(hashtag[[args$assay]], slot = "counts"))
hashtags_object <- GetAssayData(hashtag[[args$assay]], slot = "counts")
hashtag_list<-hashtags_object@Dimnames


#Demultiplexing process
if(args$model != 'naive' && args$rna_available == 'TRUE'){
    rna_counts <- hashtag$nCount_RNA
    demuxmix_demul <- demuxmix(hto_counts,rna= rna_counts, model = args$model, alpha= args$alpha_demuxmix,beta= args$beta_demuxmix,maxIter=args$maxIter_demuxmix,k.hto=args$k_hto, correctTails = as.logical(args$correctTails))

}else{
    print("Executing naive mode for Demuxmix")
    demuxmix_demul <- demuxmix(hto_counts, model = args$model, alpha= args$alpha_demuxmix,beta= args$beta_demuxmix,maxIter=args$maxIter_demuxmix,k.hto=args$k_hto)
}

demuxmix_classify <- dmmClassify(demuxmix_demul)
sumary_demuxmix <- summary(demuxmix_demul)
demuxmix_classify$Barcode <- hashtag_list[[2]]

res_dt <- as.data.table(demuxmix_classify)
res_dt[, Classification := Type]
res_dt[Classification == "multiplet", Classification := "doublet"]
res_dt[Classification == "uncertain", Classification := "negative"]
 
res_dt$HTO <- gsub(',', '_', res_dt$HTO)

write.csv(params, paste0(args$outputdir, "/params.csv"))
write.csv(res_dt, paste0(args$outputdir, "/", args$assignmentOutDemuxmix, "_assignment_demuxmix.csv"), row.names=FALSE)
write.csv(sumary_demuxmix, paste0(args$outputdir, "/", args$assignmentOutDemuxmix, "_summary_results.csv"), row.names=FALSE)



