#!/usr/bin/env Rscript
library(Seurat)
library(DropletUtils)
library(argparse)

# for emptyDrops
parser <- ArgumentParser("Parameters for Empty Drops cell identification")
parser$add_argument("--raw_hto_matrix_dir", help = "A numeric/integer matrix-like object containing UMI counts prior to any filtering or cell calling. Rows correspond to HTOs and columns correspond to cell barcodes. ")
parser$add_argument("--lower",help="A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.", default = 100, type = "double")
parser$add_argument("--niters", help = "An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations.", default = 10000, type="integer")
parser$add_argument("--testAmbient",help="A logical scalar indicating whether results should be returned for barcodes with totals less than or equal to lower.", action='store_true')
parser$add_argument("--ignore",help="A numeric scalar specifying the lower bound on the total UMI count, at or below which barcodes will be ignored.", default = NULL, type = "double")
parser$add_argument("--alpha",help="A numeric scalar specifying the scaling parameter for the Dirichlet-multinomial sampling scheme.", default = NULL, type = "double")
parser$add_argument("--round",help="Logical scalar indicating whether to check for non-integer values in m and, if present, round them for ambient profile estimation.", action='store_true')
parser$add_argument("--byRank",help="An integer scalar parametrizing an alternative method for identifying assumed empty droplets. If set, this is used to redefine lower and any specified value for lower is ignored.", default = NULL, type = "integer")
parser$add_argument("--isCellFDR",help="Threshold to filter the cells.", default = 0.01, type = "double")

parser$add_argument("--objectOutEmptyDrops",help="Prefix name for the emptyDrops RDS file", default = "emptyDroplets")
parser$add_argument("--assignmentOutEmptyDrops",help="prefex name for emptyDrops assignment CSV file", default = "emptyDroplets")

#for hashedDrops
parser$add_argument("--ambient",help="A numeric vector of length equal to nrow(matrix), specifying the relative abundance of each HTO in the ambient solution.", default = NULL, type = "double")
parser$add_argument("--minProp",help="Numeric scalar to be used to infer the ambient profile when ambient=NULL,", default = 0.05, type = "double")
parser$add_argument("--pseudoCount",help="A numeric scalar specifying the minimum pseudo-count when computing logfold changes.", default = 5, type = "double")
parser$add_argument("--constantAmbient",help="Logical scalar indicating whether a constant level of ambient contamination should be used to estimate LogFC2 for all cells.", action='store_true')
parser$add_argument("--doubletNmads",help="A numeric scalar specifying the number of median absolute deviations (MADs) to use to identify doublets.", default = 3, type = "double")
parser$add_argument("--doubletMin",help="A numeric scalar specifying the minimum threshold on the log-fold change to use to identify doublets.", default = 2, type = "double")
parser$add_argument("--doubletMixture",help="Logical scalar indicating whether to use a 2-component mixture model to identify doublets.", action='store_true')
parser$add_argument("--confidentNmads",help="A numeric scalar specifying the number of MADs to use to identify confidently assigned singlets.", default = 3, type = "double")
parser$add_argument("--confidenMin",help="A numeric scalar specifying the minimum threshold on the log-fold change to use to identify singlets.", default = 2, type = "double")
parser$add_argument("--combinations",help="An integer matrix specifying valid combinations of HTOs. Each row corresponds to a single sample and specifies the indices of rows in x corresponding to the HTOs used to label that sample.", default = NULL, type = "integer")
parser$add_argument("--objectOutHashedDrops",help="Prefix name for the hashedDrops RDS file", default = "hashedDrops")
parser$add_argument("--assignmentOutHashedDrops",help="prefex name for hashedDrops assignment CSV file", default = "hashedDrops")
parser$add_argument("--outputdir", help='Output directory')

args <- parser$parse_args()

hto <- Read10X(data.dir = args$raw_hto_matrix_dir)

emptyDrops_out <- emptyDrops(hto, lower = args$lower, niters = args$niters,  
                             test.ambient = args$testAmbient, ignore = args$ignore,
                             alpha = args$alpha, round = args$round,
                             by.rank = args$byRank)

print("------------------- emptyDrops finished ---------------------------------")


print("-------- Following Files are saved in folder hashedDrops_out ------------")
print(paste0(args$objectOutEmptyDrops, ".rds"))
print(paste0(args$assignmentOutEmptyDrops, ".csv"))
write.csv(emptyDrops_out, paste0(args$outputdir, "/", args$assignmentOutEmptyDrops, ".csv"))
saveRDS(emptyDrops_out, file=paste0(args$outputdir, "/", args$objectOutEmptyDrops, ".rds"))

print("------------------- filtering empty droplets ----------------------------")
is.cell <- emptyDrops_out$FDR <= args$isCellFDR
colors <- ifelse(is.cell, "red", "black")
png(paste0(args$outputdir, "/", "plot_emptyDrops.png"))
plot(emptyDrops_out$Total, -emptyDrops_out$LogProb, col=colors, xlab="Total UMI count", ylab="-Log Probability")
dev.off()

hashedDrops_out <- hashedDrops(hto[,which(is.cell)], min.prop = args$minProp, ambient = args$ambient, pseudo.count = args$pseudoCount, constant.ambient = args$constantAmbient, doublet.nmads = args$doubletNmads, doublet.min = args$doubletMin, doublet.mixture = args$doubletMixture, confident.nmads = args$confidentNmads, confident.min = args$confidenMin, combinations = args$combinations)

print("------------------- hashedDrops finished ---------------------------------")

ignore <- args$ignore
if(is.null(ignore)){
    ignore <- "NULL"
}

alpha <- args$alpha
if(is.null(alpha)){
    alpha <- "NULL"
}

byRank <- args$byRank
if(is.null(byRank)){
    byRank <- "NULL"
}

ambient <- args$ambient
if(is.null(ambient)){
    ambient <- "NULL"
}

minProp <- args$minProp
if(is.null(minProp)){
    minProp <- "NULL"
}

combinations <- args$combinations
if(is.null(combinations)){
   combinations <- "NULL"
}
Argument <- c("raw_hto_matrix_dir", "lower", "niters", "testAmbient", "ignore", "alpha", "round", "byRank", "isCellFDR", "ambient", "minProp", "pseudoCount", "constantAmbient", "doubletNmads", "doubletMin", "doubletMixture", "confidentNmads", "confidenMin", "combinations")
Value <- c(args$raw_hto_matrix_dir, args$lower, args$niters, args$testAmbient, ignore, alpha, args$round, byRank, args$isCellFDR, ambient, minProp, args$pseudoCount, args$constantAmbient, args$doubletNmads, args$doubletMin, args$doubletMixture, args$confidentNmads, args$confidenMin, combinations)

params <- data.frame(Argument, Value)

print("-------- Following Files are saved in folder hashedDrops_out ------------")
print(paste0(args$objectOutHashedDrops, ".rds"))
print(paste0(args$assignmentOutHashedDrops, "_res.csv"))
print(paste0(args$objectOutHashedDrops, "_LogFC.png"))
print("params.csv")
write.csv(params, paste0(args$outputdir, "/params.csv"))
write.csv(hashedDrops_out, paste0(args$outputdir, "/", args$assignmentOutHashedDrops, "_res.csv"))
saveRDS(hashedDrops_out, file=paste0(args$outputdir, "/", args$objectOutHashedDrops, ".rds"))

colors <- ifelse(hashedDrops_out$Confident, "black", ifelse(hashedDrops_out$Doublet, "red", "grey"))
png(paste0(args$outputdir, "/", "plot_hashed.png"))
plot(hashedDrops_out$LogFC, hashedDrops_out$LogFC2, col=colors,
                xlab="Log-fold change between best and second HTO",
                ylab="Log-fold change between second HTO and ambient")
dev.off()
