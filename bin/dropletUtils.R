#!/usr/bin/env Rscript
library(Seurat)
library(DropletUtils)
library(argparse)

# for emptyDrops
parser <-
  ArgumentParser("Parameters for Empty Drops cell identification")
parser$add_argument("--raw_hto_matrix_dir",
  help = "A numeric/integer matrix-like object containing UMI counts prior to any filtering or cell calling. Rows correspond to HTOs and columns correspond to cell barcodes. "
)
parser$add_argument("--lower",
  default = 100,
  type = "double",
  help = "A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets."
)
parser$add_argument("--niters",
  default = 10000,
  type = "integer",
  help = "An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations."
)
parser$add_argument("--testAmbient",
  action = "store_true",
  help = "A logical scalar indicating whether results should be returned for barcodes with totals less than or equal to lower."
)
parser$add_argument("--ignore",
  default = NULL,
  type = "double",
  help = "A numeric scalar specifying the lower bound on the total UMI count, at or below which barcodes will be ignored."
)
parser$add_argument("--alpha",
  default = NULL,
  type = "double",
  help = "A numeric scalar specifying the scaling parameter for the Dirichlet-multinomial sampling scheme."
)
parser$add_argument("--round",
  action = "store_true",
  help = "Logical scalar indicating whether to check for non-integer values in m and, if present, round them for ambient profile estimation."
)
parser$add_argument("--byRank",
  default = NULL,
  type = "integer",
  help = "An integer scalar parametrizing an alternative method for identifying assumed empty droplets. If set, this is used to redefine lower and any specified value for lower is ignored."
)
parser$add_argument("--isCellFDR",
  default = 0.01,
  type = "double",
  help = "Threshold to filter the cells."
)
parser$add_argument("--objectOutEmptyDrops",
  default = "emptyDroplets",
  help = "Prefix name for the emptyDrops RDS file"
)
parser$add_argument("--assignmentOutEmptyDrops",
  default = "emptyDroplets",
  help = "prefex name for emptyDrops assignment CSV file"
)

# for hashedDrops
parser$add_argument("--ambient",
  action = "store_true",
  help = "Whether to use the relative abundance of each HTO in the ambient solution from emtpyDrops, set TRUE only when test_ambient is TRUE."
)
parser$add_argument("--minProp",
  default = 0.05,
  type = "double",
  help = "Numeric scalar to be used to infer the ambient profile when ambient=NULL,"
)
parser$add_argument("--pseudoCount",
  default = 5,
  type = "double",
  help = "A numeric scalar specifying the minimum pseudo-count when computing logfold changes."
)
parser$add_argument("--constantAmbient",
  action = "store_true",
  help = "Logical scalar indicating whether a constant level of ambient contamination should be used to estimate LogFC2 for all cells."
)
parser$add_argument("--doubletNmads",
  default = 3,
  type = "double",
  help = "A numeric scalar specifying the number of median absolute deviations (MADs) to use to identify doublets."
)
parser$add_argument("--doubletMin",
  default = 2,
  type = "double",
  help = "A numeric scalar specifying the minimum threshold on the log-fold change to use to identify doublets."
)
parser$add_argument("--doubletMixture",
  action = "store_true",
  help = "Logical scalar indicating whether to use a 2-component mixture model to identify doublets."
)
parser$add_argument(
  "--confidentNmads",
  default = 3,
  type = "double",
  help = "A numeric scalar specifying the number of MADs to use to identify confidently assigned singlets."
)
parser$add_argument("--confidenMin",
  default = 2,
  type = "double",
  help = "A numeric scalar specifying the minimum threshold on the log-fold change to use to identify singlets."
)
parser$add_argument("--combinations",
  default = NULL,
  type = "integer",
  help = "An integer matrix specifying valid combinations of HTOs. Each row corresponds to a single sample and specifies the indices of rows in x corresponding to the HTOs used to label that sample."
)
parser$add_argument("--objectOutHashedDrops",
  default = "hashedDrops",
  help = "Prefix name for the hashedDrops RDS file"
)
parser$add_argument("--assignmentOutHashedDrops",
  default = "hashedDrops",
  help = "prefex name for hashedDrops assignment CSV file"
)
parser$add_argument("--outputdir", help = "Output directory")
parser$add_argument("--gene_col",
  help = "Specify which column of genes.tsv or features.tsv to use for gene names; default is 2",
  type = "integer",
  default = 2
)
parser$add_argument("--runEmptyDrops", action="store_false",
                     help = "Executes emptyDrops function only when desired, recomended only for raw data")


args <- parser$parse_args()

hto <-
  Read10X(
    data.dir = args$raw_hto_matrix_dir,
    gene.column = args$gene_col
  )
combinations_transformed <-
  ifelse(tolower(args$combinations) == "null", NULL, args$combinations)

if (args$runEmptyDrops == TRUE) {
  rna <-
    Read10X(
      data.dir = args$raw_rna_matrix_dir,
      gene.column = args$gene_col
    )
  print("------------------- executing emptyDrops ---------------------------------")
  ignore_transformed <-
    ifelse(tolower(args$ignore) == "null", NULL, args$ignore)
  emptyDrops_out <- emptyDrops(
    rna,
    lower = args$lower,
    niters = args$niters,
    test.ambient = args$testAmbient,
    ignore = NULL,
    alpha = args$alpha,
    round = args$round,
    by.rank = args$byRank
  )
  write.csv(
    emptyDrops_out,
    paste0(args$outputdir, "/", args$assignmentOutEmptyDrops, ".csv")
  )
  saveRDS(emptyDrops_out,
    file = paste0(args$outputdir, "/", args$objectOutEmptyDrops, ".rds")
  )

  print("------------------- filtering empty droplets ----------------------------")
  is.cell <- emptyDrops_out$FDR <= args$isCellFDR
  colors <- ifelse(is.cell, "red", "black")
  png(paste0(args$outputdir, "/", "plot_emptyDrops.png"))
  plot(
    emptyDrops_out$Total, -emptyDrops_out$LogProb,
    col = colors,
    xlab = "Total UMI count",
    ylab = "-Log Probability"
  )
  dev.off()
  

  if (args$ambient == TRUE) {
    hashedDrops_out <-
      hashedDrops(
        hto[, which(is.cell)],
        min.prop = args$minProp,
        ambient = metadata(emptyDrops_out)$ambient,
        pseudo.count = args$pseudoCount,
        constant.ambient = args$constantAmbient,
        doublet.nmads = args$doubletNmads,
        doublet.min = args$doubletMin,
        doublet.mixture = args$doubletMixture,
        confident.nmads = args$confidentNmads,
        confident.min = args$confidenMin,
        combinations = combinations_transformed
      )
  } else {
    hashedDrops_out <-
      hashedDrops(
        hto[, which(is.cell)],
        min.prop = args$minProp,
        pseudo.count = args$pseudoCount,
        constant.ambient = args$constantAmbient,
        doublet.nmads = args$doubletNmads,
        doublet.min = args$doubletMin,
        doublet.mixture = args$doubletMixture,
        confident.nmads = args$confidentNmads,
        confident.min = args$confidenMin,
        combinations = combinations_transformed
      )
  }
} else {
  hashedDrops_out <-
    hashedDrops(
      hto,
      min.prop = args$minProp,
      pseudo.count = args$pseudoCount,
      constant.ambient = args$constantAmbient,
      doublet.nmads = args$doubletNmads,
      doublet.min = args$doubletMin,
      confident.nmads = args$confidentNmads,
      confident.min = args$confidenMin
    )
}

ignore <- args$ignore
if (is.null(ignore)) {
  ignore <- "NULL"
}

alpha <- args$alpha
if (is.null(alpha)) {
  alpha <- "NULL"
}

byRank <- args$byRank
if (is.null(byRank)) {
  byRank <- "NULL"
}

minProp <- args$minProp
if (is.null(minProp)) {
  minProp <- "NULL"
}

combinations <- args$combinations
if (is.null(combinations)) {
  combinations <- "NULL"
}

Argument <- c(
  "raw_hto_matrix_dir",
  "lower",
  "niters",
  "testAmbient",
  "ignore",
  "alpha",
  "round",
  "byRank",
  "isCellFDR",
  "ambient",
  "minProp",
  "pseudoCount",
  "constantAmbient",
  "doubletNmads",
  "doubletMin",
  "doubletMixture",
  "confidentNmads",
  "confidenMin",
  "combinations"
)
Value <- c(
  args$raw_hto_matrix_dir,
  args$lower,
  args$niters,
  args$testAmbient,
  ignore,
  alpha,
  args$round,
  byRank,
  args$isCellFDR,
  args$ambient,
  minProp,
  args$pseudoCount,
  args$constantAmbient,
  args$doubletNmads,
  args$doubletMin,
  args$doubletMixture,
  args$confidentNmads,
  args$confidenMin,
  combinations
)

params <- data.frame(Argument, Value)

print("-------- Following Files are saved in folder hashedDrops_out ------------")
print(paste0(args$objectOutHashedDrops, ".rds"))
print(paste0(args$assignmentOutHashedDrops, "_res.csv"))
print("params.csv")
write.csv(params, paste0(args$outputdir, "/params.csv"))
write.csv(
  hashedDrops_out,
  paste0(
    args$outputdir,
    "/",
    args$assignmentOutHashedDrops,
    "_res.csv"
  )
)
saveRDS(hashedDrops_out,
  file = paste0(args$outputdir, "/", args$objectOutHashedDrops, ".rds")
)

colors <- ifelse(hashedDrops_out$Confident,
  "black",
  ifelse(hashedDrops_out$Doublet, "red", "grey")
)
png(paste0(args$outputdir, "/", "plot_hashed.png"))
if (sum(is.na(hashedDrops_out$LogFC2)) != length(hashedDrops_out$LogFC2)) {
  print(paste0(args$objectOutHashedDrops, "_LogFC.png"))
  plot(
    hashedDrops_out$LogFC,
    hashedDrops_out$LogFC2,
    col = colors,
    xlab = "Log-fold change between best and second HTO",
    ylab = "Log-fold change between second HTO and ambient"
  )
  dev.off()
}
