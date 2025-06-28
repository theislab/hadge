#!/usr/bin/env Rscript

################################################
################################################
## USE PARAMETERS FROM NEXTFLOW               ##
################################################
################################################

# cast parameters from nextflow
rna_matrix <- '$rna_matrix'
hto_matrix <- '$hto_matrix'
sel_method <- '$sel_method'
ndelim <- '$ndelim'
n_features <- as.numeric('$n_features')
assay <- '$assay'
margin <- as.numeric('$margin')
norm_method <- '$norm_method'
gene_col <- as.numeric('$gene_col')
prefix <- '$prefix'

# check if the files exist
if (! file.exists(hto_matrix)){
    stop(paste0(hto_matrix, ' is not a valid file'))
}

if (! file.exists(rna_matrix)){
    stop(paste0(rna_matrix, ' is not a valid file'))
}

################################################
################################################
## Load libraries                             ##
################################################
################################################

library(Seurat)

################################################
################################################
## Main Process                               ##
################################################
################################################

# Read 10X data
umi <- Read10X(data.dir = rna_matrix, gene.column = gene_col)
counts <- Read10X(data.dir = hto_matrix, gene.column = gene_col)

# Select cell barcodes detected by both RNA and HTO
joint.bcs <- intersect(colnames(umi), colnames(counts))
# Subset RNA and HTO counts by joint cell barcodes
umi <- umi[, joint.bcs]
counts <- counts[, joint.bcs]

# Setup Seurat object
hashtag <- CreateSeuratObject(counts = umi, names.delim = ndelim)
# Normalize RNA data with log normalization
hashtag <- NormalizeData(hashtag)
# Find and scale variable features
hashtag <- FindVariableFeatures(hashtag, selection.method = sel_method)
hashtag <- ScaleData(hashtag, features = VariableFeatures(hashtag))
# Add HTO data as a new assay independent from RNA
hashtag[[assay]] <- CreateAssayObject(counts = counts)
# Normalize HTO data
hashtag <- NormalizeData(hashtag, assay = assay, normalization.method = norm_method, margin = margin)

################################################
################################################
## SAVING RESULTS                             ##
################################################
################################################

# Save preprocessed Seurat object
saveRDS(hashtag, file = paste0(prefix, "_preprocessed.rds"))

# Save parameters
Argument <- c(
  "hto_matrix",
  "rna_matrix", 
  "sel_method",
  "ndelim",
  "n_features",
  "assay",
  "margin",
  "norm_method",
  "gene_col"
)

Value <- c(
  hto_matrix,
  rna_matrix,
  sel_method,
  ndelim,
  n_features,
  assay,
  margin,
  norm_method,
  gene_col
)

params <- data.frame(Argument, Value)
write.csv(params, paste0(prefix, "_params_preprocessing.csv")) 

################################################
################################################
## SAVE VERSIONS                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
seurat.version <- as.character(packageVersion('Seurat'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    seurat:', seurat.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################