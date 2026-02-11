#! /usr/bin/env Rscript

library(DropletUtils)

mtx_dir <- "${input_mtx_dir}"


sce <- read10xCounts(mtx_dir) # Read to SingleCellExperiment object

print(sce)

# Convert to matrix
count_matrix <- counts(sce)
rownames(count_matrix) <- rownames(sce)
colnames(count_matrix) <- colData(sce)\$Barcode

if ("${write_csv}" == "true") {
    write.csv(as.matrix(count_matrix), file = "${prefix}.csv", row.names = TRUE)
}

# TODO demuxem: remove if demuxEM issue is solved (https://github.com/theislab/hadge/issues/81)
# Write to h5 file
# write10xCounts(
#   path        = "${prefix}.h5",
#   x           = counts(sce),
#   barcodes    = colData(sce)\$Barcode,
#   gene.id     = rownames(sce),
#   gene.symbol = if (!is.null(rowData(sce)\$Symbol)) rowData(sce)\$Symbol else rownames(sce),
#   gene.type   = if (!is.null(rowData(sce)\$Type))   rowData(sce)\$Type   else rep("Gene Expression", nrow(sce)),
#   type        = "HDF5",
#   version     = "3",           # <-- ensures /matrix layout instead of /unknown
#   overwrite   = TRUE
# )
write10xCounts("${prefix}.h5", count_matrix, type = "HDF5")

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
dropletutils.version <- as.character(packageVersion('DropletUtils'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-dropletutils:', dropletutils.version)
    ),
'versions.yml')

############################################
