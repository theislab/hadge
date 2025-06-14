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

# Write to h5 file
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
