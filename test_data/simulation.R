library(data.table)
library(tidyr)
library(Matrix)
library(Seurat)
library(DropletUtils)
# use demuxlet tutorial dataset to simulate hashing and rna expression matrix
test_data_barcode <- fread("jurkat_293t_demuxlet.best")
test_data_barcode <- test_data_barcode[,c("BARCODE", "BEST")]
write(x=test_data_barcode$BARCODE, "barcodes.tsv")
test_data_barcode <- separate(test_data_barcode, col = BEST, into = c("classification", "identity"), sep = "\\-")
test_data_barcode_singlet <- test_data_barcode[test_data_barcode$classification == "SNG",]$BARCODE
test_data_barcode_singlet_jurkat <- test_data_barcode[test_data_barcode$classification == "SNG" &
                                                   test_data_barcode$identity=="jurkat",]$BARCODE
test_data_barcode_singlet_293T <- test_data_barcode[test_data_barcode$classification == "SNG" &
                                                        test_data_barcode$identity=="293T_RTG",]$BARCODE
test_data_barcode_doublet <- test_data_barcode[test_data_barcode$classification == "DBL",]$BARCODE

set.seed(10000)
# Simulation is adapted from the tutorial of DropletUtils
# Simulating empty droplets:
nbarcodes <- 500
nhto <- 2
ngene <- 300
count_matrix_hto <- matrix(rpois(nbarcodes*nhto, 20), nrow=nhto)
count_matrix_rna <- matrix(rpois(nbarcodes*ngene, 20), nrow=ngene)

# Simulating cells:
ncells <- 400
ndoub <- ncells/10

cell_index <- sample(nhto, ncells, replace=TRUE)
cell_index[(ndoub+1):ncells] <- sort(cell_index[(ndoub+1):ncells])
count_matrix_hto[cbind(cell_index, seq_len(ncells))] <- 1000

# Simulating doublets:
doublet_index <- (cell_index[1:ndoub] + 1) %% nrow(count_matrix_hto)
doublet_index[doublet_index==0] <- nrow(count_matrix_hto)
count_matrix_hto[cbind(doublet_index, seq_len(ndoub))] <- 500
doublet <- test_data_barcode_doublet[1:ndoub]

# Assign donor identity to cells
cell_index <- cell_index[(ndoub+1):ncells]
singlet_jurkat <- test_data_barcode_singlet_jurkat[1:length(cell_index[cell_index==1])]
singlet_293T <- test_data_barcode_singlet_293T[1:length(cell_index[cell_index==2])]
non_negative <- c(doublet, singlet_jurkat, singlet_293T)
negative <- setdiff(test_data_barcode$BARCODE, non_negative) 

length(singlet_jurkat)
length(singlet_293T)
length(doublet)
length(negative)

# Merge all droplets
singlet_jurkat_dt <- data.frame(barcode=singlet_jurkat, identity="jurkat")
singlet_293T_dt <- data.frame(barcode=singlet_293T, identity="293T_RTG")
doublet_dt <- data.frame(barcode=doublet, identity="doublet")
negative_dt <- data.frame(barcode=negative, identity="negative")
barcodes <- rbind(doublet_dt, singlet_jurkat_dt, singlet_293T_dt, negative_dt)

# Save sparse matrix in 10x format
rownames(count_matrix_hto) <- paste0("hto_", 1:2)
colnames(count_matrix_hto) <- barcodes$barcode
sparse.count_matrix_hto <- Matrix(count_matrix_hto, sparse = T)
write10xCounts("hto",sparse.count_matrix_hto)

rownames(count_matrix_rna) <- paste0("gene_", 1:ngene)
colnames(count_matrix_rna) <- barcodes$barcode
sparse.count_matrix_rna <- Matrix(count_matrix_rna, sparse = T)
write10xCounts("rna",sparse.count_matrix_rna)

# Check output
rna_data <- Read10X("rna")
hto_data <- Read10X("hto")