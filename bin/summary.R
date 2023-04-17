#!/usr/bin/env Rscript
library(data.table)
library(stringr)
library(argparse)
library(dplyr)
library(R.utils)

parser <- ArgumentParser("Parameters for summarizing results")
parser$add_argument("--gene_demulti", help = "Folder containing output files of genetic demultiplexing pipeline", default = NULL)
parser$add_argument("--hash_demulti", help = "Folder containing output files of hashing demultiplexing pipeline", default = NULL)
args <- parser$parse_args()
assignment_gene <- list.files(args$gene_demulti, pattern = "_assignment_all.csv", full.names = TRUE)
assignment_hash <- list.files(args$hash_demulti, pattern = "_assignment_all.csv", full.names = TRUE)
assignment <- append(assignment_gene, assignment_hash)

ifelse(!dir.exists("summary"), dir.create("summary"), FALSE)

assignment_all <- lapply(assignment, function(x){
  assign <- fread(x, header = TRUE)
  assign
}) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
assignment_all[assignment_all == "DBL"] <- "doublet"
assignment_all[assignment_all == "AMB"] <- "negative"
#write.table(assignment_all, "summary/assignment_all_genetic_and_hash.csv", row.names=F, quote=F, sep='\t')
write.table(assignment_all, "summary/assignment_all_genetic_and_hash.csv", row.names=F, quote=F)

classification_gene <- list.files(args$gene_demulti, pattern = "_classification_all.csv", full.names = TRUE)
classification_hash <- list.files(args$hash_demulti, pattern = "_classification_all.csv", full.names = TRUE)
classification <- append(classification_gene, classification_hash)
classification_all <- lapply(classification, function(x){
  classi <- fread(x, header = TRUE)
  classi
}) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
#classification_all <- classification_all[, -c("count_sin", "count_dou", "count_neg")]
classification_all[classification_all == "SNG"] <- "singlet"
classification_all[classification_all == "DBL"] <- "doublet"
classification_all[classification_all == "AMB"] <- "negative"
#write.table(classification_all, "summary/classification_all_genetic_and_hash.csv", row.names=F, quote=F, sep='\t')
write.table(classification_all, "summary/classification_all_genetic_and_hash.csv", row.names=F, quote=F)
