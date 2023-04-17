#!/usr/bin/env Rscript
library(data.table)
library(stringr)
library(argparse)
library(dplyr)
library("R.utils")

parser <- ArgumentParser("Parameters for comparing parameters")
parser$add_argument("--demuxem", help = "Folder containing output files of demuxem", default = NULL)
parser$add_argument("--htodemux", help = "Folder containing output files of htodemux", default = NULL)
parser$add_argument("--multiseq", help = "Folder containing output files of multiseq", default = NULL)
parser$add_argument("--hashsolo", help = "Folder containing output files of hashsolo", default = NULL)
parser$add_argument("--solo", help = "Folder containing output files of solo", default = NULL)
parser$add_argument("--hashedDrops", help = "Folder containing output files of hashedDrops", default = NULL)
parser$add_argument("--select", help = "Select the singlets detected by the specified number of tools or fewer.", default = 1)
args <- parser$parse_args()

demuxem_summary <- function(demuxem_res) {
  assign <- lapply(demuxem_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_obs.csv", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir, header = TRUE)
    colnames(obs_res)[1] <- "Barcode"
    obs_res[obs_res$demux_type == "unknown",]$assignment <- "negative"
    demuxem_assign <- obs_res[, c("Barcode", "assignment")]
    colnames(demuxem_assign)[2] <- basename(x)
    demuxem_assign
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  assign$Barcode <- paste0(assign$Barcode, '-1')
  write.csv(assign, "demuxem_assignment.csv", row.names=FALSE)
  
  classi <- lapply(demuxem_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_obs.csv", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir, header = TRUE)
    colnames(obs_res)[1] <- "Barcode"
    demuxem_classi <- obs_res[, c("Barcode", "demux_type")]
    colnames(demuxem_classi)[2] <- basename(x)
    demuxem_classi
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  classi$Barcode <- paste0(classi$Barcode, '-1')
  classi[classi == "unknown"] <- "negative"
  write.csv(classi, "demuxem_classification.csv", row.names=FALSE, quote=FALSE)
  
  params <- lapply(demuxem_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "demuxem_params.csv", row.names=FALSE, quote=FALSE)
  
}

hashsolo_summary <- function(hashsolo_res){
  assign <- lapply(hashsolo_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_res.csv", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir, header = TRUE)
    hashsolo_assign <- obs_res[, c("V1", "Classification")]
    colnames(hashsolo_assign) <- c("Barcode", basename(x))
    hashsolo_assign
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  assign[assign == "Doublet"] <- "doublet"
  assign[assign == "Negative"] <- "negative"
  write.csv(assign, "hashsolo_assignment.csv", row.names=FALSE, quote=FALSE)
  
  classi <- lapply(hashsolo_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_res.csv", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir, header = TRUE)
    hashsolo_classi <- obs_res[, c("V1", "most_likely_hypothesis")]
    hashsolo_classi$most_likely_hypothesis[hashsolo_classi$most_likely_hypothesis == 0]  <- "negative"
    hashsolo_classi$most_likely_hypothesis[hashsolo_classi$most_likely_hypothesis == 1]  <- "singlet"
    hashsolo_classi$most_likely_hypothesis[hashsolo_classi$most_likely_hypothesis == 2]  <- "doublet"
    colnames(hashsolo_classi) <- c("Barcode", basename(x))
    hashsolo_classi
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(classi, "hashsolo_classification.csv", row.names=FALSE, quote=FALSE)
  
  params <- lapply(hashsolo_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "hashsolo_params.csv", row.names=FALSE, quote=FALSE)
}

hasheddrops_summary <- function(hasheddrops_res){
  classi <- lapply(hasheddrops_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_res.csv", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir, header = TRUE)
    obs_res$Classification = ifelse(obs_res$Confident, "singlet", ifelse(obs_res$Doublet, "doublet", "negative"))
    colnames(obs_res)[1] <- 'Barcode'
    hasheddrops_classi <- obs_res[,c('Barcode', 'Classification')]
    colnames(hasheddrops_classi) <- c("Barcode", basename(x))
    hasheddrops_classi
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(classi, "hasheddrops_classification.csv", row.names=FALSE, quote=FALSE)
  
  params <- lapply(hasheddrops_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    params_res <- params_res[,2:3]
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "hasheddrops_params.csv", row.names=FALSE, quote=FALSE)
}

multiseq_summary <- function(multiseq_res){
  assign <- lapply(multiseq_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_res.csv", full.names = TRUE)[1]
    multiseq_assign <- fread(obs_res_dir, header = TRUE)
    colnames(multiseq_assign) = c("Barcode", basename(x))
    multiseq_assign
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  assign[assign == "Doublet"] <- "doublet"
  assign[assign == "Negative"] <- "negative"
  write.csv(assign, "multiseq_assignment.csv", row.names=FALSE, quote=FALSE)
  
  classi <- assign[,-1]
  classi[classi != "doublet" & classi != "negative"] <- "singlet"
  classi$Barcode <- assign$Barcode
  classi <- classi %>% select(order(colnames(classi)))
  write.csv(classi, "multiseq_classification.csv", row.names=FALSE, quote=FALSE)
  
  params <- lapply(multiseq_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    params_res <- params_res[,2:3]
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "multiseq_params.csv", row.names=FALSE, quote=FALSE)
}

htodemux_summary <- function(htodemux_res){
  assign <- lapply(htodemux_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_assignment_htodemux.csv", full.names = TRUE)[1]
    htodemux_assign <- fread(obs_res_dir, header = TRUE)
    colnames(htodemux_assign) = c("Barcode", basename(x))
    htodemux_assign
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  assign[assign == "Doublet"] <- "doublet"
  assign[assign == "Negative"] <- "negative"
  write.csv(assign, "htodemux_assignment.csv", row.names=FALSE, quote=FALSE)
  
  classi <- lapply(htodemux_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_classification_htodemux.csv", full.names = TRUE)[1]
    htodemux_classi <- fread(obs_res_dir, header = TRUE)
    colnames(htodemux_classi) = c("Barcode", basename(x))
    htodemux_classi
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  
  classi_barcode <- classi$Barcode
  classi <- classi[, -1]
  classi[classi == "Singlet"] <- "singlet"
  classi[classi == "Doublet"] <- "doublet"
  classi[classi == "Negative"] <- "negative"
  classi$Barcode <- classi_barcode
  classi <- classi %>% select(order(colnames(classi)))
  write.csv(classi, "htodemux_classification.csv", row.names=FALSE, quote=FALSE)
  
  params <- lapply(htodemux_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    params_res <- params_res[,2:3]
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "htodemux_params.csv", row.names=FALSE, quote=FALSE)
}

solo_summary <- function(solo_res){
  classi <- lapply(solo_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_res.csv", full.names = TRUE)[1]
    solo_classi <- fread(obs_res_dir, header = TRUE)
    colnames(solo_classi) = c("Barcode", basename(x))
    solo_classi$Barcode <- gsub("-0", "", solo_classi$Barcode)
    solo_classi
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(classi, "solo_classification.csv", row.names=FALSE, quote=FALSE)
  
  params <- lapply(solo_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "solo_params.csv", row.names=FALSE, quote=FALSE)
}


if (!is.null(args$hashedDrops)){
  hashedDrops_res <- substring(args$hashedDrops, 1, nchar(args$hashedDrops)-1)
  hashedDrops_res <- str_split(hashedDrops_res, pattern=':')[[1]]
  hasheddrops_summary(hashedDrops_res)
  print("hashedDrops result found")
}
if (!is.null(args$demuxem)){
  demuxem_res <- substring(args$demuxem, 1, nchar(args$demuxem)-1)
  demuxem_res <- str_split(demuxem_res, pattern=':')[[1]]
  demuxem_summary(demuxem_res)
  print("DemuxEM result found")
}
if (!is.null(args$hashsolo)){
  hashsolo_res <- substring(args$hashsolo, 1, nchar(args$hashsolo)-1)
  hashsolo_res <- str_split(hashsolo_res, pattern=':')[[1]]
  hashsolo_summary(hashsolo_res)
  print("HashSolo result found")
}
if (!is.null(args$multiseq)){
  multiseq_res <- substring(args$multiseq, 1, nchar(args$multiseq)-1)
  multiseq_res <- str_split(multiseq_res, pattern=':')[[1]]
  multiseq_summary(multiseq_res)
  print("MultiSeqDemux result found")
}
if (!is.null(args$htodemux)){
  htodemux_res <- substring(args$htodemux, 1, nchar(args$htodemux)-1)
  htodemux_res <- str_split(htodemux_res, pattern=':')[[1]]
  htodemux_summary(htodemux_res)
  print("HTODemux result found")
}
if (!is.null(args$solo)){
  solo_res <- substring(args$solo, 1, nchar(args$solo)-1)
  solo_res <- str_split(solo_res, pattern=':')[[1]]
  solo_summary(solo_res)
  print("solo result found")
}

assignment <- list.files(".", pattern = "_assignment.csv", full.names = TRUE)
assignment_all <- lapply(assignment, function(x){
  assign <- fread(x, header = TRUE)
  assign
}) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
write.csv(assignment_all, "hash_assignment_all.csv", row.names=FALSE)

classification <- list.files(".", pattern = "_classification.csv", full.names = TRUE)
classification_all <- lapply(classification, function(x){
  classi <- fread(x, header = TRUE)
  classi
}) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)

write.csv(classification_all, "hash_classification_all.csv", row.names=FALSE, quote=FALSE)
#classification_all$count_dou <- apply(classification_all[,-1], 1, function(x) sum(x == 'doublet', na.rm=TRUE))
#classification_all$count_neg <- apply(classification_all[,-1], 1, function(x) sum(x == 'negative', na.rm=TRUE))
classification_all$count_sin <- apply(classification_all[,-1], 1, function(x) sum(x == 'singlet', na.rm=TRUE))


# if(args$select != 0 ){
write.table(classification_all[classification_all$count_sin<=args$select,]$Barcode, file='selected_barcodes.tsv', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
# }else{
    # write.table(classification_all[,1], file='selected_barcodes.tsv', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
# }
