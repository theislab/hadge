#!/usr/bin/env Rscript
library(data.table)
library(stringr)
library(argparse)
library(dplyr)

parser <- ArgumentParser("Parameters for comparing parameters")
parser.add_argument("--demuxem", help="Folder containing output files of demuxem", default=None)
parser.add_argument("--htodemux", help="Folder containing output files of htodemux", default=None)
parser.add_argument("--multiseq", help="Folder containing output files of multiseq", default=None)
parser.add_argument("--hashsolo", help="Folder containing output files of hashsolo", default=None)
parser.add_argument("--solo", help="Folder containing output files of solo", default=None)
parser.add_argument("--hashedDrops", help="Folder containing output files of hashedDrops", default=None)
parser$add_argument("--gmm_demux", help = "Folder containing output files of GMM Demux", default = None)
parser$add_argument("--demuxmix", help = "Folder containing output files of Demuxmix", default = None)
parser$add_argument("--bff", help = "Folder containing output files of BFF", default = None)
parser.add_argument("--generate_anndata", help="Generate anndata", action='store_true')
parser.add_argument("--read_mtx", help="10x-Genomics-formatted mtx directory", default=None)
args = parser.parse_args()


demuxem_summary <- function(demuxem_res) {
  assign <- lapply(demuxem_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_obs.csv", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir, header = TRUE)
    colnames(obs_res)[1] <- "Barcode"
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
  assign <- lapply(hasheddrops_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_res.csv", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir, header = TRUE)
    obs_res$Classification <- ifelse(obs_res$Confident, "singlet", ifelse(obs_res$Doublet, "doublet", "negative"))
    obs_res$Best <- ifelse(!obs_res$Classification %in% c("doublet", "negative"), obs_res$Best, obs_res$Classification)
    colnames(obs_res)[1] <- 'Barcode'
    hasheddrops_assign <- obs_res[, c('Barcode', 'Best')]
    colnames(hasheddrops_assign) <- c("Barcode", basename(x))
    hasheddrops_assign
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(assign, "hasheddrops_assignment.csv", row.names=FALSE, quote=FALSE)

  classi <- lapply(hasheddrops_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_res.csv", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir, header = TRUE)
    obs_res$Classification <- ifelse(obs_res$Confident, "singlet", ifelse(obs_res$Doublet, "doublet", "negative"))
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

gmm_summary<-function(gmmDemux_res){
  classi <- lapply(gmmDemux_res, function(x){
    obs_res_dir <- list.files(x, pattern = "GMM_full.csv", full.names = TRUE)[1]
    gmm_classi <- fread(obs_res_dir, header = TRUE)
    #read classification file generated by GMM
    classification_config <- list.files(x, pattern = "GMM_full.config", full.names = TRUE)[1]
    classif_file <- fread(classification_config,header = TRUE)
    #both files as data tables
    gmm_dt <- as.data.table(gmm_classi)
    classification_dt <- as.data.table(classif_file)
    #create three types of classification according to type of assignment
    classification_dt$V2 <- grepl("-", classification_dt$V2)  # FALSE
    classification_dt[, Classification := ifelse( grepl("-", classification_dt$V2), "doublet", "singlet")]
    classification_dt[1,]$Classification <- "negative"
    gmm_dt[,Classification:="negative"]
    #check which case of classification we have in the results according with classification file
    result <- sapply(gmm_dt$Cluster_id, function(x) any(x == classification_dt$V1))
    gmm_dt$Classification <- ifelse(result, classification_dt$Classification, NA)
    gmm_dt<-gmm_dt[,c("V1", "Classification")]
    #Change column name for result and add gmm_name
    colnames(gmm_dt) = c("Barcode", basename(x))
    #trim -1 from barcodes
    gmm_dt$Barcode <- gsub("-1", "", gmm_dt$Barcode)
    gmm_dt
  })%>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(classi, "gmm_classification.csv", row.names=FALSE, quote=FALSE)
  assign <- lapply(gmmDemux_res, function(x){
    obs_res_dir <- list.files(x, pattern = "GMM_full.csv", full.names = TRUE)[1]
    gmm_classi <- fread(obs_res_dir, header = TRUE)
    #read classification file generated by GMM
    classification_config <- list.files(x, pattern = "GMM_full.config", full.names = TRUE)[1]
    classif_file <- fread(classification_config,header = TRUE)
    #both files as data tables
    gmm_dt <- as.data.table(gmm_classi)
    classification_dt <- as.data.table(classif_file)
    result <- sapply(gmm_dt$Cluster_id, function(x) any(x == classification_dt$V1))
    gmm_dt$Assignment <- ifelse(result, classification_dt$V2, NA)
    gmm_dt<-gmm_dt[,c("V1", "Assignment")]
    #Change column name for result and add gmm_name
    colnames(gmm_dt) = c("Barcode", basename(x))
    gmm_dt$Barcode <- gsub("-1", "", gmm_dt$Barcode)
    gmm_dt
  })%>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(assign, "GMM_assignment.csv", row.names=FALSE, quote=FALSE)
  params <- lapply(gmmDemux_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "GMM_params.csv", row.names=FALSE, quote=FALSE)
}
demuxmix_summary<-function(demuxmix_res){
  classi <- lapply(demuxmix_res, function(x){
    obs_res_dir <- list.files(x, pattern = "assignment_demuxmix.csv", full.names = TRUE)[1]
    demuxmix_classi <- fread(obs_res_dir, header = TRUE)
    dt_classi <- as.data.table(demuxmix_classi)
    dt_classi<-dt_classi[,c("Barcode","Classification")]
    colnames(dt_classi) = c("Barcode", basename(x))
    dt_classi
  })%>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(classi, "demuxmix_classification.csv", row.names=FALSE, quote=FALSE)
  assign <- lapply(demuxmix_res, function(x){
    obs_res_dir <- list.files(x, pattern = "assignment_demuxmix.csv", full.names = TRUE)[1]
    demuxmix_asign <- fread(obs_res_dir, header = TRUE)
    dt_assign <- as.data.table(demuxmix_asign)
    dt_assign<-dt_assign[,c("Barcode","HTO")]
    colnames(dt_assign) = c("Barcode", basename(x))
    dt_assign
  })%>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(assign, "demuxmix_assignment.csv", row.names=FALSE, quote=FALSE)
}
bff_summary<-function(bff_res){
classi <- lapply(bff_res, function(x){
  obs_res_dir <- list.files(x, pattern = "bff.csv", full.names = TRUE)[1]
  bff_classi <- fread(obs_res_dir, header = TRUE)
  data_bff <- as.data.table(bff_classi)
  dt_classi<-data_bff[,c("cellbarcode","consensuscall.global")]
  print(dt_classi)
  dt_classi[consensuscall.global == "Singlet", consensuscall.global:="singlet"]
  dt_classi[consensuscall.global == "Doublet", consensuscall.global:="doublet"]
  dt_classi[consensuscall.global == "Negative", consensuscall.global:="negative"]
  colnames(dt_classi) = c("Barcode", basename(x))
  dt_classi
  })%>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .) 
  write.csv(classi, "bff_classification.csv", row.names=FALSE, quote=FALSE)
  assign <- lapply(bff_res, function(x){
  obs_res_dir <- list.files(x, pattern = "bff.csv", full.names = TRUE)[1]
  bff_assign <- fread(obs_res_dir, header = TRUE)
  data_bff <- as.data.table(bff_assign)
  dt_assign<-data_bff[,c("cellbarcode","consensuscall")]
  dt_assign[consensuscall=="Doublet", consensuscall:="doublet"]
  colnames(dt_assign) = c("Barcode", basename(x))
  dt_assign
  })%>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(assign, "bff_assignment.csv", row.names=FALSE, quote=FALSE)   
}


if (!is.null(args$hashedDrops)){
  hashedDrops_res <- str_split(args$hashedDrops, pattern=':')[[1]]
  hasheddrops_summary(hashedDrops_res)
  print("hashedDrops result found")
}
if (!is.null(args$demuxem)){
  demuxem_res <- str_split(args$demuxem, pattern=':')[[1]]
  demuxem_summary(demuxem_res)
  print("DemuxEM result found")
}
if (!is.null(args$hashsolo)){
  hashsolo_res <- str_split(args$hashsolo, pattern=':')[[1]]
  hashsolo_summary(hashsolo_res)
  print("HashSolo result found")
}
if (!is.null(args$multiseq)){
  multiseq_res <- str_split(args$multiseq, pattern=':')[[1]]
  multiseq_summary(multiseq_res)
  print("MultiSeqDemux result found")
}
if (!is.null(args$htodemux)){
  htodemux_res <- str_split(args$htodemux, pattern=':')[[1]]
  htodemux_summary(htodemux_res)
  print("HTODemux result found")
}
if (!is.null(args$solo)){
  solo_res <- str_split(args$solo, pattern=':')[[1]]
  solo_summary(solo_res)
  print("solo result found")
}
if (!is.null(args$gmm_demux)){
  gmmDemux_res <- substring(args$gmm_demux, 1, nchar(args$gmm_demux)-1)
  gmmDemux_res <- str_split(gmmDemux_res, pattern=':')[[1]]
  gmm_summary(gmmDemux_res)
  print("GMMDemux result found")
}
if (!is.null(args$demuxmix)){
  demuxmix_res <- substring(args$demuxmix, 1, nchar(args$demuxmix)-1)
  demuxmix_res <- str_split(demuxmix_res, pattern=':')[[1]]
  demuxmix_summary(demuxmix_res)
  print("Demuxmix result found")
}

if (!is.null(args$bff)){
  bff_res <- substring(args$bff, 1, nchar(args$bff)-1)
  bff_res <- str_split(bff_res, pattern=':')[[1]]
  bff_summary(bff_res)
  print("BFF result found")
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
