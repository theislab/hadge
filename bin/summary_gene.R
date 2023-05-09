#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("R.utils"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("stringr"))

#create parser object
parser <- ArgumentParser("Parameters for comparing parameters")
parser$add_argument("--demuxlet", help = "Folder containing output files of Demuxlet", default = NULL)
parser$add_argument("--freemuxlet", help = "Folder containing output files of Freemuxlet", default = NULL)
parser$add_argument("--vireo", help = "Folder containing output files of Vireo", default = NULL)
parser$add_argument("--souporcell", help = "Folder containing output files of Souporcell", default = NULL)
parser$add_argument("--scsplit", help = "Folder containing output files of scSplit", default = NULL)

args <- parser$parse_args()

# analysis demuxlet data ----------------------------------------------------------------------------------------------------------------------------------------------------------
demuxlet_summary <- function(demuxlet_res) {
  assign <- lapply(demuxlet_res, function(x){
    obs_res_dir <- list.files(x, pattern = ".best", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir)
    obs_res <- obs_res[,c(2,5,6)]
    obs_res$Assignment <- sapply(obs_res$BEST.GUESS,function(x){
      splitlist = strsplit(x,",")[[1]]
      if (splitlist[[1]] == splitlist[[2]]){
        splitlist[[2]]}
      else{
        "NA"}
    })
    obs_res[Assignment =="NA",]$Assignment = "doublet"
    obs_res[DROPLET.TYPE=="AMB",]$Assignment = "negative"
    colnames(obs_res)[1] <- "Barcode"
    demuxlet_assign <- obs_res[, c("Barcode", "Assignment")]
    colnames(demuxlet_assign)[2] <- basename(x)
    demuxlet_assign
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(assign, "demuxlet_assignment.csv", row.names=FALSE, quote=FALSE)
  
  classi <- assign[,-1]
  classi[classi != "negative" & classi != "doublet"] <- "singlet"
  classi$Barcode <- assign$Barcode
  classi <- classi %>% select(order(colnames(classi)))
  
  write.csv(classi, "demuxlet_classification.csv", row.names=FALSE, quote=FALSE)
  
  classi_sum <- data.frame(table(t(classi[,-1])))
  colnames(classi_sum) <- c("Classification", "Frequency")
  write.csv(classi_sum, "demuxlet_classification_summary.csv", row.names=FALSE, quote=FALSE)
  assign_sum <- data.frame(table(t(assign[,-1])))
  colnames(assign_sum) <- c("Assignment", "Frequency")
  write.csv(assign_sum, "demuxlet_assignment_summary.csv", row.names=FALSE, quote=FALSE)
  
  params <- lapply(demuxlet_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "demuxlet_params.csv", row.names=FALSE, quote=FALSE)
  
}

# analysis freemuxlet data ------------------------------------------------------------------------------------------------

freemuxlet_summary <- function(freemuxlet_res) {
  assign <- lapply(freemuxlet_res, function(x){
    obs_res_dir <- list.files(x, pattern = ".clust1.samples", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir)
    obs_res <- obs_res[,c(2,5,6)]
    obs_res$Assignment <- sapply(obs_res$BEST.GUESS,function(x){
      splitlist = strsplit(x,",")[[1]]
      if (splitlist[[1]] == splitlist[[2]]){
        splitlist[[2]]}
      else{
        "NA"}
    })
    obs_res[Assignment =="NA",]$Assignment = "doublet"
    obs_res[DROPLET.TYPE=="AMB",]$Assignment = "negative"
    colnames(obs_res)[1] <- "Barcode"
    freemuxlet_assign <- obs_res[, c("Barcode", "Assignment")]
    colnames(freemuxlet_assign)[2] <- basename(x)
    freemuxlet_assign
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(assign, "freemuxlet_assignment.csv", row.names=FALSE, quote=FALSE)
  
  classi <- assign[,-1]
  classi[classi != "negative" & classi != "doublet"] <- "singlet"
  classi$Barcode <- assign$Barcode
  classi <- classi %>% select(order(colnames(classi)))
  
  write.csv(classi, "freemuxlet_classification.csv", row.names=FALSE, quote=FALSE)
  
  classi_sum <- data.frame(table(t(classi[,-1])))
  colnames(classi_sum) <- c("Classification", "Frequency")
  write.csv(classi_sum, "freemuxlet_classification_summary.csv", row.names=FALSE, quote=FALSE)
  assign_sum <- data.frame(table(t(assign[,-1])))
  colnames(assign_sum) <- c("Assignment", "Frequency")
  write.csv(assign_sum, "freemuxlet_assignment_summary.csv", row.names=FALSE, quote=FALSE)
  
  params <- lapply(freemuxlet_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "freemuxlet_params.csv", row.names=FALSE, quote=FALSE)
  
}

# analysis vireo data ----------------------------------------------------------------------------------------------------------------------------------------------------------

vireo_summary <- function(vireo_res) {
  assign <- lapply(vireo_res, function(x){
    obs_res_dir <- list.files(x, pattern = "donor_ids.tsv", full.names = TRUE, recursive = TRUE)[1]
    obs_res <- fread(obs_res_dir)
    obs_res <- obs_res[,c(1,2)]
    obs_res[obs_res =="unassigned"] = "negative"
    colnames(obs_res) <- c("Barcode", basename(x))
    obs_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(assign, "vireo_assignment.csv", row.names=FALSE, quote=FALSE)
  
  classi <- assign[,-1]
  classi[classi != "negative" & classi != "doublet"] <- "singlet"
  classi$Barcode <- assign$Barcode
  classi <- classi %>% select(order(colnames(classi)))
  write.csv(classi, "vireo_classification.csv", row.names=FALSE, quote=FALSE)
  
  classi_sum <- data.frame(table(t(classi[,-1])))
  colnames(classi_sum) <- c("Classification", "Frequency")
  write.csv(classi_sum, "vireo_classification_summary.csv", row.names=FALSE, quote=FALSE)
  assign_sum <- data.frame(table(t(assign[,-1])))
  colnames(assign_sum) <- c("Assignment", "Frequency")
  write.csv(assign_sum, "vireo_assignment_summary.csv", row.names=FALSE, quote=FALSE)
  
  params <- lapply(vireo_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "vireo_params.csv", row.names=FALSE, quote=FALSE)
  
}

# analysis souporcell data ----------------------------------------------------------------------------------------------------------------------------------------------------------
souporcell_summary <- function(souporcell_res) {
  assign <- lapply(souporcell_res, function(x){
    obs_res_dir <- list.files(x, pattern = "clusters.tsv", full.names = TRUE, recursive = TRUE)[1]
    obs_res <- fread(obs_res_dir)
    obs_res <- obs_res[,1:3]
    obs_res[obs_res$status == "doublet",]$assignment = "doublet"
    obs_res[obs_res$status == "unassigned",]$assignment = "negative"
    obs_res <- obs_res[, -2]
    colnames(obs_res) <- c("Barcode", basename(x))
    obs_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(assign, "souporcell_assignment.csv", row.names=FALSE, quote=FALSE)
  
  classi <- assign[,-1]
  classi[classi != "negative" & classi != "doublet"] <- "singlet"
  classi$Barcode <- assign$Barcode
  classi <- classi %>% select(order(colnames(classi)))
  write.csv(classi, "souporcell_classification.csv", row.names=FALSE, quote=FALSE)

  
  classi_sum <- data.frame(table(t(classi[,-1])))
  colnames(classi_sum) <- c("Classification", "Frequency")
  write.csv(classi_sum, "souporcell_classification_summary.csv", row.names=FALSE, quote=FALSE)
  assign_sum <- data.frame(table(t(assign[,-1])))
  colnames(assign_sum) <- c("Assignment", "Frequency")
  write.csv(assign_sum, "souporcell_assignment_summary.csv", row.names=FALSE, quote=FALSE)
  
  params <- lapply(souporcell_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "souporcell_params.csv", row.names=FALSE, quote=FALSE)
  
  
}

#analysis scSplit data ----------------------------------------------------------------------------------------------------------------------------------------------------------

scsplit_summary <- function(scsplit_res) {
  assign <- lapply(scsplit_res, function(x){
    obs_res_dir <- list.files(x, pattern = "scSplit_result.csv", full.names = TRUE, recursive = TRUE)[1]
    obs_res <- fread(obs_res_dir)
    obs_res <- transform(obs_res, new=do.call(rbind, strsplit(Cluster, '-', fixed=TRUE)), stringsAsFactors=F)
    obs_res <- obs_res[,-2]
    colnames(obs_res) = c("Barcode", "Classification", "Assignment")
    obs_res[obs_res$Classification == "DBL"]$Assignment = "doublet"
    obs_res <- obs_res[, -2]
    colnames(obs_res)[2] <- basename(x)
    obs_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(assign, "scsplit_assignment.csv", row.names=FALSE, quote=FALSE)
  
  classi <- assign[,-1]
  classi[classi != "negative" & classi != "doublet"] <- "singlet"
  classi$Barcode <- assign$Barcode
  classi <- classi %>% select(order(colnames(classi)))
  write.csv(classi, "scsplit_classification.csv", row.names=FALSE, quote=FALSE)
  
  classi_sum <- data.frame(table(t(classi[,-1])))
  colnames(classi_sum) <- c("Classification", "Frequency")
  write.csv(classi_sum, "scsplit_classification_summary.csv", row.names=FALSE, quote=FALSE)
  assign_sum <- data.frame(table(t(assign[,-1])))
  colnames(assign_sum) <- c("Assignment", "Frequency")
  write.csv(assign_sum, "scsplit_assignment_summary.csv", row.names=FALSE, quote=FALSE)
  
  params <- lapply(scsplit_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "scsplit_params.csv", row.names=FALSE, quote=FALSE)
  
}
 
if (!is.null(args$demuxlet)){
  demuxlet_res <- str_split(args$demuxlet, pattern=':')[[1]]
  demuxlet_summary(demuxlet_res)
  print("Demuxlet result found")
}
if (!is.null(args$freemuxlet)){
  freemuxlet_res <- str_split(args$freemuxlet, pattern=':')[[1]]
  freemuxlet_summary(freemuxlet_res)
  print("Freemuxlet result found")
}
if (!is.null(args$vireo)){
  vireo_res <- str_split(args$vireo, pattern=':')[[1]]
  vireo_summary(vireo_res)
  print("Vireo result found")
}
if (!is.null(args$scsplit)){
  scsplit_res <- str_split(args$scsplit, pattern=':')[[1]]
  scsplit_summary(scsplit_res)
  print("scSplit result found")
}
if (!is.null(args$souporcell)){
  souporcell_res <- str_split(args$souporcell, pattern=':')[[1]]
  souporcell_summary(souporcell_res)
  print("Souporcell result found")
}

assignment <- list.files(".", pattern = "_assignment.csv", full.names = TRUE)
assignment_all <- lapply(assignment, function(x){
  assign <- fread(x, header = TRUE)
  assign
}) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
write.csv(assignment_all, "genetic_assignment_all.csv", row.names=FALSE, quote=FALSE)

classification <- list.files(".", pattern = "_classification.csv", full.names = TRUE)
classification_all <- lapply(classification, function(x){
  classi <- fread(x, header = TRUE)
  classi
}) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
write.csv(classification_all, "genetic_classification_all.csv", row.names=FALSE, quote=FALSE)
