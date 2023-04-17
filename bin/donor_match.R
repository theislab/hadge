#!/usr/bin/env Rscript
options(future.globals.maxSize = 4000 * 1024^5)
library(pheatmap)
library(data.table)
library(magrittr)
library(ComplexUpset)
library(ggplot2)
library(tidyverse)
library(vcfR)
library(argparse)

parser <- ArgumentParser("Parameters for donor matching based on Pearson coorelation. ")
parser$add_argument("--result_csv", help = "The path to the csv file or the directory containing the demultiplexing assignment. ")
parser$add_argument("--barcode", help = "The path to the tsv file containing the white list of barcodes. ", default = NULL)
parser$add_argument("--method1", help="The Name of the first method in the assignment file to map donors. ")
parser$add_argument("--method2", help="The Name of the second method in the assignment file. We will use its donor names as reference.")
parser$add_argument("--findVariants", help='How to find representative variants for each donor. ', default='default')
parser$add_argument("--cell_genotype", help="The path to the VCF file containing the genotype of the cells.")
parser$add_argument("--variant_count", help='The Minimal count of a variant for filtering', type="integer", default = 0)
parser$add_argument("--variant_pct", help="The Minimal percentage of a variant for filtering", type="double", default = 0.5)
parser$add_argument("--vireo_parent_dir", help="The parent folder for searching the output of the best vireo trial, e.g. 
                      donor genotype and discriminatory variants. ")
parser$add_argument("--outputdir", help='Output directory. ')

args <- parser$parse_args()

result_csv <- fread("/Users/xichenwu/hagen/work/f5/d8c5c0b4a96cccf7761320017aa893/gene_summary/genetic_assignment_all.csv", stringsAsFactors=F, sep='\t')
#result_csv <- "/Users/xichenwu/hagen/result/parallel/summary"
method1 <- "vireo_1"
method2 <- "vireo_2"
#cell_genotype <- "/Users/xichenwu/genotype_match/cellSNP.cells.vcf.gz"
#cell_genotype <- "/Users/xichenwu/openpipeline/resources_test/vireo_test_data/cells.cellSNP.vcf.gz"
#method1 <- "vireo_1"
#method2 <- "multiseq_1"
#barcode_whitelist <- fread("/Volumes/wxicu/gx12/rna/filtered_feature_bc_matrix/barcodes.tsv", header = F, stringsAsFactors = F)$V1

convert2binary <- function(result_csv, method_name){
  method_assign <- result_csv %>% select("Barcode", method_name)
  donor_id <- setdiff(unique(method_assign[[method_name]]), c(NA, 'negative', 'doublet'))
  method_assign <- method_assign[method_assign[[method_name]] %in% donor_id,]
  method_assign_binary = data.frame(model.matrix(~ method_assign[[method_name]]-1, data=method_assign))
  names(method_assign_binary) = sort(donor_id)
  rownames(method_assign_binary) <- method_assign$Barcode
  #method_assign_binary <- method_assign_binary %>% select(Barcode, everything())
  return(method_assign_binary)
}

result_csv <- NULL
if(file.exists(args$result_csv) && !dir.exists(args$result_csv)){
  result_csv <- fread(args$result_csv, stringsAsFactors=F, sep='\t')
}
if(dir.exists(args$result_csv)){
  result_csv <- list.files(args$result_csv, pattern = "assignment_all", full.names = T)
  result_csv <- fread(result_csv[1], stringsAsFactors=F, sep='\t')
}
if (!is.null(args$barcode)){
  barcode_whitelist <- fread(args$barcode, header = F, stringsAsFactors = F)$V1
  result_csv <- result_csv[result_csv$Barcode %in% barcode_whitelist, ]
}
method1_all <- colnames(result_csv)[startsWith(colnames(result_csv), args$method1)]
method2_all <- colnames(result_csv)[startsWith(colnames(result_csv), args$method2)]
best_score <- 0
best_method1 <- NULL
best_method2 <- NULL

if (is.null(method1_all) | is.null(method2_all)){
  stop("The specified method is not found in the CSV file!")
}

for (method1 in method1_all){
  for (method2 in method2_all){
    outputdir <- file.path(args$outputdir, paste0(method1, "_vs_", method2))
    ifelse(!dir.exists(outputdir), dir.create(outputdir), FALSE)
    
    method1_res <- convert2binary(result_csv, method1)
    method2_res <- convert2binary(result_csv, method2)
    
    intersect_barcode <- intersect(rownames(method1_res), rownames(method2_res))
    method1_res <- method1_res[rownames(method1_res) %in% intersect_barcode,]
    method2_res<- method2_res[rownames(method2_res) %in% intersect_barcode,]
    
    correlation_res <- apply(method1_res, 2, function(x){
      apply(method2_res, 2, function(y){
        return(cor.test(x, y)[["estimate"]][["cor"]])
      })}
    )
    write.csv(correlation_res, file.path(outputdir, "correlation_res.csv"))
    
    match_score <- 0
    geno_match <- as.data.frame(matrix(nrow = ncol(correlation_res), ncol = 3))
    colnames(geno_match) <- c("Method1_ID","Method2_ID","Correlation")
    geno_match$Method1_ID <- colnames(correlation_res)
    for (id in geno_match$Method1_ID){
      if (max(correlation_res[,id]) == max(correlation_res[which.max(correlation_res[,id]),])){
        geno_match[which(geno_match$Method1_ID == id),2:3] <- c(rownames(correlation_res)[which.max(correlation_res[,id])], max(correlation_res[,id]))
        match_score <- match_score + max(correlation_res[,id])
      } 
      else {
        geno_match[which(geno_match$Cluster1_ID == id)] <- c("unassigned", NA)
      }
    }
    
    write.table(geno_match[,1:2], file.path(outputdir, "donor_match.csv"), row.names = F, col.names = F, sep = " ", quote = F)
    if (match_score > best_score){
      write.table(geno_match[,1:2], file.path(args$outputdir, "donor_match.csv"), row.names = F, col.names = F, sep = " ", quote = F)
    }
    best_method1 <- ifelse(match_score > best_score, method1, best_method1)
    best_method2 <- ifelse(match_score > best_score, method2, best_method2)
    best_score <- ifelse(match_score > best_score, match_score, best_score)
    
    newCols <- colorRampPalette(grDevices::rainbow(nrow(geno_match)))
    annoCol <- newCols(nrow(geno_match))
    names(annoCol) <- colnames(correlation_res)
    annoCol <- list(category = annoCol)
    correlation_res <- correlation_res[order(as.numeric(row.names(correlation_res))), ]
    pheatmap(correlation_res, treeheight_row=F, treeheight_col=F, display_numbers=T, angle_col="45",
             number_color = "white", fontsize = 12, cluster_rows = F, cluster_cols = F, 
             filename=file.path(outputdir, "concordance_heatmap.png"))
    
    result_merge <- select(result_csv, "Barcode", method1, method2)
    result_merge_new <- result_merge
    
    for (i in 1:nrow(geno_match)){
      result_merge_new[[method1]] <- replace(result_merge_new[[method1]], result_merge[[method1]] == geno_match$Method1_ID[i], 
                                                  geno_match$Method2_ID[i])
    }
    
    write.csv(result_merge_new, file.path(outputdir, "all_assignment_after_match.csv"), row.names = F)
    result_merge_new <- result_merge_new[result_merge_new$Barcode %in% intersect_barcode,]
    write.csv(result_merge_new, file.path(outputdir, "intersect_assignment_after_match.csv"), row.names = F)
    
  }
}

print(best_score)
print(best_method1)
print(best_method2)

if (args$findVariants == 'True' | args$findVariants == 'default'){
  if (startsWith(best_method1, 'vireo')){
    write.table(best_method1, file.path(args$outputdir, "best_method_vireo.txt"), sep = "\t", row.names = F, col.names = F, quote = F)
    
  } else{
    stop("Vireo is not the best method1 for donor matching!")
  }
  outputdir <- file.path(args$outputdir, paste0(best_method1, "_vs_", best_method2))
  result_merge_new <- fread(file.path(outputdir, "intersect_assignment_after_match.csv"), header = T)
  result_merge_new$match <- result_merge_new[[best_method1]] == result_merge_new[[best_method2]]
  matched <- result_merge_new[result_merge_new$match,]
  unmatched <- result_csv[!result_csv$Barcode %in% matched$Barcode, ]$Barcode
  cell_genotype_vcf <- read.vcfR(args$cell_genotype)
  cell_genotype_vcf_gt <- extract.gt(cell_genotype_vcf, element = "GT", as.numeric=TRUE)
  # cells <- colnames(cell_genotype_vcf_gt)
  # matched <- matched[matched$Barcode %in% cells,]
  # unmatched <- unmatched[unmatched %in% cells]
  donors <- sort(unique(matched[[best_method1]]))
  
  representative_variant_list <- vector(mode='list', length=length(donors))
  representative_variant_list <- setNames(representative_variant_list, donors)
  
  for (donorid in donors){
    matched_barcode <- matched[matched[[best_method1]]==donorid]$Barcode
    matched_gt_list <- cell_genotype_vcf_gt[, matched_barcode]
    matched_gt_list <- matched_gt_list[rowSums(is.na(matched_gt_list)) != ncol(matched_gt_list), ]
    matched_gt <- as.data.frame(matrix(nrow = nrow(matched_gt_list)))
    matched_gt$ref <- rowSums(matched_gt_list==0,na.rm = T)
    matched_gt$alt <- rowSums(matched_gt_list!=0,na.rm = T)
    matched_gt$V1 <- rownames(matched_gt_list)
    matched_gt$count <-  matched_gt$ref + matched_gt$alt
    matched_gt$pct <- matched_gt$alt / (matched_gt$ref+matched_gt$alt)
    matched_gt$dominant <- ifelse(matched_gt$pct>0.5, 1, 0)
    matched_gt <- matched_gt[(matched_gt$pct>=args$variant_pct |matched_gt$pct<=(1-args$variant_pct)),]
    matched_gt <- matched_gt[matched_gt$count>=args$variant_count,]
    
    unmatched_gt_list <- cell_genotype_vcf_gt[, unmatched]
    unmatched_gt_list <- unmatched_gt_list[rownames(unmatched_gt_list) %in% matched_gt$V1,]
    unmatched_gt_list <- unmatched_gt_list[rowSums(is.na(unmatched_gt_list)) != ncol(unmatched_gt_list), ]
    unmatched_gt_list <- cbind(rownames(unmatched_gt_list), unmatched_gt_list)
    unmatched_gt_list <- melt(data.table(unmatched_gt_list), id.vars = 'V1')
    unmatched_gt_list <- unmatched_gt_list[!is.na(unmatched_gt_list$value),]
    colnames(unmatched_gt_list) <- c("variant", "cell", "allele")
    
    write.csv(matched_gt, file.path(outputdir, paste0(donorid, "_matched_gt.csv")), row.names = F)
    write.csv(unmatched_gt_list, file.path(outputdir, paste0(donorid, "_unmatched_gt.csv")), row.names = F)
    
    informative_variants_cells <- merge(matched_gt, unmatched_gt_list, by.x = c("V1", "dominant"), 
                                        by.y=c("variant", "allele"))
    colnames(informative_variants_cells)[1] <- "variant"
    num_informative_variants <- informative_variants_cells %>% group_by(cell) %>% summarise(matched=n())
    if (nrow(unmatched_gt_list[!unmatched_gt_list$cell %in% num_informative_variants$cell,])>0){
      print(unmatched_gt_list[!unmatched_gt_list$cell %in% num_informative_variants$cell,])
    }
    representative_variant_list[[donorid]] <- unique(informative_variants_cells$variant)
    write.table(unique(informative_variants_cells$variant), file.path(outputdir, paste0(donorid, "_informative_variants.csv")), 
                row.names = F, col.names = F)
    
  }
  
  representative_variant <- Reduce(union, representative_variant_list)
  representative_variant_all <- data.frame(matrix(nrow = length(representative_variant), ncol = length(donors)+1)) 
  colnames(representative_variant_all) <- c("Variant", donors)
  representative_variant_all[["Variant"]] <- representative_variant
  for (donorid in donors){
    representative_variant_all[representative_variant_all$Variant %in% representative_variant_list[[donorid]],donorid]<- TRUE
  }
  representative_variant_all[is.na(representative_variant_all)] <- FALSE
  
  
  upset <- ComplexUpset::upset(representative_variant_all, donors, width_ratio = 0.3, height_ratio = 0.9, 
                               stripes = 'white', max_degree=1,
                               name="Number of variants occurring in each donors", 
                               set_sizes=(upset_set_size() + 
                                            geom_text(aes(label = ..count.., size=3), hjust = 0, 
                                                      stat = 'count', color = 'white', size = 2.5) +
                                            theme(axis.text.x = element_text(angle = 90),
                                                  text = element_text(size = 10))),
                               base_annotations=list('Intersection size'= intersection_size() )
  )
  ggsave(file.path(args$outputdir, "representative_variants_upset.png"))
  representative_variant_single <- representative_variant_all[rowSums(representative_variant_all[,-1]) == 1,]
  representative_variant_single <- separate(representative_variant_single, col= "Variant", 
                                            into = c("chr", "pos"), sep = "_")
  write.table(representative_variant_single[,c("chr", "pos")], quote = F, col.names=F,sep='\t',
              file.path(args$outputdir, "representative_variants_default.csv"), row.names = F)
}

if (args$findVariants == 'True' | args$findVariants == 'vireo'){
  if (startsWith(best_method1, 'vireo')){
    write.table(best_method1, file.path(args$outputdir, "best_method_vireo.txt"), sep = "\t", row.names = F, col.names = F, quote = F)
    
  } else{
    stop("Vireo is not the best method1 for donor matching!")
  }
  
  vireo_result_dir <- file.path(args$vireo_parent_dir, best_method1)

  representative_variant <- list.files(vireo_result_dir, "filtered_variants.tsv", full.names = T, recursive = T)[1]
  representative_variant <- fread(representative_variant)
  representative_variant <- separate(representative_variant, col= "variants", 
                                            into = c("chr", "pos"), sep = "_", extra = 'drop')
  write.table(representative_variant[,c("chr", "pos")], quote = F, col.names=F,sep='\t',
              file.path(args$outputdir, "representative_variants_vireo.csv"), row.names = F)
  
  
  } 
