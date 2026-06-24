#!/usr/bin/env Rscript

################################################
################################################
## Functions to handle Nextflow input         ##
################################################
################################################

string_to_null <- function(x, val = "[]") if (x == val) NULL else x

check_files <- function(args) {

  files <- c(
    "best_intersect_assignment_after_match",
    "cell_genotype",
    "variants_vireo",
    "result_csv"
  )

  for (f in files) {
    if (!is.null(args[[f]])){
      if(!file.exists(args[[f]])){
        stop(sprintf("'%s' does not exist.", f))
      }
      if(dir.exists(args[[f]])){
        stop(sprintf("'%s' is a directory but must be a file.", f))
      }
    }
  }
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

args <- list(
    # File inputs
    best_intersect_assignment_after_match = '$best_intersect_assignment_after_match',
    cell_genotype = string_to_null('$cell_genotype'),
    variants_vireo = string_to_null('$variants_vireo'),
    result_csv = string_to_null('$demultiplexing_result'), # only to have all barcodes

    # second in puts
    variant_count = as.numeric('$variant_count'),
    variant_pct = as.numeric('$variant_pct'),

    # others
    prefix = '$prefix' # Prefix name for output files.
)

check_files(args)

# Configure output precision
options(digits=5)

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(data.table)
library(ComplexUpset)
library(tidyverse)
library(vcfR)

################################################
################################################
## Main Process                               ##
################################################
################################################

# set TRUE to see print outputs for debugging
debugging <- FALSE

if (!is.null(args\$cell_genotype)) {

  result_csv <-
    fread(
      args\$result_csv,
      stringsAsFactors = FALSE,
      na.strings = c(NA_character_, "")
    )

  result_merge_new <-
    fread(file.path(args\$best_intersect_assignment_after_match),
      header = T
    )

  best_method1 <- names(result_merge_new)[2]
  best_method2 <- names(result_merge_new)[3]

  result_merge_new\$match <-
    result_merge_new[[best_method1]] == result_merge_new[[best_method2]]
  matched <- result_merge_new[result_merge_new\$match, ]
  unmatched <-
    result_csv[!result_csv\$Barcode %in% matched\$Barcode, ]\$Barcode
  cell_genotype_vcf <- read.vcfR(args\$cell_genotype)
  cell_genotype_vcf_gt <-
    extract.gt(cell_genotype_vcf,
      element = "GT",
      as.numeric = TRUE
    )
  donors <- sort(unique(matched[[best_method1]]))

  representative_variant_list <-
    vector(mode = "list", length = length(donors))
  representative_variant_list <-
    setNames(representative_variant_list, donors)

  for (donorid in donors) {
    matched_barcode <-
      matched[matched[[best_method1]] == donorid]\$Barcode
    matched_gt_list <- cell_genotype_vcf_gt[, matched_barcode]
    matched_gt_list <-
      matched_gt_list[rowSums(is.na(matched_gt_list)) != ncol(matched_gt_list), ]
    matched_gt <-
      as.data.frame(matrix(nrow = nrow(matched_gt_list)))

    # how many cells show homozygous reference
    matched_gt\$ref <- rowSums(matched_gt_list == 0, na.rm = TRUE)
    # how many cells show heterozygous or homozygous alternative
    matched_gt\$alt <- rowSums(matched_gt_list != 0, na.rm = TRUE)
    # chromosome positions
    matched_gt\$V1 <- rownames(matched_gt_list)
    matched_gt\$count <- matched_gt\$ref + matched_gt\$alt
    matched_gt\$pct <-
      matched_gt\$alt / (matched_gt\$ref + matched_gt\$alt)
    matched_gt\$dominant <- ifelse(matched_gt\$pct > 0.5, 1, 0)

    # variant_pct has to be in a range between [0,5;1[
    # 0.9 for example means we only keep variabts with higher than 90% or lower than 10% frequency
    matched_gt <- matched_gt[(matched_gt\$pct >= args\$variant_pct |
      matched_gt\$pct <= (1 - args\$variant_pct)), ]
    matched_gt <-
      matched_gt[matched_gt\$count >= args\$variant_count, ]

    unmatched_gt_list <- cell_genotype_vcf_gt[, unmatched]
    unmatched_gt_list <-
      unmatched_gt_list[rownames(unmatched_gt_list) %in% matched_gt\$V1, ]
    unmatched_gt_list <-
      unmatched_gt_list[rowSums(is.na(unmatched_gt_list)) != ncol(unmatched_gt_list), ]
    unmatched_gt_list <-
      cbind(rownames(unmatched_gt_list), unmatched_gt_list)
    unmatched_gt_list <-
      melt(data.table(unmatched_gt_list), id.vars = "V1")
    unmatched_gt_list <-
      unmatched_gt_list[!is.na(unmatched_gt_list\$value), ]
    colnames(unmatched_gt_list) <- c("variant", "cell", "allele")


    outputdir <- file.path(donorid)
    if (!dir.exists(outputdir)) {
      dir.create(outputdir, recursive = TRUE)
    }

    write.csv(matched_gt,
      file.path(outputdir, paste0(args\$prefix, "_", donorid, "_matched_gt.csv")),
      row.names = FALSE
    )
    write.csv(unmatched_gt_list,
      file.path(outputdir, paste0(args\$prefix, "_", donorid, "_unmatched_gt.csv")),
      row.names = FALSE
    )

    informative_variants_cells <-
      merge(
        matched_gt,
        unmatched_gt_list,
        by.x = c("V1", "dominant"),
        by.y = c("variant", "allele")
      )
    colnames(informative_variants_cells)[1] <- "variant"
    num_informative_variants <- informative_variants_cells %>%
      group_by(cell) %>%
      summarise(matched = n())
    if (nrow(unmatched_gt_list[!unmatched_gt_list\$cell %in% num_informative_variants\$cell, ]) > 0 | debugging) {
      print(unmatched_gt_list[!unmatched_gt_list\$cell %in% num_informative_variants\$cell, ])
    }
    representative_variant_list[[donorid]] <-
      list(unique(informative_variants_cells\$variant))
    write.table(
      unique(informative_variants_cells\$variant),
      file.path(outputdir,
        paste0(args\$prefix, "_", donorid, "_informative_variants.csv")
      ),
      row.names = FALSE,
      col.names = FALSE
    )
  }

  representative_variant <-
    rbindlist(representative_variant_list, idcol = "donor")
  colnames(representative_variant)[2] <- "variant"
  representative_variant_df <-
    dcast(data = representative_variant, variant ~ donor, length)
  write.csv(
    representative_variant_df,
    file.path(paste0(args\$prefix, "_all_representative_variants.csv"))
  )

  upset <- ComplexUpset::upset(
    representative_variant_df,
    donors,
    width_ratio = 0.45,
    height_ratio = 0.9,
    stripes = "white",
    max_degree = 1,
    name = "Number of donor-specific variants",
    set_sizes = (
      upset_set_size() +
        geom_text(
          aes(label = ..count.., size = 3),
          hjust = -0.1,
          stat = "count",
          color = "white",
          size = 2.3
        ) +
        theme(
          axis.text.x = element_text(angle = 90),
          text = element_text(size = 10)
        )
    ),
    base_annotations = list("Intersection size" = intersection_size())
  )
  ggsave(file.path(paste0(args\$prefix, "_donor_specific_variants_upset.png")))

  # variants that are unique to exactly one donor
  representative_variant_single <-
    representative_variant_df[rowSums(representative_variant_df[, -1]) == 1, ]
  representative_variant_single <-
    separate(
      representative_variant_single,
      col = "variant",
      into = c("chr", "pos"),
      sep = "_"
    )
  write.table(
    representative_variant_single[, c("chr", "pos")],
    quote = FALSE,
    col.names = FALSE,
    sep = "\t",
    row.names = FALSE,
    file.path(paste0(args\$prefix, "_donor_specific_variants.csv"))
  )
}

if (!is.null(args\$variants_vireo)) {

  representative_variant <- fread(args\$variants_vireo)
  representative_variant <- separate(
    representative_variant,
    col = "variants",
    into = c("chr", "pos"),
    sep = "_",
    extra = "drop"
  )
  write.table(
    representative_variant[, c("chr", "pos")],
    quote = FALSE,
    col.names = FALSE,
    sep = "\t",
    row.names = FALSE,
    file.path(paste0(args\$prefix, "_vireo_variants.csv"))
  )
}

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
data_table.version <- as.character(packageVersion('data.table'))
complexUpset.version <- as.character(packageVersion('ComplexUpset'))
tidyverse.version <- as.character(packageVersion('tidyverse'))
vcfR.version <- as.character(packageVersion('vcfR'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    r-complexupset:', complexUpset.version),
        paste('    r-data.table:', data_table.version),
        paste('    r-tidyverse:', tidyverse.version),
        paste('    r-vcfr:', vcfR.version)
    ),
'versions.yml')
