#!/usr/bin/env Rscript

################################################
################################################
## Functions to handle Nextflow input         ##
################################################
################################################

#' Check for Non-Empty, Non-Whitespace String
#'
#' This function checks if the input is non-NULL and contains more than just whitespace.
#' It returns TRUE if the input is a non-empty, non-whitespace string, and FALSE otherwise.
#'
#' @param input A variable to check.
#' @return A logical value: TRUE if the input is a valid, non-empty, non-whitespace string; FALSE otherwise.

is_valid_string <- function(input) {
    !is.null(input) && nzchar(trimws(input))
}

# Helper function for NULL condition
string_to_null <- function(x, val = "[]") if (x == val) NULL else x
null_to_string <- function(x, val = "NULL") if (is.null(x)) val else x

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

string_to_logical <- function(input) {
  if (input == "FALSE") {
    FALSE
  } else if (input == "TRUE") {
    TRUE
  } else {
    stop(paste0(input, " is not a valid logical. Use 'FALSE' or 'TRUE'."))
  }
}

################################################
################################################
## Functions for the script                   ##
################################################
################################################

convert2binary <- function(result_csv, method_name, min_cell) {
    #' Convert categorical donor assignments from a method into a binary (one-hot encoded) matrix of cells vs donors.
    #' Filters out invalid labels ("negative", "doublet", NA) if at least two different singlets assigments exist.
    #' Returns NULL if the number of valid cells is below the specified threshold.

  method_assign <- result_csv %>% select(all_of(c("Barcode", method_name)))
  donor_id <- setdiff(
    unique(method_assign[[method_name]]),
    c(NA, "negative", "doublet")
  )
  method_assign <-
    method_assign[method_assign[[method_name]] %in% donor_id, ]
  if (nrow(method_assign) < min_cell) {
    return(NULL)
  }
  if (length(unique(method_assign[[method_name]])) == 1) {
    method_assign_binary <-
      as.data.frame(matrix(0, nrow = nrow(result_csv), ncol = 1),
        row.names = result_csv\$Barcode
      )
    colnames(method_assign_binary) <-
      c(unique(method_assign[[method_name]]))
    method_assign_binary[rownames(method_assign_binary) %in% method_assign\$Barcode, ] <-
      1
  } else {
    method_assign_binary <-
      data.frame(model.matrix(~ method_assign[[method_name]] - 1, data = method_assign))
    names(method_assign_binary) <- sort(donor_id)
    rownames(method_assign_binary) <- method_assign\$Barcode
  }
  return(method_assign_binary)
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# Set defaults and classes
args <- list(
    # File inputs
    best_intersect_assignment_after_match = '$best_intersect_assignment_after_match',
    cell_genotype = '$cell_genotype',
    variants_vireo = '$variants_vireo',
    result_csv = '$demultiplexing_result', # only to have all barcodes

    # second in puts
    variant_count = as.numeric('$variant_count'),
    variant_pct = as.numeric('$variant_pct'),
    # others
    prefix = '$prefix' # Prefix name for output files.
)
opt_types <- lapply(args, class)

# Apply parameter overrides
args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{
        # Handle special cases for logicals
        if (opt_types[[ao]] == "logical") {
            opt[[ao]] <- string_to_logical(args_opt[[ao]])
        } else if (! is.null(opt[[ao]])){
            # Preserve classes from defaults where possible
            opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        } else {
            opt[[ao]] <- args_opt[[ao]]
        }
    }
}

# Configure output precision
options(digits=5)

# Check if file exists
# TODO check if files exist (not necessary for now)
# if (! file.exists(seuratObj)){
#     stop(paste0(seuratObj, ' is not a valid file'))
# }

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(pheatmap)
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

# TODO findVariants = true not implemented yet

if (args\$cell_genotype != "[]") {
  result_csv <- NULL
  if (file.exists(args\$result_csv) && !dir.exists(args\$result_csv)) {
    result_csv <-
      fread(
        args\$result_csv,
        stringsAsFactors = FALSE,
        na.strings = c(NA_character_, "")
      )
  }

  result_merge_new <-
    fread(file.path(args\$best_intersect_assignment_after_match),
      header = T
    )

  best_method1 <- names(result_merge_new)[2]
  best_method2 <- names(result_merge_new)[3]

  # outputdir <-
  #   file.path(args\$outputdir, paste0(best_method1, "_vs_", best_method2))
  # outputdir_variant <- file.path(outputdir, "variant_filtering")
  # ifelse(!dir.exists(outputdir_variant),
  #   dir.create(outputdir_variant),
  #   FALSE
  # )

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
    matched_gt\$ref <- rowSums(matched_gt_list == 0, na.rm = TRUE)
    matched_gt\$alt <- rowSums(matched_gt_list != 0, na.rm = TRUE)
    matched_gt\$V1 <- rownames(matched_gt_list)
    matched_gt\$count <- matched_gt\$ref + matched_gt\$alt
    matched_gt\$pct <-
      matched_gt\$alt / (matched_gt\$ref + matched_gt\$alt)
    matched_gt\$dominant <- ifelse(matched_gt\$pct > 0.5, 1, 0)
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
    if (nrow(unmatched_gt_list[!unmatched_gt_list\$cell %in% num_informative_variants\$cell, ]) > 0) {
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
    file.path(paste0(args\$prefix, "_all_representative_variant_df.csv"))
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
    file.path(paste0(args\$prefix, "_donor_match_representative_variants.csv"))
  )
}

if (args\$variants_vireo != "[]") {

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
    file.path(paste0(args\$prefix, "_vireo_representative_variants.csv"))
  )
}

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
pheatmap.version <- as.character(packageVersion('pheatmap'))
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
        paste('    r-pheatmap:', pheatmap.version),
        paste('    r-tidyverse:', tidyverse.version),
        paste('    r-vcfr:', vcfR.version)
    ),
'versions.yml')
