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
    result_csv = '$demultiplexing_result',
    barcode = '$barcode_whitelist',
    ndonor = as.numeric('$meta.n_samples'),
    cell_genotype = '$cell_genotype',
    vireo_parent_dir = '$vireo_parent_dir',

    # second in puts
    method1 = string_to_null('$match_donor_method1'),
    method2 = string_to_null('$match_donor_method2'),
    findVariants = as.logical('$findVariants'),
    variant_count = as.numeric('$variant_count'),
    variant_pct = as.numeric('$variant_pct'),

    # others
    prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'), # Prefix name for output files.
    outputdir = ""
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

# read assignment_all csv
result_csv <- NULL
min_cell <- 0
if (file.exists(args\$result_csv) && !dir.exists(args\$result_csv)) {
  result_csv <-
    fread(
      args\$result_csv,
      stringsAsFactors = FALSE,
      na.strings = c(NA_character_, "")
    )
}

# remove barcode that are not in the whitelist
if (!is.null(args\$barcode)) {
  barcode_whitelist <- fread(args\$barcode,
    header = FALSE,
    stringsAsFactors = FALSE
  )\$V1
  result_csv <-
    result_csv[result_csv\$Barcode %in% barcode_whitelist, ]
}


# finds all columns in the CSV that contain at least one real donor label (not “negative” or “doublet”), and returns their column names
colname_with_singlet <-
  colnames(result_csv %>% select_if(~ any(. != "negative" &
    . != "doublet")))
colname_with_singlet <-
  colname_with_singlet[colname_with_singlet != "Barcode"]

if (length(colname_with_singlet) < 2) {
  stop("Please choose more methods to run donor matching!")
}

hashing_methods <-
  c(
    "demuxem",
    "htodemux",
    "multiseq",
    "hashsolo",
    "hasheddrops",
    "bff_raw",
    "bff_cluster",
    "bff_consensuscall",
    "gmmdemux"
  )
genetic_methods <-
  c("demuxlet", "freemuxlet", "vireo", "scsplit", "souporcell")


if (!is.null(args\$method1) && !is.null(args\$method2)) {

    for (m in c("method1", "method2")) {
    method <- get(paste0("args\$", m))
    if (!any(startsWith(colname_with_singlet, method))) {
        warning(sprintf(
        "⚠️ %s ('%s') couldn't find at least one singlet. Ensure that the method you chose has at least one real donor label (not 'negative' or 'doublet') in one of the tasks.",
        tools::toTitleCase(m), method
        ))
    }
    }

    method1_all <- colname_with_singlet[startsWith(colname_with_singlet, args\$method1)]
    method2_all <- colname_with_singlet[startsWith(colname_with_singlet, args\$method2)]

} else {

  # get all column names that are genetic
  genetics_all <-
    Filter(function(x) {
      any(sapply(genetic_methods, function(y) {
        grepl(y, x)
      }))
    }, colname_with_singlet)

  # get all column names that are hashing
  hashing_all <-
    Filter(function(x) {
      any(sapply(hashing_methods, function(y) {
        grepl(y, x)
      }))
    }, colname_with_singlet)


  # Build pairs of methods that we want to compare in the for-loop

  # Match between genetics- and hashing-based methods
  if (length(hashing_all) > 0 && length(genetics_all) > 0) {
    all_methods_pair <-
      expand.grid(genetics = genetics_all, hashing = hashing_all)
    method1_all <- as.character(all_methods_pair\$genetics)
    method2_all <- as.character(all_methods_pair\$hashing)
  }

  # Compare only within hashing methods
  else if (length(hashing_all) > 0) {
    method_pair <- combn(hashing_all, 2)
    method1_all <- method_pair[1, ]
    method2_all <- method_pair[2, ]
  }

  # Compare only within genetics methods
  else if (length(genetics_all) > 0) {
    method_pair <- combn(genetics_all, 2)
    method1_all <- method_pair[1, ]
    method2_all <- method_pair[2, ]
  }
}

best_result <- 0
best_method1 <- "None"
best_method2 <- "None"
num_trial <- 1

result_record <- data.frame(
  best_method1 = character(),
  best_method2 = character(),
  score = numeric(),
  matched_donor = numeric(),
  remain_na = logical(),
  stringsAsFactors = FALSE
)

if (is.null(method1_all) || is.null(method2_all)) {
  stop("No method was found in the CSV file!")
}

for (i in 1:length(method1_all)) {

  # extract the pair of methods we would like to compare now
  method1 <- method1_all[i]
  method2 <- method2_all[i]

  # put the hashing method on the second place
  if (grepl(paste(hashing_methods, collapse = "|"), method1) &&
    (grepl(paste(genetic_methods, collapse = "|"), method2))) {
    hash_method <- method1
    method1 <- method2
    method2 <- hash_method
  }

  if(debugging){
    print(paste0("Comapring ", method1, " and ", method2))
  }

  outputdir <- file.path(paste0(method1, "_vs_", method2))
  if (!dir.exists(outputdir)) {
    dir.create(outputdir, recursive = TRUE)
  }

  filename_prefix <- paste0(args\$prefix,"_",method1, "_vs_", method2)

  method1_res <- convert2binary(result_csv, method1, min_cell)
  method2_res <- convert2binary(result_csv, method2, min_cell)
  if (is.null(method1_res) || is.null(method2_res)) {
    next
  }

  # Extract barcodes classified as singlets by both methods.
  # This meaning of intersect is not true for  edge cases
  # where a method assigned only one singlet label (see convert2binary if-statement).
  intersect_barcode <-
    intersect(rownames(method1_res), rownames(method2_res))
  if (length(intersect_barcode) == 0) {
    next
  }
  method1_res <-
    method1_res[rownames(method1_res) %in% intersect_barcode, , drop = FALSE]
  method2_res <-
    method2_res[rownames(method2_res) %in% intersect_barcode, , drop = FALSE]

  # correlation matrix with donor x donor
  correlation_res <- try(
    {
      apply(method1_res, 2, function(x) {
        apply(method2_res, 2, function(y) {
          return(cor.test(x, y)[["estimate"]][["cor"]])
        })
      })
    },
    silent = TRUE
  )
  # Skip this method pair if correlation calculation failed
  if (inherits(correlation_res, "try-error")) {
    cat("Failed to calculate phi coefficient")
    next
  }

  if (is.vector(correlation_res)) {
    correlation_res <- t(as.data.frame(correlation_res))
    rownames(correlation_res) <- colnames(method2_res)
  }
  write.csv(
    correlation_res,
    file.path(outputdir, paste0(filename_prefix, "_correlation_res.csv"))
  )

  match_score <- 0
  matched_donor <- 0
  geno_match <-
    as.data.frame(matrix(nrow = ncol(correlation_res), ncol = 3))
  colnames(geno_match) <- c("Method1", "Method2", "Correlation")
  geno_match\$Method1 <- colnames(correlation_res)

  for (id in geno_match\$Method1) {
    # Checks if the max value is finite and
    # if the current donor pair is a mutual best match between both methods
    if (!is.infinite(-max(correlation_res[, id], na.rm = TRUE)) &&
      max(correlation_res[, id], na.rm = TRUE) ==
        max(correlation_res[which.max(correlation_res[, id]), ], na.rm = TRUE)) {
      geno_match[which(geno_match\$Method1 == id), 2:3] <-
        c(
          rownames(correlation_res)[which.max(correlation_res[, id])],
          max(correlation_res[, id], na.rm = TRUE)
        )
      match_score <-
        match_score + max(correlation_res[, id], na.rm = TRUE)
      matched_donor <- matched_donor + 1
    } else {
      geno_match[which(geno_match\$Cluster1_ID == id)] <-
        c("unassigned", NA)
    }
  }

  write.table(
    geno_match[, 1:2],
    file.path(outputdir, paste0(filename_prefix, "_donor_match.csv")),
    row.names = FALSE,
    col.names = FALSE,
    sep = " ",
    quote = FALSE
  )

  # save concordance heatmap
  if (!all(is.na(correlation_res))) {
    newCols <- colorRampPalette(grDevices::rainbow(nrow(geno_match)))
    annoCol <- newCols(nrow(geno_match))
    names(annoCol) <- colnames(correlation_res)
    annoCol <- list(category = annoCol)
    correlation_res <-
      correlation_res[!is.na(row.names(correlation_res)), , drop = FALSE]
    correlation_res <-
      correlation_res[order(as.numeric(row.names(correlation_res))), , drop = FALSE]
    pheatmap(
      correlation_res,
      treeheight_row = FALSE,
      treeheight_col = FALSE,
      display_numbers = TRUE,
      angle_col = "45",
      number_color = "white",
      fontsize = 12,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      width = 7,
      height = 5,
      filename = file.path(outputdir, paste0(filename_prefix, "_concordance_heatmap.png"))
    )
  }

  if (grepl(paste(hashing_methods, collapse = "|"), method2) &&
    grepl(paste(genetic_methods, collapse = "|"), method1)) {
    remain_na <- (matched_donor != args\$ndonor)
    match_score <- match_score / args\$ndonor

    if (match_score > best_result && !remain_na) {
      write.table(
        geno_match[, 1:2],
        file.path(paste0(args\$prefix,"_best_donor_match.csv")),
        row.names = FALSE,
        col.names = FALSE,
        sep = " ",
        quote = FALSE
      )
      best_method1 <- method1
      best_method2 <- method2
      best_result <- match_score

    }

    new_record <-
      c(method1, method2, match_score, matched_donor, remain_na)
    result_record[num_trial, ] <- new_record
    num_trial <- num_trial + 1

    result_merge <- select(result_csv, "Barcode", method1, method2)
    result_merge_new <- result_merge

    # replace donor ID's from the genetic assignment with HTO of hashing
    for (i in 1:nrow(geno_match)) {
      result_merge_new[[method1]] <- replace(
        result_merge_new[[method1]],
        result_merge[[method1]] == geno_match\$Method1[i],
        geno_match\$Method2[i]
      )
    }

    # only retain barcodes that are in the intersect (definition of intersect see above)
    result_merge_new_intersect <-
      result_merge_new[result_merge_new\$Barcode %in% intersect_barcode, ]

    write.csv(
      result_merge_new,
      file.path(outputdir,paste0(filename_prefix, "_all_assignment_after_match.csv")),
      row.names = FALSE
    )

    write.csv(
      result_merge_new_intersect,
      file.path(outputdir, paste0(filename_prefix, "_intersect_assignment_after_match.csv")),
      row.names = FALSE
    )

    if (best_result == match_score) {
      write.csv(
        result_merge_new,
        file.path(paste0(args\$prefix,"_best_all_assignment_after_match.csv")),
        row.names = FALSE
      )

      write.csv(
        result_merge_new_intersect,
        file.path(paste0(args\$prefix,"_best_intersect_assignment_after_match.csv")),
        row.names = FALSE
      )
    }
  }
}

# TODO what is if there is more than one best match between methods?
if (best_method1 != "None" && best_method2 != "None" && debugging) {
  print(
    paste0(
      "Best method pair: ",
      best_method1,
      " and ",
      best_method2,
      " with score ",
      best_result
    )
  )
  print("------------------------------------------------------------------")
}

if (nrow(result_record) > 1) {
  write.csv(result_record,
    row.names = FALSE,
    file.path(paste0(args\$prefix,"_score_record.csv"))
  )
}

# TODO findVariants = true not implemented yet

if (args\$findVariants == "True" || args\$findVariants == "default") {
  if (startsWith(best_method1, "vireo")) {
    write.table(
      best_method1,
      file.path(args\$outputdir, "best_method_vireo.txt"),
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  } else {
    stop("Vireo is not the best method for donor matching!")
  }
  outputdir <-
    file.path(args\$outputdir, paste0(best_method1, "_vs_", best_method2))
  outputdir_variant <- file.path(outputdir, "variant_filtering")
  ifelse(!dir.exists(outputdir_variant),
    dir.create(outputdir_variant),
    FALSE
  )
  result_merge_new <-
    fread(file.path(outputdir, "intersect_assignment_after_match.csv"),
      header = T
    )
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

    write.csv(matched_gt,
      file.path(outputdir_variant, paste0(donorid, "_matched_gt.csv")),
      row.names = FALSE
    )
    write.csv(unmatched_gt_list,
      file.path(outputdir_variant, paste0(donorid, "_unmatched_gt.csv")),
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
      file.path(
        outputdir_variant,
        paste0(donorid, "_informative_variants.csv")
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
    file.path(args\$outputdir, "all_representative_variant_df.csv")
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
  ggsave(file.path(args\$outputdir, "donor_specific_variants_upset.png"))
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
    file.path(args\$outputdir, "donor_specific_variants.csv")
  )
}

if (args\$findVariants == "True" || args\$findVariants == "vireo") {
  if (startsWith(best_method1, "vireo")) {
    write.table(
      best_method1,
      file.path(args\$outputdir, "best_method_vireo.txt"),
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  } else {
    stop("Vireo is not the best method1 for donor matching, variants can not be filtered!")
  }

  vireo_result_dir <- file.path(args\$vireo_parent_dir, best_method1)

  representative_variant <-
    list.files(
      vireo_result_dir,
      "filtered_variants.tsv",
      full.names = TRUE,
      recursive = TRUE
    )[1]
  representative_variant <- fread(representative_variant)
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
    file.path(args\$outputdir, "representative_variants_vireo.csv")
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

################################################
################################################
################################################
################################################
