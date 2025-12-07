#!/usr/bin/env Rscript

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
## Functions to handle Nextflow input         ##
################################################
################################################

string_to_null <- function(x, val = "[]") if (x == val) NULL else x

check_files <- function(args) {

  files <- c(
    "barcode",
    "result_csv"
  )

  for (f in files) {
    if(!file.exists(args[[f]])){
      stop(sprintf("'%s' does not exist.", f))
    }
    if(dir.exists(args[[f]])){
      stop(sprintf("'%s' is a directory but must be a file.", f))
    }
  }
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

    # second in puts
    method1 = string_to_null('$match_donor_method1'),
    method2 = string_to_null('$match_donor_method2'),

    # others
    prefix = '$prefix', # Prefix name for output files.
    outputdir = ""
)

check_files(args)

# Configure output precision
options(digits=5)

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(pheatmap)
library(data.table)
library(tidyverse)

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

# if there is more than one best match between methods the last best_method will be printed
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

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
pheatmap.version <- as.character(packageVersion('pheatmap'))
data_table.version <- as.character(packageVersion('data.table'))
tidyverse.version <- as.character(packageVersion('tidyverse'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    r-data.table:', data_table.version),
        paste('    r-pheatmap:', pheatmap.version),
        paste('    r-tidyverse:', tidyverse.version),
    ),
'versions.yml')
