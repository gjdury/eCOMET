########################################################################################
# Define functions for data import and normalization
########################################################################################
#' Import mzmine feature data and metadata to create a mmo object
#' @description
#' This function reads mzmine feature data and metadata from specified directories,
#' to initiate a mmo object containing feature data and metadata
#'
#' @param mzmine_dir Path to the mzmine feature data CSV file
#' @param metadata_dir Path to the metadata CSV file. The metadata file should contain columns: sample, group, mass
#' @return A mmo object containing the feature data and metadata

GetMZmineFeature <- function(mzmine_dir, metadata_dir, group_col){
  mmo <- list()
  data <- read.csv(mzmine_dir)
  #create feature column
  data <- data |> mutate(feature = paste(.data$mz, .data$rt, sep = '_'))
  area_columns <- grep("datafile.+\\.mzML.area", names(data), value = TRUE)
  feature_df <- data |> select(.data$id, .data$feature, all_of(area_columns))
  feature_df$id <- gsub(" ", "", feature_df$id)
  mmo$feature_data <- feature_df
  metadata <- read.csv(metadata_dir)
  metadata$sample <- paste('datafile.', metadata$sample, '.area', sep = '')
  mmo$metadata <- metadata
  mmo$pairwise <- data.frame(feature = mmo$feature_data$feature, id = mmo$feature_data$id)
  if (missing(group_col) || !(group_col %in% colnames(metadata))) {
    stop("group_col must be provided and must exist in the metadata file.")
  }  
  mmo$metadata$group <- as.factor(mmo$metadata[[group_col]])
  print('MMO object created.')
  print(paste0('Feature number: ', nrow(mmo$feature_data)))
  print(paste0(nrow(mmo$metadata), ' samples in ', length(unique(mmo$metadata$group)), ' groups'))
  return(mmo)
}


#' Switch the group column in the mmo object
#'
#' @description
#' This function switches the group column in the metadata of the mmo object to a new specified column.
#' The new group column must exist in the metadata file.
#' @param mmo The mmo object
#' @param new_group_col The name of the new group column in the metadata file
#' @return The mmo object with the updated group column
SwitchGroup <- function(mmo, new_group_col) {
  if (missing(new_group_col) || !(new_group_col %in% colnames(mmo$metadata))) {
    stop("new_group_col must be provided and must exist in the metadata file.")
  }
  mmo$metadata$group <- as.factor(mmo$metadata[[new_group_col]])
  print(paste0('Group column switched to ', new_group_col))
  print(paste0(length(unique(mmo$metadata$group)), ' groups in total'))
  return(mmo)
}

#' Adding annotation from SIRIUS to the mmo object
#'
#' @description
#' This function reads SIRIUS structure identification and formula summary files,
#' and adds the annotations to the mmo object.
#' @param mmo The mmo object
#' @param canopus_structuredir Path to the SIRIUS structure_identification.tsv file
#' @param canopus_formuladir Path to the SIRIUS canopus_formula_summary.tsv file
#' @return The mmo object with SIRIUS annotations added
AddSiriusAnnot <- function(mmo, canopus_structuredir, canopus_formuladir){
  structure_identifications <- readr::read_tsv(canopus_structuredir, show_col_types = FALSE)
  structure_identifications$mappingFeatureId <- gsub(" ", "", structure_identifications$mappingFeatureId)
  canopus_formula_summary <- readr::read_tsv(canopus_formuladir, show_col_types = FALSE)
  canopus_formula_summary$mappingFeatureId <- gsub(" ", "", canopus_formula_summary$mappingFeatureId)
  siriused_ids <- unique(union(structure_identifications$mappingFeatureId, canopus_formula_summary$mappingFeatureId))
  sirius_df <- mmo$feature_data |> select(.data$id, .data$feature)
  sirius_df <- sirius_df |>
  left_join(structure_identifications, by = c("id" = "mappingFeatureId")) |>
  left_join(canopus_formula_summary, by = c("id" = "mappingFeatureId"))
  mmo$sirius_annot <- sirius_df
  print('SIRIUS annotation added to mmo$sirius_annot')
  return(mmo)
}

#' Add custom annotations to an mmo object
#'
#' @description
#' Match features to a custom DB by m/z (ppm) and RT (minutes) tolerances and
#' attach a list-column of candidate compound IDs per feature.
#'
#' @param mmo An `mmo` object created by `GetMZmineFeature()`.
#' @param DB_file CSV path with at least columns `compound`, `mz`, `rt`.
#' @param mztol m/z tolerance in ppm (default 5).
#' @param rttol RT tolerance in minutes (default 0.5).
#' @return The same `mmo` object with `mmo$custom_annot` (id, feature, custom_annot).
#' @export
#'



AddCustomAnnot <- function(mmo, DB_file, mztol = 5, rttol = 0.5) {
  DB <- read.csv(DB_file, stringsAsFactors = FALSE)

  DB <- dplyr::mutate(
    DB,
    mz = as.numeric(.data$mz),
    rt = as.numeric(.data$rt)
  )

  # Parse "feature" (formatted like "mz_rt") into numeric mz/rt columns
  feature_annot <- mmo$feature_data |>
    tidyr::separate(
      col = "feature",
      into = c("mz", "rt"),
      sep = "_",
      remove = FALSE,
      convert = TRUE
    )

  # For each feature, collect matching DB compounds within tolerances
  annotated_features <- dplyr::mutate(
    feature_annot,
    custom_annot = purrr::map2(.data$mz, .data$rt, function(mzi, rti) {
      dplyr::filter(
        DB,
        abs(mzi - .data$mz) / mzi * 1e6 <= mztol,
        abs(rti - .data$rt) <= rttol
      ) |>
        dplyr::pull(.data$compound)
    })
  )

  mmo$custom_annot <- dplyr::select(
    annotated_features,
    .data$id, .data$feature, .data$custom_annot
  )

  message("Custom annotation added to mmo$custom_annot using ", DB_file)
  mmo
}


#' Replace zero and NA values in the mmo object
#'
#' This function replaces zero and NA values in the feature data of the mmo object.
#'
#' @param mmo The mmo object
#' @param method The method to use for replacement. Options are 'one' (replace with 1) or 'half_mean' (replace with half of the smallest non-zero value in the row)
ReplaceZero <- function(mmo, method = 'one') {
  df <- mmo$feature_data
  df[] <- apply(df, 1, function(row) {
    # Convert the row to numeric, ignoring non-numeric columns
    numeric_row <- as.numeric(row[-c(1, 2)])  # Skip 'id' and 'feature' columns
    # Get the smallest non-zero, non-NA value in the row
    smallest_value <- min(numeric_row[numeric_row > 0], na.rm = TRUE)
    # Replace 0 and NA with half of the smallest_value
    row[-c(1, 2)] <- sapply(numeric_row, function(x) {
      if (is.na(x) || x == 0) {
        if (method == 'one') {
          return(1)
        } else if (method == 'half_mean') {
          return(smallest_value / 2)
        }
      } else {
        return(x)
      }
    })

    return(row)
  }) |>
    t() |>
    as.data.frame()  # Convert back to dataframe
  mmo$feature_data <- df
  print(paste('Missing values were filled with', method))
  return(mmo)
}
#' Use mass in the metadata file to normalize the peak area
#'
#' This function normalizes the peak area in the feature data of the mmo object by the mass of each sample, provided in the metadata.
#' The feature data is replaced by (original value * mean mass) / sample mass.
#'
#' @param mmo The mmo object
#' @return The mmo object with normalized feature data (mmo$feature_data)
MassNormalization <- function(mmo){
  normalized_df <- mmo$feature_data
  metadata <- mmo$metadata
  mean_mass <- mean(mmo$metadata$mass)
  for (sample_col in colnames(mmo$feature_data)[-c(1,2)]) {
    sample_metadata <- metadata[metadata$sample == sample_col, ]
    mass <- sample_metadata$mass
    normalized_df[[sample_col]] <- as.numeric(mmo$feature_data[[sample_col]])*mean_mass/mass
  }
  mmo$feature_data <- normalized_df
  print("Peak area are normalized by sample mass")
  return(mmo)
}

#' Log-normalize the peak area in the mmo object
#'
#' This function applies log2 transformation to the peak area in the feature data of the mmo object.
#' @param mmo The mmo object
#' @return The mmo object with log-normalized feature data (mmo$log)
LogNormalization <- function(mmo){
  feature_data_only <- mmo$feature_data[,-(1:2)]
  log_data <- log2(feature_data_only)
  log_df <- cbind(mmo$feature_data[, 1:2], log_data)
  mmo$log <- log_df
  print('Log-normalized values were added to mmo$log')
  return(mmo)
}

#' Mean-center the peak area in the mmo object
#'
#' This function applies mean-centering to the peak area in the feature data of the mmo object.
#'
#' @param mmo The mmo object
#' @return The mmo object with mean-centered feature data (mmo$meancentered)
MeancenterNormalization <- function(mmo){
  feature_data_only <- mmo$feature_data[,-(1:2)]
  mean_centered_data <- t(apply(feature_data_only, 1, function(x) x - mean(x, na.rm = TRUE)))
  mean_centered_df <- cbind(mmo$feature_data[, 1:2], mean_centered_data)
  mmo$meancentered <- mean_centered_df
  print('Meancentered values were added to mmo$meancentered')
  return(mmo)
}

#' Z-normalize the peak area in the mmo object
#'
#' This function applies Z-score normalization to the peak area in the feature data of the mmo object.
#'
#' @param mmo The mmo object
#' @return The mmo object with Z-normalized feature data (mmo$zscore)
ZNormalization <- function(mmo){
  feature_data_only <- mmo$feature_data[,-(1:2)]
  zscore_df <- t(apply(feature_data_only, 1, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
  zscore_df <- cbind(mmo$feature_data[, 1:2], zscore_df)
  mmo$zscore <- zscore_df
  print('Z-score values were added to mmo$zscore')
  return(mmo)
}


#' Add chemical distance matrices to the mmo object
#'
#' This function reads cosine, DREAMS, and MS2DeepScore molecular networking outputs from mzmine,
#' then transform the similarity to distance and adds the dissimilarity matrices to the mmo object.
#'
#' @param mmo The mmo object
#' @param cos_dir Path to the cosine similarity CSV file
#' @param dreams_dir Path to the DREAMS similarity CSV file
#' @param m2ds_dir Path to the MS2DeepScore similarity CSV file
#' @return The mmo object with dissimilarity matrices added (mmo$cos.dissim, mmo$dreams.dissim, mmo$m2ds.dissim)
AddChemDist <- function(mmo, cos_dir = NULL, dreams_dir = NULL, m2ds_dir = NULL) {
  .require_pkg("data.table")
  add_dissim_matrix <- function(mmo, sim_dir, slot_name) {
    # Fast read
    sim <- data.table::fread(sim_dir, col.names = c("cluster1", "cluster2", "metric", "similarity", "etc"))
    sim$dissimilarity <- 1 - sim$similarity

    # Ensure cluster IDs are clean characters
    sim$cluster1 <- trimws(as.character(sim$cluster1))
    sim$cluster2 <- trimws(as.character(sim$cluster2))
    
    # Build mapping
    clusters <- unique(c(sim$cluster1, sim$cluster2))
    cluster_index <- setNames(seq_along(clusters), clusters)
    
    # Map clusters to integer indices
    i <- cluster_index[sim$cluster1]
    j <- cluster_index[sim$cluster2]
    
    # Check for problems
    if (anyNA(i) || anyNA(j)) {
      bad <- unique(c(sim$cluster1[is.na(i)], sim$cluster2[is.na(j)]))
      stop("Cluster IDs not found in mapping: ", paste(bad, collapse = ", "))
    }
    
    # Preallocate full dense matrix
    n <- length(clusters)
    dissim_mat <- matrix(1, nrow = n, ncol = n)
    diag(dissim_mat) <- 0
    
    # Fill dissimilarities
    dissim_mat[cbind(i, j)] <- sim$dissimilarity
    dissim_mat[cbind(j, i)] <- sim$dissimilarity
    
    # Add names
    dimnames(dissim_mat) <- list(clusters, clusters)
    
    # Reduce memory footprint
    mode(dissim_mat) <- "single"
    
    mmo[[slot_name]] <- dissim_mat
    message(slot_name, " added to mmo")
    return(mmo)
  }
  
  if (!is.null(cos_dir))    mmo <- add_dissim_matrix(mmo, cos_dir,    "cos.dissim")
  if (!is.null(dreams_dir)) mmo <- add_dissim_matrix(mmo, dreams_dir, "dreams.dissim")
  if (!is.null(m2ds_dir))   mmo <- add_dissim_matrix(mmo, m2ds_dir,   "m2ds.dissim")
  
  if (is.null(cos_dir) && is.null(dreams_dir) && is.null(m2ds_dir)) {
    stop("Please provide at least one valid directory.")
  }
  
  return(mmo)
}



#' Reorder samples in the mmo object based on group order
#' 
#' This function reorders the samples in the mmo object based on a specified group order.
#' The function updates the order of samples in the feature data, log-normalized data, z-score data, and mean-centered data.
#' 
#' @param mmo The mmo object
#' @param group_order A vector specifying the desired order of groups
#' @return The mmo object with reordered samples
ReorderGroups <- function(mmo, group_order) {
  metadata <- mmo$metadata
  # Get sample names in the specified group order
  ordered_samples <- unlist(lapply(group_order, function(g) metadata$sample[metadata$group == g]))
  # Reorder columns: id, feature, then ordered samples
  mmo$feature_data <- mmo$feature_data |>
    dplyr::select(.data$id, .data$feature, all_of(ordered_samples))
  mmo$log <- mmo$log |>
    dplyr::select(.data$id, .data$feature, all_of(ordered_samples))
  mmo$zscore <- mmo$zscore |>
    dplyr::select(.data$id, .data$feature, all_of(ordered_samples))
  mmo$meancentered <- mmo$meancentered |>
    dplyr::select(.data$id, .data$feature, all_of(ordered_samples))
  return(mmo)
}

########################################################################################
# Define functions for supporting analysis
########################################################################################

#' Get feature data from the mmo object, with normalization options
#'
#' This function retrieves the feature data from the mmo object based on the specified normalization method.
#'
#' @param mmo The mmo object
#' @param normalization The normalization method to use. Options are 'None', 'Log', 'Meancentered', or 'Z'
#' @return The feature data corresponding to the specified normalization method
GetNormFeature <- function(mmo, normalization = 'None'){
  if (normalization == 'None'){
    feature <- mmo$feature_data
  } else if (normalization == 'Log'){
    feature <- mmo$log
  } else if (normalization == 'Meancentered'){
    feature <- mmo$meancentered
  } else if (normalization == 'Z'){
    feature <- mmo$zscore
  } else {
    print('The normalization should be None, Log, Meancentered, or Z')
  }
  return(feature)
}


#' Get the distance matrix from the mmo object based on the specified distance metric
#'
#' This function retrieves the distance matrix from the mmo object based on the specified distance metric.
#'
#' @param mmo The mmo object
#' @param distance The distance metric to use. Options are 'dreams', 'cosine', or 'm2ds'
#' @return The distance matrix corresponding to the specified distance metric
GetDistanceMat <- function(mmo, distance = 'dreams'){
  if (distance == 'dreams'){
    distance_matrix <- mmo$dreams.dissim
  } else if (distance == 'cosine'){
    distance_matrix <- mmo$cos.dissim
  } else if (distance == 'm2ds'){
    distance_matrix <- mmo$m2ds.dissim
  }
  return(distance_matrix)
}

#' Convert feature names to IDs in the mmo object
#'
#' This function converts feature names to their corresponding IDs in the mmo object.
#'
#' @param mmo The mmo object
#' @param feature_names A vector of feature names to convert
#' @return A vector of feature IDs corresponding to the input feature names
FeatureToID <- function(mmo, feature_names) {
  feature_data <- mmo$feature_data
  feature_ids <- feature_data |>
    filter(.data$feature %in% feature_names) |>
    select(.data$feature, .data$id)
  # Match the order of feature_names
  feature_ids <- feature_ids$id[match(feature_names, feature_ids$feature)]
  return(feature_ids)
}

#' Convert feature IDs to names in the mmo object
#'
#' This function converts feature IDs to their corresponding names in the mmo object.
#'
#' @param mmo The mmo object
#' @param feature_ids A vector of feature IDs to convert
#' @return A vector of feature names corresponding to the input feature IDs
IDToFeature <- function(mmo, feature_ids) {
  feature_data <- mmo$feature_data
  feature_names <- feature_data |>
    filter(.data$id %in% feature_ids) |>
    select(.data$id, .data$feature)
  # Match the order of feature_ids
  feature_names <- feature_names$feature[match(feature_ids, feature_names$id)]
  return(feature_names)
}

#' Calculate group means from the mmo object
#'
#' This function calculates the mean feature values for each group in the mmo object.
#'
#' @param mmo The mmo object
#' @param normalization The normalization method to use. Options are 'None', 'Log', 'Meancentered', or 'Z'
#' @param filter_feature Boolean to filter features based on a provided list (default: FALSE)
#' @param feature_list A vector of feature names to filter (default: NULL)
#' @param filter_group Boolean to filter groups based on a provided list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
#' @return A data frame containing the mean feature values for each group
GetGroupMeans <- function(mmo, normalization = 'None', filter_feature = FALSE, feature_list = NULL, filter_group = FALSE, group_list = NULL) {
  feature_data <- GetNormFeature(mmo, normalization = normalization)
  metadata <- mmo$metadata

  # Melt the feature data to long format
  long_feature_data <- feature_data |>
    tidyr::pivot_longer(
      cols = -c(.data$id, .data$feature),                 # all columns except id and feature
      names_to = "sample",                    # old column names go here
      values_to = "feature_value"             # values go here
    )

  colnames(long_feature_data) <- c('id', 'feature', 'sample', 'feature_value')

  # Merge with metadata to get group information
  merged_data <- merge(long_feature_data, metadata[, c('sample', 'group')], by = 'sample')
  if (filter_group == TRUE){
    merged_data <- merged_data |> filter(.data$group %in% group_list)
  }
  # Calculate group means
  group_means <- merged_data |>
    dplyr::group_by(.data$group, .data$id) |>
    dplyr::summarise(mean_value = mean(.data$feature_value, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = .data$group,
      values_from = .data$mean_value
    )

  if (filter_feature) {
    group_means <- group_means |>
      dplyr::filter(.data$id %in% FeatureToID(mmo, feature_list))
  }
  return(group_means)
}

#' Calculate log2 fold change for a given control group
#'
#' This function calculates the log2 fold change for each group compared to a specified control group.
#' @param group_means A data frame containing the mean feature values for each group
#' @param control_group The name of the control group to compare against
#' @return A data frame with log2 fold change values for each group compared to the control group
GetLog2FoldChange <- function(group_means, control_group) {
  control_means <- group_means[[control_group]]
  fold_change <- group_means |>
    mutate(across(-.data$id, ~ log2(. / control_means)))

  return(fold_change)
}


#' Perform ANOVA and Tukey's HSD test on the mmo object
#'
#' This function performs ANOVA and Tukey's HSD test on the feature data of the mmo object,
#'
#' @param df The data frame containing the feature data and metadata
#' @param formula The formula for the ANOVA test, e.g., "feature ~ group"
#' @return A list containing the ANOVA results, Tukey's HSD results, Tukey's significance letters, and Dunnett's test results
anova_tukey_dunnett <- function(df, formula) {
  .require_pkg("DescTools")
  .require_pkg("multcompView")
  aov_res <- aov(as.formula(formula), data = df)
  tukey_res <- TukeyHSD(aov_res)
  tukey_sig <- multcompView::multcompLetters4(aov_res, tukey_res)
  dunnett_res <- DescTools::DunnettTest(as.formula(formula), data = df)
  return(list(aov_res = aov_res, tukey_res = tukey_res, tukey_sig = tukey_sig, dunnett_res = dunnett_res))
}

#' Write ANOVA and Tukey's HSD results to a CSV file
#'
#' This function writes the results of ANOVA and Tukey's HSD test to a CSV file.
#'
#' @param anova_data A list containing the results of ANOVA and Tukey's HSD test
#' @param outdir The output directory where the results will be saved
#' @param way The type of ANOVA test to perform. Options are 'oneway' or 'twoway'
write_anova <- function(anova_data, outdir, way='oneway'){
  way_num <- switch(way, oneway = 1, twoway = 3)
  # Perform ANOVA and Tukey HSD
  aov_res <- anova_data$aov_res
  tukey_res <- anova_data$tukey_res
  tukey_sig <- anova_data$tukey_sig[way_num]
  dunnett_res <- anova_data$dunnett_res

  # Save ANOVA and Tukey HSD results
  anova_df <- as.data.frame(summary(aov_res)[[1]])
  anova_df$Comparison <- rownames(anova_df)
  tukey_df <- as.data.frame(tukey_res[way_num])
  tukey_df$Comparison <- rownames(tukey_df)
  sig_letter <- as.data.frame(unlist(tukey_sig))
  sig_letter$Comparison <- rownames(sig_letter)
  dunnett_df <- as.data.frame(dunnett_res[[1]])
  dunnett_df$comp <- rownames(dunnett_df)

  # Create a combined results data frame
  combined_df <- dplyr::bind_rows(
    dplyr::tibble(Test = "ANOVA", anova_df),
    dplyr::tibble(Test = "Tukey", tukey_df),
    dplyr::tibble(Test = 'sig', sig_letter),
    dplyr::tibble(Test = 'Dunnett', dunnett_df)
  )
  readr::write_csv(combined_df, file = outdir)
}

#' Perform PERMANOVA and pairwise comparisons
#' 
#' This function performs PERMANOVA on the given data and metadata, with options for filtering groups.
#' It also conducts post-hoc pairwise comparisons and adjusts p-values for multiple testing.
#' 
#' @param data A data frame or distance matrix for PERMANOVA
#' @param metadata A data frame containing sample metadata, including a 'group' column
#' @param mode The mode of the input data: 'data' for raw data or 'distance' for a distance matrix
#' @param filter_group Boolean to filter groups based on a provided list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
#' @param permutations The number of permutations for the PERMANOVA test (default: 5000)
#' @return A list containing the PERMANOVA results, raw pairwise comparison results, and a matrix of adjusted p-values
permanova_stat <- function(data, metadata, mode, filter_group = FALSE, group_list = NULL, permutations = 5000){ 
  .require_pkg("vegan")
  .require_pkg("stringr")
  .require_pkg("pairwiseAdonis")
  #Perform PERMANOVA
  if (mode == 'data'){
    dist_mat <- dist(data)
  } else if (mode == 'distance'){
    dist_mat <- as.dist(data)
  } else {
    stop("Invalid mode for adonis2. Please use 'data' or 'distance'")
  }
  dist_mat <- dist(data)
  metadata <- metadata[match(rownames(data), metadata$sample), ]
  if (filter_group == TRUE){
    metadata <- metadata |> filter(.data$group %in% group_list)
    dist_mat <- as.dist(as.matrix(dist_mat)[metadata$sample, metadata$sample])
  }
  permanova_res <- vegan::adonis2(dist_mat ~ group, data = metadata, permutations = permutations)
  #Post-hoc
  pairwise_permanova <- pairwiseAdonis::pairwise.adonis2(dist_mat ~ group, data = metadata, permutations = permutations)
  # Multiple test correction for pairwise p-values
  pvals <- sapply(pairwise_permanova[-1], function(x) x[1,5])
  pvals_adj <- p.adjust(pvals, method = "BH")
  # Add adjusted p-values to pairwise_permanova
  for (i in seq_along(pvals_adj)) {
    pairwise_permanova[[i+1]]$padj <- pvals_adj[i]
  }
  #Organize posthoc results for visual
  pairwise_matrix <- matrix(NA, nrow = length(unique(metadata$group)), ncol = length(unique(metadata$group)))
  rownames(pairwise_matrix) <- unique(metadata$group)
  colnames(pairwise_matrix) <- unique(metadata$group)
  for (i in 2:length(pairwise_permanova)){
    group1 <- str_split(names(pairwise_permanova)[i], '_vs_')[[1]][1]
    group2 <- str_split(names(pairwise_permanova)[i], '_vs_')[[1]][2]
    pval <- pairwise_permanova[[i]][1,6]
    pairwise_matrix[group2, group1] <- pval
  }
  return(list(permanova_res = permanova_res, pairwise_raw = pairwise_permanova, pairwise_matrix = pairwise_matrix))
}

########################################################################################
# Define functions for pairwise comparison and visualization
########################################################################################

#' Perform pairwise comparison between two groups in the mmo object
#'
#' This function performs pairwise comparison between two groups in the mmo object,
#' calculating log2 fold change and adjusted p-values for given comparison of two groups.
#' The function adds the results to the mmo$pairwise data frame.
#'
#' @param mmo The mmo object
#' @param group1 The name of the nominator group
#' @param group2 The name of the denominator group
#' @param correction The method for multiple comparison correction. Options are 'BH', 'holm', 'bonferroni', etc. Inherits from p.adjust()
PairwiseComp <- function(mmo, group1, group2, correction = 'BH'){
  feature <- mmo$feature_data
  metadata <- mmo$metadata
  #Get sample names
  group1_samples <- metadata |> filter(.data$group == group1) |> pull(sample)
  group2_samples <- metadata |> filter(.data$group == group2) |> pull(sample)
  #Get data from the samples
  group1_data <- feature |> select(.data$id, .data$feature, all_of(group1_samples))
  group2_data <- feature |> select(.data$id, .data$feature, all_of(group2_samples))
  #Make empty column
  log2FC <- numeric(nrow(feature))
  pval <- numeric(nrow(feature))
  #Pairwise comparison
  for (i in 1:nrow(feature)){
    group1_value <- as.numeric(group1_data[i, -c(1,2)])
    group2_value <- as.numeric(group2_data[i, -c(1,2)])

    group1_mean <- mean(group1_value, na.rm = TRUE)
    group2_mean <- mean(group2_value, na.rm = TRUE)
    log2FC[i] <- log2(group2_mean/group1_mean)

    pval[i] <- tryCatch(
      expr = {
        p <- t.test(group1_value, group2_value, na.rm = TRUE)$p.value
        p
      },
      error = function(e)
      {
        return(1)
      }
    )

    # ttest <- t.test(group1_value, group2_value, na.rm = TRUE)

    # pval[i] <- ttest$p.value
  }
  padj <- p.adjust(pval, method = correction)
  #Store in results
  results <- data.frame(
    log2FC= log2FC,
    padj = padj
  )
  names(results) <- c(paste(group1, "vs", group2, "log2FC", sep = "_"), paste(group1, "vs", group2, "padj", sep = "_"))
  #Add pairwise results to the mmo object
  mmo$pairwise <- cbind(mmo$pairwise, results)
  print(paste(group2, '/', group1, 'comparison was completed'))
  return(mmo)
}


#' Get lists of DAMs (Differentially Accumulated Metabolites) for each comparison
#'
#' This function generates lists of upregulated and downregulated DAMs for each pairwise comparison in the mmo object.
#'
#' @param mmo The mmo object with pairwise comparison matrix
#' @param fc_cutoff The threshold of log2 fold change to be considered significant (default: 0.5849625, which is log2(1.5))
#' @param pval_cutoff The threshold of adjusted p-value to be considered significant (default: 0.05)
#' @return A list containing two lists: DAMs_up and DAMs_down
GetDAMs <- function(mmo, fc_cutoff = 0.5849625, pval_cutoff = 0.05) {
  # Generate the list of comparisons automatically by looking up mmo$pairwise
  comparison_columns <- colnames(mmo$pairwise)
  log2FC_columns <- grep("log2FC", comparison_columns, value = TRUE)
  comparisons <- unique(sub("log2FC", "", log2FC_columns))
  comparisons <- sub("_$", "", comparisons)  # Remove trailing underscore from comparisons
  # Make list of DAMs for up and downregulation for each comparison
  DAMs_up <- list()
  DAMs_down <- list()
  for (comp in comparisons) {
    group1 <- strsplit(comp, "_vs_")[[1]][1]
    group2 <- strsplit(comp, "_vs_")[[1]][2]
    DAMs_up[[paste(comp, "up", sep = ".")]] <- filter(mmo$pairwise, get(paste(comp, "log2FC", sep = "_")) > fc_cutoff & get(paste(comp, "padj", sep = "_")) < pval_cutoff)$feature
    DAMs_down[[paste(comp, "down", sep = ".")]] <- filter(mmo$pairwise, get(paste(comp, "log2FC", sep = "_")) < -fc_cutoff & get(paste(comp, "padj", sep = "_")) < pval_cutoff)$feature
  }
  names(DAMs_up) <- paste(comparisons, "up", sep = ".")
  names(DAMs_down) <- paste(comparisons, "down", sep = ".")
  return(list(DAMs_up = DAMs_up, DAMs_down = DAMs_down))
}

# Plot log2FC and -log(pval) from pairwise comparison.  PairwiseComp(mmo, 'group1', 'group2') should be precended
# mmo : mmo object with pairwise comparison matrix
# topk : the number of features to sshow labels
# pthr : the threshold of adjusted p-value to be considered significant
# outdir : plot name
# height : plot height (<20)
# width  : plot width (<20)

#' Volcano plot for visualizing differential metabolite analysis results
#'
#' This function generates a volcano plot using data from mmo$pairwise,
#' highlighting upregulated and downregulated features based on log2 fold change and adjusted p-value
#'
#' @param mmo The mmo object with pairwise comparison matrix
#' @param comp The comparison to visualize, e.g., 'group1_vs_group2
#' @param topk The number of top features to label in the plot (default: 10)
#' @param pthr The threshold of adjusted p-value to be considered significant (default: 0.05)
#' @param outdir The output file path for the volcano plot (default: 'volcano.png')
#' @param height The height of the output plot in inches (default: 5)
#' @param width The width of the output plot in inches (default: 5)
VolcanoPlot <- function(mmo, comp, topk = 10, pthr = 0.05, outdir = 'volcano.png', height = 5, width = 5){
  .require_pkg("ggrepel")
  VolData <- mmo$pairwise |> select(.data$feature,all_of(c(paste(comp, 'log2FC', sep = '_'), paste(comp, 'padj', sep = '_'))))
  colnames(VolData) <- c('feature', 'log2FC', 'padj')
  VolData <- VolData |>
    mutate(
      Expression = dplyr::case_when(log2FC >= 1 & padj <= pthr ~ "Up-regulated",
                            log2FC <= -1 & padj <= pthr ~ "Down-regulated",
                            TRUE ~ "Not significant")
      )

  top_features <- dplyr::bind_rows(
    VolData |>
      filter(.data$Expression =='Up-regulated')  |>
      arrange(dplyr::desc(abs(.data$log2FC)), .data$padj) |>
      head(topk),
    VolData |>
      filter(.data$Expression == 'Down-regulated') |>
      arrange(dplyr::desc(abs(.data$log2FC)), .data$padj) |>
      head(topk)
  )

  volcano <- ggplot(VolData, aes(x = .data$log2FC, y = -log(.data$padj, 10))) +
    geom_point(aes(color = .data$Expression), size = 0.4)+
    xlab(expression("log"[2]*"FC")) +
    ylab(expression("-log"[10]*"FDR"))+
    scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    theme_classic()+
    ggrepel::geom_label_repel(data = top_features,
                    mapping = aes(.data$log2FC, -log(.data$padj,10), label = .data$feature),
                    size = 2)

  volcano
  ggsave(outdir, height = height, width = width)
}

########################################################################################
# Define functions for multivariate analysis
########################################################################################


#' PCA plot
#' 
#' This function performs PCA analysis and generates a PCA plot with optional filtering of features and groups.
#' It also conducts PERMANOVA and saves the results to CSV files.
#'
#' @param mmo The mmo object with feature data and metadata
#' @param color A vector of colors for the groups in the plot
#' @param outdir The output file path for the PCA plot (default: 'PCA')
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'Z')
#' @param filter_feature Boolean to filter features by feature_list (default: FALSE)
#' @param feature_list A vector of feature names to filter (default: NULL)
#' @param filter_group Boolean to filter groups by group_list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
#' @param label Boolean to indicate whether to label points with sample names (default: TRUE)
#' @return A PCA plot saved as a PDF file and PERMANOVA results saved as CSV files
#'
PCAplot <- function(mmo, color, outdir = 'PCA', normalization = 'Z', filter_feature = FALSE, feature_list = NULL, filter_group = FALSE, group_list = NULL, label = TRUE){
  .require_pkg("ggrepel")

  metadata <- mmo$metadata
  feature <- GetNormFeature(mmo, normalization)
  if (filter_feature == TRUE){
    feature <- feature |> filter(.data$feature %in% feature_list)
  }
  # Perform PCA on normalized feature data
  feature_data_pca <- feature[, -(1:2)]
  feature_data_pca <- t(feature_data_pca) # samples as rows, features as columns
  pca_res <- prcomp(feature_data_pca, scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x)
  pca_df$group <- metadata$group[match(rownames(pca_df), metadata$sample)]

  if (filter_group == TRUE){
    pca_df <- pca_df |> filter(.data$group %in% group_list)
  }
  if (label == TRUE){
    plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = rownames(pca_df))) +
      geom_point(size = 3) +
      geom_label_repel(aes(label = rownames(pca_df)), size = 3) +
      theme_classic() +
      labs(x = "PC1", y = "PC2") +
      scale_color_manual(values = custom_colors)+
      stat_ellipse(aes(group = group), level = 0.90)
  } else {
    plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
      geom_point(size = 3) +
      # geom_label_repel(aes(label = rownames(pca_df)), size = 3) +
      theme_classic() +
      labs(x = "PC1", y = "PC2") +
      scale_color_manual(values = custom_colors)+
      stat_ellipse(aes(group = group), level = 0.90)
  }
  plot
  ggsave(paste0(outdir, '.pdf'), width = 6, height = 6)

  permanova <- permanova_stat(feature_data_pca, metadata, mode = 'data', filter_group = filter_group, group_list = group_list)
  write.csv(permanova$permanova_res, paste0(outdir, '_permanova_results.csv'))
  write.csv(as.data.frame(permanova$pairwise_raw), paste0(outdir, '_pairwise_permanova_results.csv'))
  write.csv(as.data.frame(permanova$pairwise_matrix), paste0(outdir, '_pairwise_permanova_pvalue_matrix.csv'))
}


# PLS-DA plot
# colors : for scale_color_manual
# topk : number of features to show loading
# Choose which feature area to use for calculation among 'None' (mass-normalized), 'Log','Meancentered', and 'Z'

#' PLSDAplot
#'
# This function performs PLS-DA analysis and generates a PLS-DA plot with feature loadings.
#'
#' @param mmo The mmo object with feature data and metadata
#' @param color A vector of colors for the groups in the plot
#' @param topk The number of top features to display in the plot (default: 10)
#' @param outdir The output file path for the PLS-DA plot
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'Z')
#' @param filter_feature Boolean to filter features by feature_list (default: FALSE)
#' @param feature_list A vector of feature names to filter (default: NULL)
#' @param filter_group Boolean to filter groups by group_list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
PLSDAplot <- function(mmo, color, topk = 10, outdir, normalization = 'Z', filter_feature = FALSE, feature_list = NULL, filter_group = FALSE, group_list = NULL) {
  .require_pkg("caret")
  .require_pkg("ggrepel")

  metadata <- mmo$metadata
  #Get appropriate feature by normalization parameter
  feature <- GetNormFeature(mmo, normalization)

  # All feature or filtered feature
  if (filter_feature == TRUE){
    feature <- feature |> filter(feature %in% feature_list)
  }

  X <- t(as.matrix(feature[, -(1:2)]))
  Y <- c()
  for (col in colnames(feature)[-c(1,2)]){
    Y <- append(Y, metadata[metadata$sample == col, ]$group)
  }
  Y <- as.factor(Y)

  plsda_model <- caret::plsda(X, Y, ncomp = 2)
  scores <- plsda_model$scores[, 1:2]
  plsda_df <- data.frame(Comp1 = scores[, 1], Comp2 = scores[, 2], Group = Y)
  loadings <- plsda_model$loadings
  loadings_comp1 <- loadings[, 1]
  loadings_comp2 <- loadings[, 2]
  if (filter_feature == FALSE){
    loadings_df <- data.frame(Feature = mmo$feature_data$feature,
                            Comp1_Loading = loadings_comp1,
                            Comp2_Loading = loadings_comp2)
  } else {
    loadings_df <- data.frame(Feature = feature_list,
                            Comp1_Loading = loadings_comp1,
                            Comp2_Loading = loadings_comp2)
  }


  top_features <- loadings_df |>
  mutate(abs_loading_comp1 = abs(.data$Comp1_Loading),
         abs_loading_comp2 = abs(.data$Comp2_Loading)) |>
  arrange(dplyr::desc(.data$abs_loading_comp1 + .data$abs_loading_comp2)) |>
  head(topk)
  loading_scale <- 1
  if (topk > 0){
  loading_scale <- max(abs(scores))/(4*max(abs(top_features$Comp1_Loading)))
  }

  if (filter_group == TRUE){
    plsda_df <- plsda_df |> filter(.data$Group %in% group_list)
  }

  ggplot(plsda_df, aes(x = .data$Comp1, y = .data$Comp2, color = .data$Group)) +
  geom_point(size = 3) +
  theme_classic() +
  stat_ellipse(level = 0.90) +
  ggtitle("PLS-DA Plot") +
  labs(x = "Component 1", y = "Component 2") +
  scale_color_manual(values = .data$custom_colors) +
  theme(legend.position = "right")+
  geom_segment(data = top_features,
               aes(x = 0, y = 0, xend = .data$Comp1_Loading * loading_scale, yend = .data$Comp2_Loading * loading_scale),  # Scale the arrows
               arrow = grid::arrow(length = grid::unit(0.3, "cm")), color = "grey", linewidth = 1) +
  # Add labels for the top 10 features
  ggrepel::geom_text_repel(data = top_features,
            aes(x = .data$Comp1_Loading * loading_scale, y = .data$Comp2_Loading * loading_scale, label = .data$Feature),
            color = "black", vjust = 1.5, size = 3)

  ggsave(outdir, height = 6, width = 6)
  write.csv(loadings_df, 'PLSDA_loadings.csv')
  print(paste(normalization, '-normalized feature was used'))
}

#' Generate heatmap inputs from the mmo object
#'
#' This function generates heatmap inputs from the mmo object, including fold change or mean values,
#' distance matrix, and row labels for custom-annotated features.
#'
#' @param mmo The mmo object with sirius annotation and normalized data
#' @param filter_feature Boolean to filter features by feature_list (default: FALSE)
#' @param feature_list A vector of feature names to filter (default: NULL)
#' @param filter_group Boolean to filter groups by group_list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
#' @param summarize The summarization method to use. Options are 'fold_change' or 'mean' (default: 'fold_change')
#' @param control_group The group to use as control for fold change calculation (default: 'ctrl')
#' @param normalization The normalization method to use. Options are 'None', 'Log', 'Meancentered', or 'Z'
#' @param distance The distance metric to use. Options are 'dreams', 'cosine', or 'm2ds' (default: 'dreams')
#' @return A list containing the following elements:
#' - FC_matrix: A matrix of fold change or mean values
#' - dist_matrix: A distance matrix based on the specified distance metric
#' - row_label: A vector of row labels for custom-annotated features
#' - heatmap_data: A data frame containing the heatmap data with feature IDs and values
GenerateHeatmapInputs <- function(mmo, filter_feature = FALSE, feature_list = NULL,
                                filter_group = FALSE, group_list = NULL,
                                summarize = 'fold_change', control_group = 'ctrl',
                                normalization = 'None', distance = 'dreams') {
  # 12.1.1. Get summarized data (group mean or FC)
  if (filter_group){
    group_means <- GetGroupMeans(mmo, normalization = normalization, filter_group = TRUE, group_list = group_list)
  } else {
    group_means <- GetGroupMeans(mmo, normalization = normalization)
  }
  if (summarize == 'fold_change'){
    fold_change <- GetLog2FoldChange(group_means, control_group = control_group)
    heatmap_data <- fold_change
    heatmap_data[[control_group]] <- NULL
  } else if(summarize == 'mean'){
    heatmap_data <- group_means
  }
  # 12.1.2. Filter features
  # Determine distance metric
  distance_matrix <- GetDistanceMat(mmo, distance = distance)
  heatmap_data <- heatmap_data |> filter(id %in% rownames(distance_matrix)) # remove features not in distance matrix

  # make matrix for heatmap
  FC_matrix <- as.matrix(heatmap_data[,-1])
  rownames(FC_matrix) <- heatmap_data$id
  # Reorder the rows of distance_matrix to match the order of FC_matrix_
  distance_matrix <- distance_matrix[rownames(FC_matrix), rownames(FC_matrix)]
  dist_matrix <- as.dist(distance_matrix)

  row_label <- rownames(FC_matrix)
  if (filter_feature){
    filter_list <- feature_list
    filter_id <- FeatureToID(mmo, filter_list)
    filter_id <- filter_id[filter_id %in% rownames(distance_matrix)] # remove custom-annotated but not in the distance matrix
    filter_distance <- distance_matrix[filter_id, filter_id]
    heatmap_data <- heatmap_data |> filter(id %in% filter_id)

    # make matrix for heatmap
    FC_matrix <- as.matrix(heatmap_data[,-1])
    rownames(FC_matrix) <- heatmap_data$id


    # Reorder the rows of distance_matrix to match the order of FC_matrix_
    filter_distance <- filter_distance[rownames(FC_matrix), rownames(FC_matrix)]
    dist_matrix <- as.dist(filter_distance)
    #Label custm-annotated features
    row_label <- rownames(FC_matrix)
    for (i in 1:length(rownames(FC_matrix))){
      id <- rownames(FC_matrix)[i]
      custom_annot <- mmo$custom_annot$custom_annot[mmo$custom_annot$id == id]
      if (length(custom_annot[[1]]) > 0) {
        row_label[i] <- custom_annot[[1]]
      }
    }
  }
  return(list(FC_matrix = FC_matrix, dist_matrix = dist_matrix, row_label = row_label, heatmap_data = heatmap_data))
}


################################################
#Define enrichment analysis using Canopus-predicted terms
# mmo : the mmo object with sirius annotation and normalized
# list_test : a vector containing names of features to analyze
# pthr : the threshold for adjusted p-value to be considered significant
# sig : a logical vaue to show only significant terms or not
# term_level : the level of term to use for enrichment analysis
# representation : 'greater' for overrepresentation

#' Enrichment analysis for Canopus-predicted terms
#'
#' This function performs enrichment analysis for Canopus-predicted terms on a given list of features.
#'
#' @param mmo The mmo object with sirius annotation and normalized data
#' @param list_test A vector containing names of features to analyze
#' @param pthr The threshold for adjusted p-value to be considered significant (default: 0.1)
#' @param sig A logical value indicating whether to return only significant terms (default: TRUE)
#' @param term_level The level of term to use for enrichment analysis
#'               Options are 'NPC_pathway', 'NPC_superclass', 'NPC_class', 'ClassyFire_superclass', 'ClassyFire_class',
#'              'ClassyFire_subclass', 'ClassyFire_level5', or 'ClassyFire_most_specific' (default: 'NPC_pathway')
#' @param representation The representation type for enrichment analysis. Options are 'greater' for overrepresentation (default: 'greater')
#' @return A data frame containing the enrichment results, including term level, term name, subset count, total count, fold enrichment, p-value, and adjusted p-value (FDR)
#'
CanopusLevelEnrichmentAnal <- function(mmo,list_test, pthr = 0.1, sig=TRUE, term_level = 'NPC_pathway', representation = 'greater'){
  all_feature <- mmo$sirius_annot
  subset_feature <- mmo$sirius_annot |> filter(.data$feature %in% list_test)
  # print(paste('total features:', nrow(all_feature), 'list_test features:', nrow(subset_feature)))
  # Select the appropriate term level for enrichment analysis
  if (term_level == "NPC_pathway") {
    all_feature$classifications_split <- all_feature[[32]]
    subset_feature$classifications_split <- subset_feature[[32]]
  } else if (term_level == "NPC_superclass") {
    all_feature$classifications_split <- all_feature[[34]]
    subset_feature$classifications_split <- subset_feature[[34]]
  } else if (term_level == "NPC_class") {
    all_feature$classifications_split <- all_feature[[36]]
    subset_feature$classifications_split <- subset_feature[[36]]
  } else if (term_level == "ClassyFire_superclass") {
    all_feature$classifications_split <- all_feature[[38]]
    subset_feature$classifications_split <- subset_feature[[38]]
  } else if (term_level == "ClassyFire_class") {
    all_feature$classifications_split <- all_feature[[40]]
    subset_feature$classifications_split <- subset_feature[[40]]
  } else if (term_level == "ClassyFire_subclass") {
    all_feature$classifications_split <- all_feature[[42]]
    subset_feature$classifications_split <- subset_feature[[42]]
  } else if (term_level == "ClassyFire_level5") {
    all_feature$classifications_split <- all_feature[[44]]
    subset_feature$classifications_split <- subset_feature[[44]]
  } else if (term_level == "ClassyFire_most_specific") {
    all_feature$classifications_split <- all_feature[[46]]
    subset_feature$classifications_split <- subset_feature[[46]]
  } else {
    stop("Invalid term level. Please choose a valid term level.")
  }

  total_term_counts <- table(unlist(all_feature$classifications_split))
  subset_term_counts <- table(unlist(subset_feature$classifications_split))

  total_term_counts['None'] <- sum(is.na(all_feature$classifications_split))
  subset_term_counts['None'] <- sum(is.na(subset_feature$classifications_split))

  # Perform enrichment analysis using Fisher's exact test
  enrichment_results <- sapply(names(subset_term_counts), function(term) {
    contingency_matrix <- matrix(c(
      subset_term_counts[[term]],
      sum(subset_term_counts) - subset_term_counts[[term]],
      total_term_counts[[term]],
      sum(total_term_counts) - total_term_counts[[term]]
    ), nrow = 2, byrow = TRUE)
    fisher.test(contingency_matrix, alternative = representation)$p.value
  })
  # Adjust p-values for multiple testing
  adjusted_pvalues <- p.adjust(enrichment_results, method = "fdr")
  # Create a results dataframe
  results <- data.frame(
    term_level = term_level,
    term = names(enrichment_results),
    subsetcount = as.numeric(subset_term_counts[names(enrichment_results)]),
    totalcount = as.numeric(total_term_counts[names(enrichment_results)]),
    foldenrichment = (as.numeric(subset_term_counts[names(enrichment_results)]) / length(subset_feature))/(as.numeric(total_term_counts[names(enrichment_results)]) / nrow(all_feature)),
    pval = enrichment_results,
    fdr = adjusted_pvalues
  )
  results <- results |> filter(.data$term != 'None')
  # Filter for significantly enriched terms
  significant_terms <- results |>
    filter(.data$pval < pthr) |>
    arrange(.data$pval)
  if(sig==TRUE){
    return(significant_terms)
  }else{
    return(results)
  }
}

#' Generate a plot for enrichment analysis of Canopus-predicted terms
#'
#' This function generates a plot for enrichment analysis of Canopus-predicted terms,
#' showing fold enrichment, p-value, and subset count for each term level.
#'
#' @param mmo The mmo object with sirius annotation and normalized data
#' @param feature_list A vector containing names of features to analyze
#' @param pthr The threshold for adjusted p-value to be considered significant (default: 0.05)
#' @param outdir The output file path for the enrichment plot
#' @param height The height of the output plot in inches (default: 5)
#' @param width The width of the output plot in inches (default: 5)
CanopusListEnrichmentPlot <- function(mmo, feature_list, pthr = 0.05, outdir, height = 5, width = 5){
  term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  sig.canopus <- data.frame(term = character(),  term_level = character(),subsetcount = double(), totalcount = double(), foldenrichment = double(), pval = double(), fdr = double())
  for (term_level in term_levels){
    sig.canopus <- rbind(sig.canopus, CanopusLevelEnrichmentAnal(mmo, feature_list, pthr = pthr, sig = TRUE, term_level = term_level, representation = 'greater'))
  }
  sig.canopus <- sig.canopus |> arrange(dplyr::desc(.data$foldenrichment))
  ggplot(sig.canopus, aes(x = .data$foldenrichment, y = reorder(.data$term, .data$foldenrichment), color = -log(.data$pval), size = .data$subsetcount)) +
    geom_point() +
    scale_color_gradient(low = 'grey', high = 'red') +
    theme_classic()+
    facet_grid(term_level ~ ., scales = 'free_y', space = 'free', switch = 'y')+
    xlim(0,max(sig.canopus$foldenrichment+1))
    #facet_wrap(~term_level, ncol = 1, scales = 'free_y', strip.position = 'right', shrink = TRUE)

  ggsave(outdir, height = height, width = width)
}

#' Generate a plot for enrichment analysis of Canopus-predicted terms across multiple levels
#'
#' This function generates a plot for enrichment analysis of Canopus-predicted terms across multiple levels,
#' showing fold enrichment, p-value, and subset count for each term level.
#'
#' @param mmo The mmo object with sirius annotation and normalized data
#' @param feature_list A vector containing names of features to analyze
#' @param pthr The threshold for adjusted p-value to be considered significant (default: 0.05)
#' @param outdir The output file path for the enrichment plot
#' @param height The height of the output plot in inches (default: 5)
#' @param width The width of the output plot in inches (default: 5)
#' @param topn The number of top terms to display in the plot (default: 5)
CanopusListEnrichmentPlot_2 <- function(mmo, feature_list, pthr = 0.05, outdir, height = 5, width = 5, topn = 5){
  term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  sig.canopus <- data.frame(term = character(),  term_level = character(),subsetcount = double(), totalcount = double(), foldenrichment = double(), pval = double(), fdr = double())
  for (term_level in term_levels){
    sig.canopus <- rbind(sig.canopus, CanopusLevelEnrichmentAnal(mmo, feature_list, pthr = pthr, sig = TRUE, term_level = term_level, representation = 'greater'))
  }
  sig.canopus$term <- paste(sig.canopus$term, ';', sig.canopus$term_level)
  sig.canopus <- sig.canopus |> dplyr::slice_max(order_by = -.data$pval, n = topn)
  sig.canopus <- sig.canopus |> dplyr::arrange(dplyr::desc(.data$foldenrichment))
  ggplot(sig.canopus, aes(x = .data$foldenrichment, y = reorder(.data$term, .data$foldenrichment), color = -log(.data$pval), size = .data$subsetcount)) +
    geom_point() +
    scale_color_gradient(low = 'grey', high = 'red') +
    theme_classic()+
    xlim(0,max(sig.canopus$foldenrichment+1))+
    ylab('Chemical Class')

  ggsave(outdir, height = height, width = width)
}

#' Generate a plot for enrichment analysis of Canopus-predicted terms at a specific level
#'
#' This function generates a plot for enrichment analysis of Canopus-predicted terms at a specific level,
#' showing fold enrichment, p-value, and subset count for each term.
#'
#' @param mmo The mmo object with sirius annotation and normalized data
#' @param comp.list A list to analyze, where each element is a vector of feature names
#' @param term_level The level of term to use for enrichment analysis.
#'               Options are 'NPC_pathway', 'NPC_superclass', 'NPC_class',
#'              'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass',
#'             'ClassyFire_level5', or 'ClassyFire_most_specific' (default: 'NPC_pathway')
#' @param pthr The threshold for adjusted p-value to be considered significant (default: 0.1)
#' @param representation The representation type for enrichment analysis. Options are 'greater' for overrepresentation (default: 'greater')
#' @param prefix The prefix for output files (default: 'enrichment')
#' @param height The height of the output plot in inches (default: 5)
#' @param width The width of the output plot in inches (default: 5)
CanopusLevelEnrichmentPlot <- function(mmo = mmo, comp.list, term_level = 'NPC_pathway',pthr = 0.1, representation = 'greater', prefix = 'enrichment', height = 5, width = 5){
  df.EA <- data.frame()
  sig.terms <- c()
  for(list in names(comp.list)){
    # Calculate enrichment score for all terms
    res <- CanopusLevelEnrichmentAnal(mmo = mmo, list_test = comp.list[[list]],sig=FALSE, pthr = pthr, representation = representation, term_level = term_level)
    res <- res |> mutate(comp = list)
    df.EA <- dplyr::bind_rows(df.EA, res)
    # get terms that are at least once enriched in one comparison
    res.sig <- CanopusLevelEnrichmentAnal(mmo = mmo, list_test = comp.list[[list]],sig=TRUE, pthr = pthr, representation = representation, term_level = term_level)
    sig.terms <- append(sig.terms, res.sig$term)
  }
  sig.terms <- unique(sig.terms)
  df.EA.sig <- df.EA |> filter(.data$term %in% sig.terms)
  df.EA.sig <- df.EA.sig |>
    mutate(label = cut(
      .data$pval,
        breaks = c(0,0.001, 0.01, 0.05, 1),
        labels = c("***", "**", "*", "")
    ))

  enrichment_plot <- ggplot(data = df.EA.sig, aes(x = .data$comp, y = .data$term, label = .data$label))+
    geom_point(aes(size = .data$subsetcount, color = .data$pval))+
    geom_text()+
    scale_size_area(name = 'Count', max_size = 10)+
    scale_color_gradient2(low = 'red', high = 'grey', mid = 'grey', midpoint = 0.4)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))+
    xlab('Comparisons')+
    ylab('Chemical classes')
  enrichment_plot
  write.csv(df.EA, paste(prefix, '.csv'), row.names = FALSE)
  write.csv(df.EA.sig, paste(prefix, '_sig.csv'), row.names = FALSE)
  ggsave(paste(prefix, '.pdf'), width = width, height = height)
}

#' Generate a plot for enrichment analysis of Canopus-predicted terms across all levels
#'
#' This function generates a plot for enrichment analysis of Canopus-predicted terms across all levels,
#' showing fold enrichment, p-value, and subset count for each term level.
#' @param mmo The mmo object with sirius annotation and normalized data
#' @param comp.list A list to analyze, where each element is a vector of feature names
#' @param terms The terms to analyze. Options are 'all_terms', 'NPC', 'ClassyFire', or 'custom' (default: 'all_terms')
#' @param term_levels list of custom term levels to use
#' @param pthr The threshold for adjusted p-value to be considered significant (default: 0.1)
#' @param representation The representation type for enrichment analysis. Options are 'greater' for overrepresentation (default: 'greater')
#' @param prefix The prefix for output files (default: 'enrichment')
#' @param height The height of the output plot in inches (default: 10)
#' @param width The width of the output plot in inches (default: 8)
CanopusAllLevelEnrichmentPlot <- function(mmo = mmo, comp.list, terms = 'all_terms', term_levels = NULL, pthr = 0.1, representation = 'greater', prefix = 'enrichment', height = 10, width = 8){
  df.EA <- data.frame()
  sig.terms <- c()
  if(terms == 'all_terms'){
    term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  } else if (terms == 'NPC'){
    term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway')
  } else if (terms == 'ClassyFire'){
    term_levels = c('ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  } else if (terms == 'custom'){
    term_levels = term_levels
  }
  for(term_level in term_levels){
    for(list in names(comp.list)){
      # Calculate enrichment score for all terms
      res <- CanopusLevelEnrichmentAnal(mmo = mmo, list_test = comp.list[[list]],sig=FALSE, pthr = pthr, representation = representation, term_level = term_level)
      res <- res |> mutate(comp = list)
      df.EA <- dplyr::bind_rows(df.EA, res)
      # get terms that are at least once enriched in one comparison
      res.sig <- CanopusLevelEnrichmentAnal(mmo = mmo, list_test = comp.list[[list]],sig=TRUE, pthr = pthr, representation = representation, term_level = term_level)
      sig.terms <- append(sig.terms, res.sig$term)
    }
    sig.terms <- unique(sig.terms)
    df.EA.sig <- df.EA |> filter(.data$term %in% sig.terms)
    df.EA.sig <- df.EA.sig |>
      mutate(label = cut(
        .data$pval,
          breaks = c(0,0.001, 0.01, 0.05, 1),
          labels = c("***", "**", "*", "")
      ))
  }
  enrichment_plot <- ggplot(data = df.EA.sig, aes(x = .data$comp, y = .data$term, label = .data$label))+
    geom_point(aes(size = .data$subsetcount, color = .data$pval))+
    geom_text()+
    scale_size_area(name = 'Count', max_size = 10)+
    scale_color_gradient2(low = 'red', high = 'grey', mid = 'grey', midpoint = 0.4)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))+
    xlab('Comparisons')+
    ylab('Chemical classes')+
    facet_grid(term_level ~ ., scales = 'free_y', space = 'free', switch = 'y')
  enrichment_plot
  write.csv(df.EA, paste0(prefix, '.csv'), row.names = FALSE)
  write.csv(df.EA.sig, paste0(prefix, '_sig.csv'), row.names = FALSE)
  ggsave(paste0(prefix, '.pdf'), width = width, height = height)
}


#' Metabolite Set Enrichment Analysis (MSEA)
#' This function performs Metabolite Set Enrichment Analysis (MSEA) using the fgsea package.
#' It takes a ranked list of feature scores and tests for enrichment of metabolite sets based on Canopus-predicted terms.
#' The results are saved as a CSV file and a PDF plot.
#' @param mmo The mmo object with sirius annotation and normalized data
#' @param feature_name A vector of feature names corresponding to the feature scores
#' @param feature_score A vector of feature scores (e.g., log2 fold changes)
#' @param term_level The level of term to use for enrichment analysis.
#'               Options are 'NPC_pathway', 'NPC_superclass', 'NPC_class',
#'              'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass',
#'             'ClassyFire_level5', or 'ClassyFire_most_specific' (default: 'NPC_class')
#' @param pthr The threshold for adjusted p-value to be considered significant (default: 0.05)
#' @param prefix The prefix for output files (default: 'MSEA')
#' @param width The width of the output plot in inches (default: 8)
#' @param height The height of the output plot in inches (default: 12)
#' @param sig A logical value indicating whether to return only significant terms (default: FALSE)
MSEA <- function(mmo, feature_name, feature_score, term_level = 'NPC_class', pthr = 0.05, prefix = 'MSEA', width = 8, height = 12, sig = FALSE){
  # Create a named vector of feature scores
  ranked_list <- feature_score
  names(ranked_list) <- feature_name
  ranked_list <- sort(ranked_list, decreasing = TRUE)

  # Retrieve metabolite sets based on the specified term level
  if(term_level == 'NPC_class'){
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[[36]])
  } else if (term_level == 'NPC_superclass'){
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[[34]])
  } else if (term_level == 'NPC_pathway'){
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[[32]])
  } else if (term_level == "ClassyFire_superclass") {
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[[38]])
  } else if (term_level == "ClassyFire_class") {
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[[40]])
  } else if (term_level == "ClassyFire_subclass") {
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[[42]])
  } else if (term_level == "ClassyFire_level5") {
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[[44]])
  } else if (term_level == "ClassyFire_most_specific") {
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[[46]])
  } else {
    stop("Invalid term level. Please choose a valid term level.")
  }
  msea_res <- fgsea::fgsea(pathways = metabolite_sets,
                       stats    = ranked_list,
                       minSize  = 5,   # minimum number of features in a class
                       maxSize  = 1500,
                       nPermSimple = 10000)
  msea_res <- msea_res |> arrange(padj)
  readr::write_csv(msea_res, paste0(prefix,'_', term_level,'_results.csv'))
  if (sig) {
    msea_res <- msea_res |> filter(padj < pthr)
  }
  ggplot(msea_res, aes(x = reorder(pathway, NES), y = NES)) +
    geom_point(shape = 21, aes(color = padj < 0.05, size = size, fill = -log(padj)), stroke = 1) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    scale_fill_gradient(low = "grey", high = "red") +
    scale_color_manual(values = c("TRUE" = 'black', "FALSE" = 'white')) +
    guides(shape = "none") +
    labs(x = "Metabolite Class", y = "Normalized Enrichment Score (NES)", title = "MSEA Results", color = "-log10(padj)", size = "Set Size") +
    theme_classic() +
    theme(legend.position = "top", axis.text.y = element_text(size = 6)) 
  ggsave(paste0(prefix,'_', term_level,'.pdf'), width = width, height = height)
}


#' FeaturePerformanceRegression
#'
#' This function performs regression analysis of a specific feature against a phenotype performance in the metadata.
#' It can use linear mixed models (LMM), simple linear regression (LM), or Pearson correlation.
#'
#' @param mmo The mmo object with feature data and metadata
#' @param target The name of the feature to analyze
#' @param phenotype The name of the phenotype performance in the metadata
#' @param groups A vector of group names from the metadata containing performance data
#' @param model The type of regression model to use. Options are 'lmm' for linear mixed model, 'lm' for simple linear regression, or 'pearson' for Pearson correlation (default: 'lmm')
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'Z')
#' @param output The output file path for the regression plot

FeaturePerformanceRegression <- function(mmo, target, phenotype, groups, model = 'lmm', normalization = 'Z', output){
  .require_pkg("ggrepel")
  feature <- GetNormFeature(mmo, normalization)
  metadata <- mmo$metadata

  # Get phenotype performance from the metadata, get the feature value from the feature matrix, then combine
  phenotype.df <- data.frame(sample = metadata$sample, group = metadata$group, performance = metadata[,phenotype]) |> filter(.data$group %in% groups)
  feature_df <- data.frame(sample = colnames(feature[,-(1:2)]), feature_value = as.numeric(feature[feature$feature == target, -(1:2)]))
  combined_df <- merge(phenotype.df, feature_df, by='sample')

  # Perform linear mixed model or simple linear regression
  if (model == 'lmm'){
    fit <- lme4::lmer(combined_df$performance ~ combined_df$feature_value + (1|combined_df$group))
    p_value <- summary(fit)$coefficients[2, 5]
  } else if (model == 'lm'){
    fit <- lm(combined_df$performance ~ combined_df$feature_value)
    p_value <- summary(fit)$coefficients[2, 4]
  } else if (model == 'pearson'){
    pearson <- cor.test(combined_df$performance, combined_df$feature_value)
    p_value <- pearson[[3]]
  } else {
    stop("Invalid model type. Please use 'lmm' or 'lm' or 'pearson")
  }

  # Plot the fit using ggplot
  ggplot(combined_df, aes(x = .data$feature_value, y = .data$performance, color = .data$group)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    ggrepel::geom_text_repel(aes(label = sample), size = 2.5, show.legend = FALSE) +
    theme_classic() +
    labs(title = paste("Regression of", target, "against", phenotype, "performance"),
         x = "Feature Value",
         y = "Performance") +
    theme(legend.position = "right") +
    annotate("text", x = Inf, y = Inf, label = paste("p-value:", signif(p_value, digits = 4)),
             hjust = 1.1, vjust = 1.1, size = 3, color = "black")

  ggsave(output, height = 6, width = 6)
}




#' GetPerformanceFeatureRegression
#'
#' This function performs linear regression analysis of all features against a phenotype performance in the metadata.
#' @param mmo The mmo object with feature data and metadata
#' @param phenotype The name of the phenotype performance in the metadata
#' @param groups A vector of group names from the metadata containing performance data
#' @param DAM.list A list of DAMs to tag features
#' @param comparisons A list of pairwise comparisons to add fold change columns
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @return A data frame containing regression results for each feature, including effect size, p-value, and fold change columns for specified comparisons.
#'
GetPerformanceFeatureRegression <- function(mmo, phenotype, groups, DAM.list, comparisons, normalization = 'None'){
  feature <- GetNormFeature(mmo, normalization)
  metadata <- mmo$metadata

  # phenotype.sample <- metadata |> filter(group %in% groups) |> pull(sample)
  # phenotype.area <- feature |> select(id, feature, all_of(phenotype.sample))

  performance.linreg <- data.frame(pval = double(), effect.size = double())
  phenotype.df <- data.frame(sample = metadata$sample, group = metadata$group, performance = metadata[,phenotype]) |> filter(.data$group %in% groups)

  regression_results <- data.frame(feature = character(), effect.size = numeric(), p_value = numeric(), is.Spec = logical(), stringsAsFactors = FALSE)
  for (i in 1:nrow(feature)) {
    feature_name <- feature$feature[i]
    feature_df <- data.frame(sample = colnames(feature[,-(1:2)]), feature_value = as.numeric(feature[i, -(1:2)]))
    combined_df <- merge(phenotype.df, feature_df, by='sample')

    fit <- lm(combined_df$performance ~ combined_df$feature_value)
    effect.size <- coef(fit)[2]
    p_value <- summary(fit)$coefficients[2, 4]
    tag <- "else"
    for (list_name in names(DAM.list)) {
      if (feature_name %in% DAM.list[[list_name]]) {
        tag <- list_name
      }
    }
    #is.Spec <- feature_name %in% target

    regression_results <- rbind(regression_results, data.frame(
      feature = feature_name, effect.size = effect.size, p_value = p_value, tag = tag
    ))
  }
  # Add FC columns to regression_results
  for (comparison in comparisons) {
    fc_column <- paste(comparison, "log2FC", sep = "_")
    regression_results[[fc_column]] <- mmo$pairwise[[fc_column]][match(regression_results$feature, mmo$pairwise$feature)]
  }
  return(regression_results)
}

#' GetPerformanceFeatureLMM
#'
#' This function performs linear mixed model analysis of all features against a phenotype performance in the metadata.
#'
#' @param mmo The mmo object with feature data and metadata
#' @param phenotype The name of the phenotype performance in the metadata
#' @param groups A vector of group names from the metadata containing performance data
#' @param DAM.list A list of DAMs to tag features
#' @param comparisons A list of pairwise comparisons to add fold change columns
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'Z')
#' @return A data frame containing regression results for each feature, including effect size, p-value, and fold change columns for specified comparisons.
GetPerformanceFeatureLMM <- function(mmo, phenotype, groups, DAM.list, comparisons, normalization = 'Z'){
  feature <- GetNormFeature(mmo, normalization)
  # if (normalization == 'None'){
  #   feature <- mmo$feature_data
  # } else if (normalization == 'Log'){
  #   feature <- mmo$log
  # } else if (normalization == 'Meancentered'){
  #   feature <- mmo$meancentered
  # } else if (normalization == 'Z'){
  #   feature <- mmo$zscore
  # }
  metadata <- mmo$metadata

  # get phenotype performance data
  phenotype.df <- data.frame(sample = metadata$sample, group = metadata$group, performance = metadata[,phenotype]) |> filter(.data$group %in% groups)
  #create an empty dataframe to store regression results
  regression_results <- data.frame(feature = character(), effect.size = numeric(), p_value = numeric(), is.Spec = logical(), stringsAsFactors = FALSE)
  # iterate regression analysis
  for (i in 1:nrow(feature)) {
    # for each feature, generate phenotype performance X feature value data
    feature_name <- feature$feature[i]
    feature_df <- data.frame(sample = colnames(feature[,-(1:2)]), feature_value = as.numeric(feature[i, -(1:2)]))
    combined_df <- merge(phenotype.df, feature_df, by='sample')

    # linear mixed model
    lmm_fit <- lme4::lmer(performance ~ feature_value + (1|group), data = combined_df)
    fixed_effects <- lme4::fixef(lmm_fit)
    effect.size <- fixed_effects[2]
    p_value <- summary(lmm_fit)$coefficients[2, 5]

    #tag using DAM.list
    tag <- "else"
    for (list_name in names(DAM.list)) {
      if (feature_name %in% DAM.list[[list_name]]) {
        tag <- list_name
      }
    }

    regression_results <- rbind(regression_results, data.frame(
      feature = feature_name, effect.size = effect.size, p_value = p_value, tag = tag
    ))
  }
  # Add FC columns to regression_results
  for (comparison in comparisons) {
    fc_column <- paste(comparison, "log2FC", sep = "_")
    regression_results[[fc_column]] <- mmo$pairwise[[fc_column]][match(regression_results$feature, mmo$pairwise$feature)]
  }
  return(regression_results)
}

#' GetPerformanceFeatureCorrelation
#'
#'
#' This function calculates the Pearson correlation between each feature and a specified phenotype performance in the metadata.
#' @param mmo The mmo object with feature data and metadata
#' @param phenotype The name of the phenotype performance in the metadata
#' @param groups A vector of group names from the metadata containing performance data
#' @param DAM.list A list of DAMs to tag features
#' @param comparisons A list of pairwise comparisons to add fold change columns
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @return A data frame containing correlation results for each feature, including effect size (correlation coefficient), p-value, and fold change columns for specified comparisons.
GetPerformanceFeatureCorrelation <- function(mmo, phenotype, groups, DAM.list, comparisons, normalization = 'None'){
  feature <- GetNormFeature(mmo, normalization)
  metadata <- mmo$metadata

  # phenotype.sample <- metadata |> filter(group %in% groups) |> pull(sample)
  # phenotype.area <- feature |> select(id, feature, all_of(phenotype.sample))

  performance.linreg <- data.frame(pval = double(), effect.size = double())
  phenotype.df <- data.frame(sample = metadata$sample, group = metadata$group, performance = metadata[,phenotype]) |> filter(.data$group %in% groups)

  regression_results <- data.frame(feature = character(), effect.size = numeric(), p_value = numeric(), is.Spec = logical(), stringsAsFactors = FALSE)
  for (i in 1:nrow(feature)) {
    feature_name <- feature$feature[i]
    feature_df <- data.frame(sample = colnames(feature[,-(1:2)]), feature_value = as.numeric(feature[i, -(1:2)]))
    combined_df <- merge(phenotype.df, feature_df, by='sample')
    pearson <- cor.test(combined_df$performance, combined_df$feature_value)
    pval <- pearson[[3]]
    cor <- pearson[[4]]
    tag <- "else"
    for (list_name in names(DAM.list)) {
      if (feature_name %in% DAM.list[[list_name]]) {
        tag <- list_name
      }
    }
    #is.Spec <- feature_name %in% target

    regression_results <- rbind(regression_results, data.frame(
      feature = feature_name, effect.size = cor, p_value = pval, tag = tag
    ))
  }
  # Add FC columns to regression_results
  for (comparison in comparisons) {
    fc_column <- paste(comparison, "log2FC", sep = "_")
    regression_results[[fc_column]] <- mmo$pairwise[[fc_column]][match(regression_results$feature, mmo$pairwise$feature)]
  }
  return(regression_results)
}


#' PlotFoldchangeResistanceRegression
#'
#' This function plots the regression results of a feature against a fold change in resistance, including
#' regression line, p-value, and R-squared value.
#'
#' @param performance_regression The regression results data frame containing effect size, fold change, and tag. The output from GetPerformanceFeatureRegression, GetPerformanceFeatureLMM, or GetPerformanceFeatureCorrelation.
#' @param fold_change The name of the fold change column in the performance_regression dataframe
#' @param color A vector of colors for the points in the plot
#' @param output_dir The output file path for the regression plot
PlotFoldchangeResistanceRegression <- function(performance_regression, fold_change, color, output_dir){
  ind_fit <- lm(data = performance_regression, formula = as.formula(paste("-effect.size ~", fold_change)))
  summary_fit <- summary(ind_fit)
  p_value <- summary_fit$coefficients[2, 4]
  r_squared <- summary_fit$r.squared

  ggplot(performance_regression, aes(x = !!rlang::sym(fold_change), y = -.data$effect.size)) +
    geom_point(size = 0.5, aes(color = .data$tag)) +
    geom_smooth(method = "lm", se = TRUE, color = "black", level = 0.95) +
    xlab(fold_change) +
    ylab('-effect.size') +
    scale_color_manual(values = color) +
    theme_classic() +
    annotate("text", x = Inf, y = Inf, label = paste("p-value:", round(p_value, 500), "\nR-squared:", round(r_squared, 4)),
            hjust = 1.1, vjust = 1.1, size = 3, color = "black")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")
  ggsave(output_dir, height = 6, width = 6)
}

#' PlotFoldchangeResistanceRegression_t
#'
#' This function plots the regression results of a feature against a fold change in resistance, including
#' regression line, p-value, and R-squared value. Transposed version of PlotFoldchangeResistanceRegression.
#'
#' @param performance_regression The regression results data frame containing effect size, fold change, and tag. The output from GetPerformanceFeatureRegression, GetPerformanceFeatureLMM, or GetPerformanceFeatureCorrelation.
#' @param fold_change The name of the fold change column in the performance_regression dataframe
#' @param color A vector of colors for the points in the plot
#' @param output_dir The output file path for the regression plot
PlotFoldchangeResistanceRegression_t <- function(performance_regression, fold_change, color, output_dir){
  ind_fit <- lm(data = performance_regression, formula = as.formula(paste(fold_change, "~ -effect.size")))
  summary_fit <- summary(ind_fit)
  p_value <- summary_fit$coefficients[4]
  r_squared <- summary_fit$r.squared

  ggplot(performance_regression, aes(x = -.data$effect.size, y = !!rlang::sym(fold_change))) +
    geom_point(size = 0.5, aes(color = .data$tag)) +
    geom_smooth(method = "lm", se = TRUE, color = "black", level = 0.95) +
    xlab('-effect.size') +
    ylab(fold_change) +
    scale_color_manual(values = color) +
    theme_classic() +
    annotate("text", x = Inf, y = Inf, label = paste("p-value:", round(p_value, 500), "\nR-squared:", round(r_squared, 4)),
            hjust = 1.1, vjust = 1.1, size = 3, color = "black")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")
  ggsave(output_dir, height = 6, width = 6)
}


get_quadrant <- function(x, y) {
  if (x > 0 & y > 0) return("Q1")
  if (x < 0 & y > 0) return("Q2")
  if (x < 0 & y < 0) return("Q3")
  if (x > 0 & y < 0) return("Q4")
  return("Edge")  # For points on axes
}

#' PlotFoldchangeResistanceQuad
#'
#' This function plots the fold change resistance in a quadrant plot, categorizing points into quadrants based on their effect size and fold change.
#' It also performs a binomial test to assess the distribution of points across quadrants.
#' @param performance_regression The regression results data frame containing effect size, fold change, and tag. The output from GetPerformanceFeatureRegression, GetPerformanceFeatureLMM, or GetPerformanceFeatureCorrelation.
#' @param fold_change The name of the fold change column in the performance_regression dataframe
#' @param color A vector of colors for the points in the plot
#' @param output_dir The output file path for the quadrant plot
PlotFoldchangeResistanceQuad <- function(performance_regression, fold_change, color, output_dir){
  performance_regression <- performance_regression |>
  mutate(
    quadrant = dplyr::case_when(
      -effect.size > 0 & !!rlang::sym(fold_change) > 0 ~ "Q1",
      -effect.size < 0 & !!rlang::sym(fold_change) < 0 ~ "Q3",
      -effect.size < 0 & !!rlang::sym(fold_change) > 0 ~ "Q2",
      -effect.size > 0 & !!rlang::sym(fold_change) < 0 ~ "Q4",
      TRUE ~ "Edge"             # For points on axes
    )
  )
  q_counts <- table(performance_regression$quadrant)
  q13 <- sum(q_counts[c("Q1", "Q3")], na.rm = TRUE)
  q24 <- sum(q_counts[c("Q2", "Q4")], na.rm = TRUE)
  binom_test <- binom.test(q13, q13+q24, p = 0.5, alternative = "two.sided")

  ggplot(performance_regression, aes(x = -.data$effect.size, y = !!sym(fold_change))) +
    geom_point(size = 0.5, aes(color = .data$tag)) +
    xlab('-effect.size') +
    ylab(fold_change) +
    scale_color_manual(values = color) +
    theme_classic() +
    annotate("text", x = Inf, y = Inf, label = paste("p-value:", round(binom_test[[3]], 500)),
            hjust = 1.1, vjust = 1.1, size = 3, color = "black")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")
  ggsave(output_dir, height = 6, width = 6)
}
################### Singlevariate analyses ###################

#' AnovaBarPlot
#'
#' This function generates bar plots for a specified feature across different groups in the metadata, performing ANOVA and Tukey's HSD test for post-hoc analysis.
#'
#' @param mmo The mmo object containing metadata and feature data
#' @param ID_list A list of feature IDs to analyze
#' @param outdir The output directory to save the bar plots and ANOVA results
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param filter_group A boolean indicating whether to filter the feature values by a specific group list (default: FALSE)
#' @param group_list A list of groups to filter the feature values by, if filter_group is TRUE (default: NULL)
AnovaBarPlot <- function(mmo, ID_list, outdir, normalization = 'None', filter_group = FALSE, group_list = NULL) {
  .require_pkg("ggbeeswarm")
  # Extract metadata and feature data
  metadata <- mmo$metadata
  feature_data <- GetNormFeature(mmo, normalization)

  # Iterate through each feature ID
  for (target_id in ID_list) {
    # Extract feature values and merge with metadata
    feature_values <- feature_data |>
      filter(.data$id == target_id) |>
      select(-.data$id, -.data$feature) |>
      t() |>
      as.data.frame()
    colnames(feature_values) <- "value"
    feature_values$sample <- rownames(feature_values)
    feature_values <- merge(feature_values, metadata, by = "sample")
    if (filter_group == TRUE){
      feature_values <- feature_values |> filter(.data$group %in% group_list)
    }
    # Perform ANOVA
    anova <- anova_tukey_dunnett(feature_values, 'value ~ group')
    write_anova(anova, outdir = paste0(outdir,'/', target_id, '_anova.csv'), way = 'oneway')


    # Generate bar plot
    p <- ggplot(feature_values, aes(x = .data$group, y = .data$value, fill = .data$group)) +
      geom_bar(stat = "summary", fun = "mean", position = "dodge") +
      geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.9), width = 0.2) +
      ggbeeswarm::geom_beeswarm() +
      theme_classic() +
      labs(title = paste("Feature:", target_id), x = "Group", y = "Value") +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

    # Save the plot
    ggsave(file.path(outdir, paste0(target_id, "_barplot.png")), plot = p, width = 6, height = 4)
  }
}

#' ExportFeaturesToCSV
#'
#' This function exports selected features, their annotations, and pairwise comparisons to a CSV file.
#'
#' @param mmo The mmo object containing feature data, annotations, and pairwise comparisons
#' @param feature_list A list of feature names to filter and export
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param output_dir The output directory to save the CSV file
ExportFeaturesToCSV <- function(mmo, feature_list, normalization = 'None', output_dir){
  feature <- GetNormFeature(mmo, normalization = normalization) # Get normalized feature data
  # Filter the feature data, annotation, and DA analysis for the list provided
  selected_annotations <- mmo$sirius_annot |> filter(feature %in% feature_list) |>
    select(id = 1, feature = 2, formula = 12, compound = 15, pubchem = 18, ionmass = 21, NPC_pathway = 30, NPC_superclass = 32, NPC_class = 34, ClassyFire_superclass = 36, ClassyFire_class = 38, ClassyFire_subclass = 40, ClassyFire_level5 = 42, ClassyFire_most_specific = 44)
  selected_feature <- feature |> filter(feature %in% feature_list)
  selected_pairwise <- mmo$pairwise |> filter(feature %in% feature_list)
  # Merge all
  merged_df <- merge(selected_annotations, selected_feature, by = 'feature')
  merged_df <- merge(merged_df, selected_pairwise, by = 'feature')

  write.csv(merged_df, output_dir)
}

##### Chemical Diversity indices #####
# Function to calculate unweighted Hill numbers

#' GetFunctionalHillNumber
#'
#' This function calculates the functional Hill number for a given mmo object, normalization method, and distance metric.
#'
#' @param mmo The mmo object containing feature data and metadata
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param q The order of the Hill number to calculate (default: 1)
#' @param distance The distance metric to use for calculating dissimilarity. Options are 'dreams', 'm2ds', or 'cosine' (default: 'dreams')
#' @param filter_feature A boolean indicating whether to filter the feature data by a specific list (default: FALSE)
#' @param feature_list A list of feature names to filter the feature data by, if filter_feature is TRUE (default: NULL)
#' @return A data frame containing the functional Hill number for each group in the metadata, with columns for group and hill number.
GetFunctionalHillNumber <- function(mmo, normalization = 'None',q = 1, distance = 'dreams', filter_feature = FALSE, feature_list = NULL){
  feature <- GetNormFeature(mmo, normalization = normalization)
  metadata <- mmo$metadata
  distance_matrix <- GetDistanceMat(mmo, distance = distance)
  # Scale the  distance matrix to be between 0 and 1

  if (filter_feature == TRUE){
    id_list <- FeatureToID(mmo, feature_list)
    distance_matrix <- distance_matrix[id_list, id_list]
  }
  scaled_dissimilarity <- distance_matrix / max(distance_matrix)
  # Calculate the relative proportions of each feature and reorder them to match the order of the distance matrix
  q.feature <- feature |> filter(.data$id %in% colnames(scaled_dissimilarity))
  relative_proportions <- apply(q.feature[, -(1:2)], 2, function(x) x / sum(x))
  rownames(relative_proportions) <- q.feature$id
  relative_proportions <- relative_proportions[rownames(scaled_dissimilarity), ]
  scaled_dissimilarity <- as.matrix(scaled_dissimilarity)
  raoQ <- colSums(relative_proportions * (scaled_dissimilarity %*% relative_proportions))
  # Calculate Hill
  functional_hill_number <- c()
  if (q == 1){
    mask <- relative_proportions > 0
    Plog <- ifelse(mask, relative_proportions/raoQ * log(relative_proportions/raoQ), 0)
    DP <- scaled_dissimilarity %*% relative_proportions
    vals <- 2 * colSums(Plog * DP)
    functional_hill_number <- exp(-vals)
  } else {
    Pq <- (relative_proportions/raoQ)^q
    DPq <- scaled_dissimilarity %*% Pq
    vals <- colSums(Pq*DPq)
    functional_hill_number <- vals^(1/(1-q))
  }
  names(functional_hill_number) <- colnames(relative_proportions)
  # Get the group information
  groups <- c()
  for (col in colnames(feature)[-c(1, 2)]) {
    groups <- append(groups, metadata[metadata$sample == col, ]$group)
  }

  hill_df <- data.frame(group = groups, hill_number = functional_hill_number)
  return(hill_df)
}

#' GetHillNumbers
#'
#' This function calculates the Hill numbers for a given mmo object, normalization method, and order of the Hill number.
#'
#' @param mmo The mmo object containing feature data and metadata
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param q The order of the Hill number to calculate (default: 0)
#' @param filter_feature A boolean indicating whether to filter the feature data by a specific list (default: FALSE)
#' @param feature_list A list of feature names to filter the feature data by, if filter_feature is TRUE (default: NULL)
#' @return A data frame containing the Hill number for each group in the metadata, with columns for group and hill number.
GetHillNumbers <- function(mmo, normalization = 'None', q = 0, filter_feature = FALSE, feature_list = NULL) {
  feature <- GetNormFeature(mmo, normalization = normalization)
  if (filter_feature == TRUE) {
    feature <- feature |> filter(feature %in% feature_list)
  }
  metadata <- mmo$metadata

  hill_numbers <- apply(feature[, -(1:2)], 2, function(x) {
    p <- x / sum(x)
    if (q == 0) {
      return(length(p))
    } else if (q == 1) {
      return(exp(-sum(p * log(p))))
    } else {
      return((sum(p^q))^(1 / (1 - q)))
    }
  })

  groups <- c()
  for (col in colnames(feature)[-c(1, 2)]) {
    groups <- append(groups, metadata[metadata$sample == col, ]$group)
  }

  hill_df <- data.frame(group = groups, hill_number = hill_numbers)


  return(hill_df)
}

#' GetAlphaDiversity
#'
#' This function calculates the alpha diversity for a given mmo object, order of the Hill number, normalization method, mode (weighted or unweighted), distance metric, and optional feature filtering.
#' @param mmo The mmo object containing feature data and metadata
#' @param q The order of the Hill number to calculate (default: 1)
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param mode The mode of diversity calculation. Options are 'weighted' or 'unweighted' for chemical distance(default: 'weighted')
#' @param distance The distance metric to use for calculating dissimilarity. Options are 'dreams', 'm2ds', or 'cosine' (default: 'dreams')
#' @param filter_feature A boolean indicating whether to filter the feature data by a specific list (default: FALSE)
#' @param feature_list A list of feature names to filter the feature data by, if filter_feature is TRUE (default: NULL)
#' @return A data frame containing the alpha diversity for each group in the metadata, with columns for group and alpha diversity value.
GetAlphaDiversity <- function(mmo, q = 1, normalization = 'None', mode = 'weighted', distance = 'dreams', filter_feature = FALSE, feature_list = NULL){
  if (mode == 'weighted'){
    GetFunctionalHillNumber(mmo, normalization = normalization, q = q, distance = distance, filter_feature = filter_feature, feature_list = feature_list)
  } else if (mode == 'unweighted'){
    GetHillNumbers(mmo, normalization = normalization, q = q, filter_feature = filter_feature, feature_list = feature_list)
  } else{
    print('mode should be weighted or unweighted')
  }
}

#' GetSpecializationIndex
#'
#' This function calculates the specialization index for a given mmo object, normalization method, and optional filtering by groups and features.
#'
#' @param mmo The mmo object containing feature data and metadata
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param filter_group A boolean indicating whether to filter the feature data by a specific group list (default: FALSE)
#' @param group_list A list of groups to filter the feature data by, if filter_group is TRUE (default: NULL)
#' @param filter_feature A boolean indicating whether to filter the feature data by a specific list (default: FALSE)
#' @param feature_list A list of feature names to filter the feature data by, if filter_feature is TRUE (default: NULL)
#' @return A data frame containing the specialization index for each group in the metadata, with columns for group and specialization index.
GetSpecializationIndex <- function(mmo, normalization = 'None', filter_group = FALSE, group_list = NULL, filter_feature = FALSE, feature_list = NULL){
  metadata <- mmo$metadata
  feature <- GetNormFeature(mmo, normalization)

  # All feature or filtered feature
  if (filter_feature == TRUE){
    feature <- feature %>% filter(feature %in% filter_list)
  }
  if (filter_group == TRUE){
    samples <- c()
    for (group in group_list){
      samples <- append(samples, metadata %>% filter(group == !!group) %>% pull(sample))
    }
    feature <- feature %>% select(id, feature, all_of(samples))
  }

  # Get frequency matrix
  Pij <- feature[, -(1:2)]
  rownames(Pij) <- feature$feature
  Pij <- t(Pij)
  Pij <- Pij / rowSums(Pij)
  Pij[is.na(Pij)] <- 0

  # Get shannon diversity index
  # Pij.diversity <- Pij * log(Pij, base = 2)
  # Pij.diversity[is.na(Pij.diversity)] <- 0
  # Pij.diversity <- -rowSums(Pij.diversity)

  # Get specialization index
  Pi <- colSums(Pij)/nrow(Pij) # average frequency of each metabolite
  
  Si <- t((t(Pij)/Pi) * log((t(Pij)/Pi), base = 2)) # specialized degree
  Si[is.na(Si)] = 0
  Si <- colSums(Si)/ncol(Pij) # for each feature get average specialization index
  Pij.specialization <- colSums(t(Pij)*Si) # for each sample get specialization index
  # Retrieve group for each sample from metadata
  groups <- sapply(colnames(feature)[-c(1,2)], function(sample) {
    metadata$group[metadata$sample == sample]
  })
  output <- data.frame(group = groups, specialization = Pij.specialization)
}

#' GetBetaDiversity
#'
#' This function calculates the beta diversity for a given mmo object, method (Generalized Unifrac, bray, jaccard, or CSCS), normalization method, distance metric, and optional feature filtering.
#' Then it returns a distance matrix of beta diversity values.
#'
#' @param mmo The mmo object containing feature data and metadata
#' @param method The method of beta diversity calculation. Options are 'Gen.Uni' for Generalized UniFrac, 'bray' for Bray-Curtis, 'jaccard' for Jaccard, or 'CSCS' for Compound Similarity and Chemical structural and compositional similarity (default: 'Gen.Uni')
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param distance The distance metric to use for calculating dissimilarity. Options are 'dreams', 'm2ds', or 'cosine' (default: 'dreams')
#' @param filter_feature A boolean indicating whether to filter the feature data by a specific list (default: FALSE)
#' @param feature_list A list of feature names to filter the feature data by, if filter_feature is TRUE (default: NULL)
#' @return A distance matrix of beta diversity values
GetBetaDiversity <- function(mmo, method = 'Gen.Uni', normalization = 'None', distance = 'dreams', filter_feature = FALSE, feature_list = NULL) {
  # Get compound distance and build tree for UniFrac
  scaled_dissimilarity <- GetDistanceMat(mmo, distance = distance) / max(GetDistanceMat(mmo, distance = distance))
  if (filter_feature == TRUE) {
    id_list <- FeatureToID(mmo, feature_list)
    scaled_dissimilarity <- scaled_dissimilarity[id_list, id_list]
  }
  compound_tree <- ape::as.phylo(hclust(as.dist(scaled_dissimilarity), method = "average"))

  # Get feature matrix of relative proportion
  feature <- GetNormFeature(mmo, normalization)
  if (filter_group == TRUE){
    samples <- c()
    for (group in group_list){
      samples <- append(samples, metadata |> filter(.data$group == !!group) |> pull(sample))
    }
    feature <- feature |> dplyr::select(.data$id, .data$feature, all_of(samples))
  }
  metadata <- mmo$metadata
  feature <- feature |> filter(.data$id %in% colnames(scaled_dissimilarity))
  relative_proportions <- apply(feature[, -(1:2)], 2, function(x) x / sum(x))
  rownames(relative_proportions) <- feature$id
  relative_proportions <- relative_proportions[rownames(scaled_dissimilarity), ] #reorder
  relative_proportions <- t(relative_proportions)
  # Calculate Generalized UniFrac
  if (method == 'Gen.Uni') {
    guni <- GUniFrac::GUniFrac(relative_proportions, compound_tree, alpha = c(0, 0.5, 1), verbose = TRUE)
    beta_div <- guni$unifracs
  } else if (method == 'bray') {
    beta_div <- as.matrix(vegan::vegdist(relative_proportions, method = 'bray'))
  } else if (method == 'jaccard') {
    beta_div <- as.matrix(vegan::vegdist(relative_proportions, method = 'jaccard'))
  } else if (method == 'CSCS') {
    CSS <- 1-GetDistanceMat(mmo, distance = distance)
    diag(CSS)
    q.feature <- GetNormFeature(mmo, normalization = normalization) |> filter(.data$id %in% colnames(CSS))
    relative_proportions <- apply(q.feature[, -(1:2)], 2, function(x) x / sum(x))
    CSCS_all <- t(relative_proportions) %*% CSS %*% relative_proportions

    sample_names <- colnames(relative_proportions)
    n_samples <- length(sample_names)
    CSCS_matrix <- matrix(NA, nrow = n_samples, ncol = n_samples)
    rownames(CSCS_matrix) <- sample_names
    colnames(CSCS_matrix) <- sample_names
    for (i in 1:n_samples) {
      for (j in 1:n_samples) {
        CSCS_matrix[i,j] <- CSCS_all[i, j] / max(CSCS_all[i, i], CSCS_all[j, j])
      }
    }
    beta_div <- 1-CSCS_matrix
  } else {
    stop("Invalid method. Please use 'Gen.Uni', 'bray' or 'jaccard'")
  }


  return(beta_div)
}




#' CalculateGroupBetaDistance
#'
#' This function calculates the beta diversity distance between a reference group and other groups in the metadata.
#'
#' @param mmo The mmo object containing feature data and metadata
#' @param beta_div The beta diversity distance matrix
#' @param reference_group The name of the reference group to compare against
#' @param groups A vector of group names from the metadata to calculate distances for
#' @return A data frame containing the group names, sample names, and their corresponding beta diversity distances from the reference group.
#'
CalculateGroupBetaDistance <- function(mmo, beta_div, reference_group, groups) {
  metadata <- mmo$metadata
  distances <- data.frame(group = character(), distance = numeric())

  for (group in groups) {
    if (group != reference_group) {
      group_samples <- metadata |> filter(group == !!group) |> pull(sample)
      reference_samples <- metadata |> filter(group == !!reference_group) |> pull(sample)

      for (sample in group_samples) {
        for (ref_sample in reference_samples) {
          distance <- beta_div[sample, ref_sample]
          distances <- rbind(distances, data.frame(group = group,sample = sample, distance = distance))
        }
      }
    }
  }

  return(distances)
}
