########################################################################################
# Define functions for data import and normalization
########################################################################################
#' Import mzmine feature data and metadata to create a mmo object
#' @description
#' This function reads mzmine feature data and metadata from specified directories,
#' to initiate a mmo object containing feature data and metadata
#'
#' Import mzmine feature data and metadata to create a mmo object
#'
#' @param mzmine_dir Path to the mzmine feature data CSV file
#' @param metadata_dir Path to the metadata CSV file (must include sample_col and group_col)
#' @param group_col Column name in the metadata file used for grouping samples together i.e into treatments or species.
#' @param sample_col Column in metadata file used to identify and match individual samples 
#' @param mz_col Optional m/z column name (defaults to "mz" or "row m/z")
#' @param rt_col Optional RT column name (defaults to "rt" or "row retention time")
#' @return A mmo object containing the feature data and metadata
#' @export
GetMZmineFeature <- function(mzmine_dir, metadata_dir, group_col, sample_col,
                             mz_col = NULL, rt_col = NULL) {
  mmo <- list()
  data <- read.csv(mzmine_dir, check.names = FALSE,stringsAsFactors = FALSE, na.strings = c("", "NA"))

  metadata <- read.csv(metadata_dir, check.names = FALSE)

  if (missing(group_col) || !(group_col %in% names(metadata)))
    stop("group_col must be provided and exist in the metadata file.")
  if (missing(sample_col) || !(sample_col %in% names(metadata)))
    stop("sample_col must be provided and exist in the metadata file.")

  # --- detect mz / rt (or use overrides) ---
  if (is.null(mz_col))
    mz_col <- if ("mz" %in% names(data)) "mz" else if ("row m/z" %in% names(data)) "row m/z" else stop("No m/z column found.")
  if (is.null(rt_col))
    rt_col <- if ("rt" %in% names(data)) "rt" else if ("row retention time" %in% names(data)) "row retention time" else stop("No RT column found.")

  # --- feature column ---
  data <- data |> mutate(feature = paste(.data[[mz_col]], .data[[rt_col]], sep = "_"))

  # --- area columns from EXACT metadata filenames ---
  samples <- trimws(metadata[[sample_col]])                     # e.g., "I1_C1_1.mzML"
  expected_old <- paste0("datafile:", samples, ":area")         # old export column names
  expected_new <- paste0(samples, " Peak area")                 # new export column names

  area_columns <- intersect(expected_old, names(data))
  if (length(area_columns) == 0)
    area_columns <- intersect(expected_new, names(data))
  if (length(area_columns) == 0)
    stop("No area columns matched the EXACT names from metadata[[sample_col]].")

  #make sure areas are treated numerically
  data[area_columns] <- lapply(data[area_columns], function(x) as.numeric(gsub(",", "", x)))


  # --- build feature_df (keep your original layout) ---
  if (!("id" %in% names(data))) data$id <- seq_len(nrow(data))  # minimal safety
  feature_df <- data |> select(.data$id, .data$feature, all_of(area_columns))
  feature_df$id <- gsub(" ", "", feature_df$id)

  # --- finalize metadata and output ---
  metadata$sample <- paste0("datafile.", samples, ".area")
  mmo$feature_data <- feature_df
  mmo$metadata <- metadata
  mmo$pairwise <- data.frame(feature = mmo$feature_data$feature, id = mmo$feature_data$id)
  mmo$metadata$group <- as.factor(mmo$metadata[[group_col]])

  print("MMO object created.")
  print(paste0("Feature number: ", nrow(mmo$feature_data)))
  print(paste0(nrow(mmo$metadata), " samples in ", length(unique(mmo$metadata$group)), " groups"))
  return(mmo)
}

#' Add feature_info to an existing mmo from a full_feature CSV
#'
#' @param mmo An existing mmo list (will be returned with $feature_info added)
#' @param full_feature_csv Path to the MZmine full feature table (CSV)
#' @return The same mmo with mmo$feature_info set
#' @export
AddFeatureInfo <- function(mmo, full_feature_csv) {
  required <- c(
    "id","rt","rt_range:min","rt_range:max","mz","mz_range:min","mz_range:max",
    "feature_group","ion_identities:iin_id","ion_identities:ion_identities"
  )

  ff <- read.csv(full_feature_csv, check.names = FALSE, stringsAsFactors = FALSE)

  missing <- setdiff(required, names(ff))
  if (length(missing)) {
    stop("Missing required columns in full_feature_csv: ",
         paste(missing, collapse = ", "))
  }

  # Keep only the required columns, in order
  mmo$feature_info <- ff[required]
  rownames(mmo$feature_info) <- NULL
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
#' @export
#' @examplesIf FALSE
#' mmo <- SwitchGroup(mmo, new_group_col = "genotype")
SwitchGroup <- function(mmo, new_group_col) {
  if (missing(new_group_col) || !(new_group_col %in% colnames(mmo$metadata))) {
    stop("new_group_col must be provided and must exist in the metadata file.")
  }
  mmo$metadata$group <- as.factor(mmo$metadata[[new_group_col]])
  print(paste0('Group column switched to ', new_group_col))
  print(paste0(length(unique(mmo$metadata$group)), ' groups in total'))
  print(paste0('The list of groups are: ', paste(levels(mmo$metadata$group), collapse = ', ')))
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
#' @export
#' @examplesIf FALSE
#' mmo <- AddSiriusAnnot(mmo, 
#'  canopus_structuredir = "path/to/structure_identification.tsv", 
#'  canopus_formuladir = "path/to/canopus_formula_summary.tsv"
#' )
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
#' @examplesIf FALSE
#' mmo <- AddCustomAnnot(mmo, DB_file = "path/to/custom_db.csv", mztol = 5, rttol = 0.5)
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
#' Run this function before MassNormalization(), LogNormalization(), MeancenterNormalization(), or ZNormalization().
#'
#' @param mmo The mmo object
#' @param method The method to use for replacement. Options are 'one' (replace with 1) or 'half_mean' (replace with half of the smallest non-zero value in the row)
#' @return The mmo object with replaced values in the feature data (mmo$feature_data)
#' @export
#' @examplesIf FALSE
#' mmo <- ReplaceZero(mmo, method = 'one')
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
#' Use sample mass in the metadata file to normalize the peak area
#'
#' This function normalizes the peak area in the feature data of the mmo object by the mass of each sample, provided in the metadata.
#' The feature data is replaced by (original value * mean mass) / sample mass.
#'
#' @param mmo The mmo object
#' @return The mmo object with normalized feature data (mmo$feature_data)
#' @export
#' @examplesIf FALSE
#' mmo <- MassNormalization(mmo)
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
#' Run ReplaceZero() before this function to avoid -Inf values.
#' 
#' @param mmo The mmo object
#' @return The mmo object with log-normalized feature data (mmo$log)
#' @export
#' @examplesIf FALSE
#' mmo <- LogNormalization(mmo)
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
#' @export
#' @examplesIf FALSE
#' mmo <- MeancenterNormalization(mmo)
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
#' @export
#' @examplesIf FALSE
#' mmo <- ZNormalization(mmo)
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
#' This function reads cosine, DREAMS, and MS2DeepScore molecular networking outputs from MZmine,
#' then transform the similarity to distance and adds the dissimilarity matrices to the mmo object.
#'
#' @param mmo The mmo object
#' @param cos_dir Path to the cosine similarity CSV file from MZMine (molecular networking)
#' @param dreams_dir Path to the DREAMS similarity CSV file from MZMine (molecular networking)
#' @param m2ds_dir Path to the MS2DeepScore similarity CSV file from MZMine (molecular networking)
#' @return The mmo object with dissimilarity matrices added (mmo$cos.dissim, mmo$dreams.dissim, mmo$m2ds.dissim)
#' @export
#' @examplesIf FALSE
#' mmo <- AddChemDist(mmo, 
#'  cos_dir = "path/to/cosine_similarity.csv", 
#'  dreams_dir = "path/to/dreams_similarity.csv", 
#'  m2ds_dir = "path/to/ms2deepscore_similarity.csv"
#' )
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
    cluster_index <- stats::setNames(seq_along(clusters), clusters)
    
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
#' Use this function before plotting heatmaps or other visualizations to ensure consistent group ordering.
#' 
#' @param mmo The mmo object
#' @param group_order A vector specifying the desired order of groups
#' @return The mmo object with reordered samples
#' @export
#' @examplesIf FALSE
#' mmo <- ReorderGroups(mmo, group_order = c("Control", "Treatment1", "Treatment2"))
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

#' Retrieve feature data from the mmo object, with normalization options
#'
#' This function retrieves the feature data from the mmo object based on the specified normalization method.
#'
#' @param mmo The mmo object
#' @param normalization The normalization method to use. Options are 'None', 'Log', 'Meancentered', or 'Z'
#' @return The feature data corresponding to the specified normalization method
#' @export
#' @examplesIf FALSE
#' feature_data <- GetNormFeature(mmo, normalization = 'Log')
GetNormFeature <- function(mmo, normalization){
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
#' @export
#' @examplesIf FALSE
#' distance_matrix <- GetDistanceMat(mmo, distance = 'dreams')
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
#' @export 
#' @examplesIf FALSE
#' feature_ids <- FeatureToID(mmo, feature_names = c("100.0_5.0", "150.0_10.0"))
#' feature_ids <- FeatureToID(mmo, feature_names = mmo$feature_data$feature[1:10])
#' feature_ids <- FeatureToID(mmo, 
#'  feature_names = Glucosinolates
#' ) # if Glucosinolates is a vector of feature names
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
#' @export
#' @examplesIf FALSE
#' feature_names <- IDToFeature(mmo, feature_ids = c("1219", "2250", "3360"))
#' feature_names <- IDToFeature(mmo, feature_ids = mmo$feature_data$id[1:10])
#' feature_names <- IDToFeature(mmo, 
#'  feature_ids = FeatureToID(mmo, feature_names = Glucosinolates)
#' ) # if Glucosinolates is a vector of feature names
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
#' This function calculates and returns a dataframe of mean feature values for each group in the mmo object, with options for normalization and filtering.
#' Use SwitchGroup() to change the grouping variable before running this function.
#'
#' @param mmo The mmo object
#' @param normalization The normalization method to use. Options are 'None', 'Log', 'Meancentered', or 'Z'
#' @param filter_feature Boolean to filter features based on a provided list (default: FALSE)
#' @param feature_list A vector of feature names to filter (default: NULL)
#' @param filter_group Boolean to filter groups based on a provided list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
#' @return A data frame containing the mean feature values for each group
#' @export
#' @examplesIf FALSE
#' group_means <- GetGroupMeans(mmo, normalization = 'Log')
#' group_means <- GetGroupMeans(mmo, 
#'  normalization = 'None', 
#'  filter_feature = TRUE, feature_list = Glucosinolates
#' ) # if Glucosinolates is a vector of feature names
#' group_means <- GetGroupMeans(mmo, 
#'  normalization = 'Z', 
#'  filter_group = TRUE, 
#'  group_list = c("Control", "Treatment1")
#' )
#' group_means <- GetGroupMeans(mmo, 
#'  normalization = 'Meancentered', 
#'  filter_feature = TRUE, feature_list = Glucosinolates, 
#'  filter_group = TRUE, group_list = c("Control", "Treatment1")
#' ) # if Glucosinolates is a vector of feature names   
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
#' This function calculates and returns a dataframe of log2 fold change values for each group compared to a specified control group.
#' Takes inputs from GetGroupMeans() function.
#' 
#' @param group_means A data frame containing the mean feature values for each group
#' @param control_group The name of the control group to compare against
#' @return A data frame with log2 fold change values for each group compared to the control group
#' @export
#' @examplesIf FALSE
#' fold_change <- GetLog2FoldChange(GetGroupMeans(mmo, normalization = 'Log'), control_group = 'Control')
GetLog2FoldChange <- function(group_means, control_group) {
  control_means <- group_means[[control_group]]
  fold_change <- group_means |>
    mutate(across(-.data$id, ~ log2(. / control_means)))

  return(fold_change)
}


#' Perform ANOVA and Tukey's HSD test on the mmo object
#'
#' This function performs ANOVA and Tukey's HSD test on the feature data of the mmo object,
#' Returns a list of ANOVA results, Tukey's HSD results, Tukey's significance letters, and Dunnett's test results.
#'
#' @param df The data frame containing the feature data and metadata
#' @param formula The formula for the ANOVA test, e.g., "feature ~ group"
#' @return A list containing the ANOVA results, Tukey's HSD results, Tukey's significance letters, and Dunnett's test results
#' @export
#' @examplesIf FALSE
#' anova_results <- anova_tukey_dunnett(df = merged_data, formula = "feature_value ~ group")
anova_tukey_dunnett <- function(df, formula) {
  .require_pkg("DescTools")
  .require_pkg("multcompView")
  aov_res <- aov(as.formula(formula), data = df)
  tukey_res <- TukeyHSD(aov_res)
  tukey_sig <- multcompView::multcompLetters4(aov_res, tukey_res)
  dunnett_res <- DescTools::DunnettTest(as.formula(formula), data = df)
  return(list(aov_res = aov_res, tukey_res = tukey_res, tukey_sig = tukey_sig, dunnett_res = dunnett_res))
}

#' Write results of anova_tukey_dunnett to a CSV file
#'
#' This function writes the results of ANOVA and Tukey's HSD test to a CSV file.
#'
#' @param anova_data A list containing the results of ANOVA and Tukey's HSD test
#' @param outdir The output directory where the results will be saved
#' @param way The type of ANOVA test to perform. Options are 'oneway' or 'twoway'
#' @export 
#' @examplesIf FALSE
#' write_anova(anova_data = anova_results, outdir = "anova_tukey_results.csv", way = 'oneway')
#' 
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
#' The function returns the PERMANOVA results, raw pairwise comparison results, and matrices of adjusted p-values, F values, and R square for pairwise comparisons
#' 
#' @param data A data frame or distance matrix for PERMANOVA
#' @param metadata A data frame containing sample metadata, including a 'group' column
#' @param mode The mode of the input data: 'data' for raw data or 'distance' for a distance matrix
#' @param filter_group Boolean to filter groups based on a provided list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
#' @param permutations The number of permutations for the PERMANOVA test (default: 5000)
#' @return A list containing the PERMANOVA results, raw pairwise comparison results, and matrices of adjusted p-values, F values, and R square for pairwise comparisons
#' @export
#' @examplesIf FALSE
#' permanova_results <- permanova_stat(
#'  data = feature_data, metadata = mmo$metadata, 
#'  mode = 'data', filter_group = TRUE, group_list = c("Control", "Treatment1"), 
#'  permutations = 5000
#' )
#' permanova_results <- permanova_stat(
#'  data = betadiv, metadata = mmo$metadata, 
#'  mode = 'distance', permutations = 10000
#' )
permanova_stat <- function(data, metadata, mode, filter_group = FALSE, group_list = NULL, permutations = 5000){ 
  .require_pkg("vegan")
  .require_pkg("stringr")
  .require_pkg("pairwiseAdonis")
  #Perform PERMANOVA
  if (mode == 'data'){
    dist_mat <- stats::dist(data)
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
  pairwise_p_matrix <- matrix(NA, nrow = length(unique(metadata$group)), ncol = length(unique(metadata$group)))
  rownames(pairwise_p_matrix) <- unique(metadata$group)
  colnames(pairwise_p_matrix) <- unique(metadata$group)
  for (i in 2:length(pairwise_permanova)){
    group1 <- stringr::str_split(names(pairwise_permanova)[i], '_vs_')[[1]][1]
    group2 <- stringr::str_split(names(pairwise_permanova)[i], '_vs_')[[1]][2]
    pval <- pairwise_permanova[[i]][1,6]
    pairwise_p_matrix[group2, group1] <- pval
  }
  # Organize posthoc F values for visualization
  pairwise_F_matrix <- matrix(NA, nrow = length(unique(metadata$group)), ncol = length(unique(metadata$group)))
  rownames(pairwise_F_matrix) <- unique(metadata$group)
  colnames(pairwise_F_matrix) <- unique(metadata$group)
  for (i in 2:length(pairwise_permanova)){
    group1 <- stringr::str_split(names(pairwise_permanova)[i], '_vs_')[[1]][1]
    group2 <- stringr::str_split(names(pairwise_permanova)[i], '_vs_')[[1]][2]
    Fval <- pairwise_permanova[[i]][1,4]
    pairwise_F_matrix[group2, group1] <- Fval
  }
  # Organize posthoc R^2 values for visualization
  pairwise_R2_matrix <- matrix(NA, nrow = length(unique(metadata$group)), ncol = length(unique(metadata$group)))
  rownames(pairwise_R2_matrix) <- unique(metadata$group)
  colnames(pairwise_R2_matrix) <- unique(metadata$group)
  for (i in 2:length(pairwise_permanova)){
    group1 <- stringr::str_split(names(pairwise_permanova)[i], '_vs_')[[1]][1]
    group2 <- stringr::str_split(names(pairwise_permanova)[i], '_vs_')[[1]][2]
    R2val <- pairwise_permanova[[i]][1,3]
    pairwise_R2_matrix[group2, group1] <- R2val
  }
  return(list(permanova_res = permanova_res, pairwise_raw = pairwise_permanova, pairwise_p_matrix = pairwise_p_matrix, pairwise_F_matrix = pairwise_F_matrix, pairwise_R2_matrix = pairwise_R2_matrix))
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
#' @param correction The method for multiple comparison correction. Options are 'BH', 'holm', 'bonferroni', etc. Inherits from p.adjust() (default: 'BH')
#' @return The mmo object with pairwise comparison results added to mmo$pairwise
#' @export
#' @examplesIf FALSE
#' mmo <- PairwiseComp(mmo, group1 = 'Control', group2 = 'Treatment1')
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


#' Generates lists of DAMs (Differentially Accumulated Metabolites) for each comparison in the mmo object
#'
#' This function generates lists of upregulated and downregulated DAMs for each pairwise comparison in the mmo object.
#' It uses log2 fold change and adjusted p-value thresholds to determine significance.
#' Make sure to run PairwiseComp() for all desired comparisons before using this function.
#'
#' @param mmo The mmo object with pairwise comparison matrix
#' @param fc_cutoff The threshold of log2 fold change to be considered significant (default: 0.5849625, which is log2(1.5))
#' @param pval_cutoff The threshold of adjusted p-value to be considered significant (default: 0.05)
#' @return A list containing two lists: DAMs_up and DAMs_down
#' @export
#' @examplesIf FALSE
#' dams <- GetDAMs(mmo, fc_cutoff = 0.5849625, pval_cutoff = 0.05)
#' dams_up <- dams$DAMs_up
#' dams_down <- dams$DAMs_down
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



#' Volcano plot for visualizing differential metabolite analysis results
#'
#' This function generates a volcano plot using data from mmo$pairwise (PairwiseComp(mmo, 'group1', 'group2') should be precended),
#' highlighting upregulated and downregulated features based on log2 fold change and adjusted p-value
#'
#' @param mmo The mmo object with pairwise comparison matrix
#' @param comp The comparison to visualize, e.g., 'group1_vs_group2
#' @param topk The number of top features to label in the plot (default: 10)
#' @param pthr The threshold of adjusted p-value to be considered significant (default: 0.05)
#' @param outdir The output file path for the volcano plot (default: 'volcano.png')
#' @param height The height of the output plot in inches (default: 5)
#' @param width The width of the output plot in inches (default: 5)
#' @export 
#' @examplesIf FALSE
#' VolcanoPlot(
#'  mmo, comp = 'Control_vs_Treatment1', 
#'  topk = 10, pthr = 0.05, 
#'  outdir = 'volcano_con_tre1.png', height = 5, width = 5
#' )
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


#' Plots PCA and performs PERMANOVA
#' 
#' This function performs PCA analysis and generates a PCA plot with optional filtering of features and groups.
#' It also conducts PERMANOVA and saves the results to CSV files.
#'
#' @param mmo The mmo object with feature data and metadata
#' @param color A vector of colors for the groups in the plot. Make sure the names correspond to the group names in metadata
#' @param outdir The output file path for the PCA plot (default: 'PCA')
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'Z')
#' @param filter_feature Boolean to filter features by feature_list (default: FALSE)
#' @param feature_list A vector of feature names to filter (default: NULL)
#' @param filter_group Boolean to filter groups by group_list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
#' @param label Boolean to indicate whether to label points with sample names (default: TRUE)
#' @export 
#' @examplesIf FALSE
#' PCAplot(
#'  mmo, color = c("Control" = "blue", "Treatment1" = "red", "Treatment2" = "green"), 
#'  outdir = 'PCA_plot', normalization = 'None', 
#'  filter_feature = FALSE, filter_group = FALSE, label = FALSE
#' )
#' PCAplot(
#'  mmo, color = c("Control" = "blue", "Treatment1" = "red"), 
#'  outdir = 'PCA_plot', normalization = 'Z', 
#'  filter_feature = TRUE, feature_list = Glucosinolates, 
#'  filter_group = TRUE, group_list = c("Control", "Treatment1"), label = TRUE
#' )
PCAplot <- function(mmo, color, outdir = 'PCA', normalization = 'Z', filter_feature = FALSE, feature_list = NULL, filter_group = FALSE, group_list = NULL, label = TRUE){
  .require_pkg("ggrepel")
  .require_pkg("stats")
  metadata <- mmo$metadata
  feature <- GetNormFeature(mmo, normalization)
  if (filter_feature == TRUE){
    feature <- feature |> filter(.data$feature %in% feature_list)
  }
  # Perform PCA on normalized feature data
  feature_data_pca <- feature[, -(1:2)]
  feature_data_pca <- t(feature_data_pca) # samples as rows, features as columns
  pca_res <- stats::prcomp(feature_data_pca, scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x)
  pca_df$group <- metadata$group[match(rownames(pca_df), metadata$sample)]

  if (filter_group == TRUE){
    pca_df <- pca_df |> filter(.data$group %in% group_list)
  }
  if (label == TRUE){
    plot <- ggplot(pca_df, aes(x = .data$PC1, y = .data$PC2, color = .data$group, label = rownames(pca_df))) +
      geom_point(size = 3) +
      ggrepel::geom_label_repel(aes(label = rownames(pca_df)), size = 3) +
      theme_classic() +
      labs(x = "PC1", y = "PC2") +
      scale_color_manual(values = color)+
      stat_ellipse(aes(group = .data$group), level = 0.90)
  } else {
    plot <- ggplot(pca_df, aes(x = .data$PC1, y = .data$PC2, color = .data$group)) +
      geom_point(size = 3) +
      # geom_label_repel(aes(label = rownames(pca_df)), size = 3) +
      theme_classic() +
      labs(x = "PC1", y = "PC2") +
      scale_color_manual(values = color)+
      stat_ellipse(aes(group = .data$group), level = 0.90)
  }
  plot
  ggsave(paste0(outdir, '.pdf'), width = 6, height = 6)

  permanova <- permanova_stat(feature_data_pca, metadata, mode = 'data', filter_group = filter_group, group_list = group_list)
  write.csv(permanova$permanova_res, paste0(outdir, '_permanova_results.csv'))
  write.csv(as.data.frame(permanova$pairwise_raw), paste0(outdir, '_pairwise_permanova_results.csv'))
  write.csv(as.data.frame(permanova$pairwise_p_matrix), paste0(outdir, '_pairwise_permanova_pvalue_matrix.csv'))
  write.csv(as.data.frame(permanova$pairwise_F_matrix), paste0(outdir, '_pairwise_permanova_Fvalue_matrix.csv'))
  write.csv(as.data.frame(permanova$pairwise_R2_matrix), paste0(outdir, '_pairwise_permanova_R2value_matrix.csv'))
}




#' PLS-DA plot with feature loadings
#'
# This function performs PLS-DA analysis and generates a PLS-DA plot with feature loadings.
#'
#' @param mmo The mmo object with feature data and metadata
#' @param color A vector of colors for the groups in the plot. Make sure the names correspond to the group names in metadata
#' @param topk The number of top features to display in the plot (default: 10)
#' @param outdir The output file path for the PLS-DA plot
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'Z')
#' @param filter_feature Boolean to filter features by feature_list (default: FALSE)
#' @param feature_list A vector of feature names to filter (default: NULL)
#' @param filter_group Boolean to filter groups by group_list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
#' @export 
#' @examplesIf FALSE
#' PLSDAplot(
#'  mmo, color = c("Control" = "blue", "Treatment1" = "red", "Treatment2" = "green"), 
#'  topk = 10, outdir = 'PLSDA_plot.pdf', normalization = 'Z', 
#'  filter_feature = FALSE, filter_group = FALSE
#' )
#' PLSDAplot(
#'  mmo, color = c("Control" = "blue", "Treatment1" = "red"), 
#'  topk = 5, outdir = 'PLSDA_plot.pdf', normalization = 'Log', 
#'  filter_feature = TRUE, feature_list = Glucosinolates, 
#'  filter_group = TRUE, group_list = c("Control", "Treatment1")
#' )
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
  scale_color_manual(values = color) +
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

#' Generate input files to be used for pheatmap from the mmo object
#'
#' This function generates heatmap inputs from the mmo object, including fold change or mean values,
#' distance matrix, and row labels for custom-annotated features.
#'
#' @param mmo The mmo object with sirius annotation and normalized data
#' @param filter_feature Boolean to filter features by feature_list (default: FALSE)
#' @param feature_list A vector of feature names to filter (default: NULL)
#' @param filter_group Boolean to filter groups by group_list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
#' @param summarize The summarization method to use. Options are 'fold_change' or 'mean' (default: 'mean')
#' @param control_group The group to use as control for fold change calculation (default: 'ctrl')
#' @param normalization The normalization method to use. Options are 'None', 'Log', 'Meancentered', or 'Z'
#' @param distance The distance metric to use. Options are 'dreams', 'cosine', or 'm2ds' (default: 'dreams')
#' @return A list containing the following elements:
#' - FC_matrix: A matrix of fold change or mean values
#' - dist_matrix: A distance matrix based on the specified distance metric
#' - row_label: A vector of row labels for custom-annotated features (See AddCustomAnnot()). If no custom annotation is available, feature IDs are used.
#' - heatmap_data: A data frame containing the heatmap data with feature IDs and values
#' @export
#' @examplesIf FALSE
#' # Generate heatmap inputs to visualize fold change values with log normalization and dreams distance
#' heatmap_inputs <- GenerateHeatmapInputs(
#'  mmo, summarize = 'fold_change', control_group = 'Control', 
#'  normalization = 'None', distance = 'dreams'
#' )
#' # Generate heatmap inputs to visualize mean values
#' heatmap_inputs <- GenerateHeatmapInputs(
#'  mmo, summarize = 'mean', normalization = 'None', distance = 'dreams'
#' )
#' # The resulting list contains FC_matrix, dist_matrix, row_label, and heatmap_data
#' # A heatmap can be generated using pheatmap
#' # 'clustering_distance_rows' option make the dendrogram follows chemical distances of features. 
#' #  -Delete this option to visualize the heatmap following cannonical clustering
#' pheatmap(mat = heatmap_inputs$FC_matrix, 
#'     cluster_rows = TRUE, #do not change
#'     clustering_distance_rows = heatmap_inputs$dist_matrix, 
#'     cluster_cols = TRUE, 
#'     clustering_method = "average", #UPGMA
#'     show_rownames = TRUE, 
#'     show_colnames = TRUE,
#'     cellwidth = 25,
#'     cellheight = 0.05,
#'     treeheight_row = 100,
#'     fontsize_row = 3,
#'     fontsize_col = 15,
#'     scale = 'none',
#'     annotation_names_row = TRUE,
#'     labels_row = heatmap_inputs$row_label,
#'     )
GenerateHeatmapInputs <- function(mmo, filter_feature = FALSE, feature_list = NULL,
                                filter_group = FALSE, group_list = NULL,
                                summarize = 'mean', control_group = 'ctrl',
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

#' PlotNPCStackedBar
#' 
#' This function generates a stacked bar plot showing the count of features in each group categorized by NPC_pathway.
#' It uses the mmo object with sirius annotation and normalized data.
#' Make sure you don't run ReplaceZero() before using this function, as it may remove presence/absence information.
#' 
#' @param mmo The mmo object with sirius annotation and normalized data
#' @param group_col The column name in metadata to use for grouping samples
#' @param output_file The output file path for the stacked bar plot (e.g., 'NPC_stacked_bar.png')
#' @param width The width of the output plot
#' @param height The height of the output plot
#' @export
#' @examplesIf FALSE
#' PlotNPCStackedBar(
#'  mmo, group_col = 'treatment', 
#'  output_file = 'NPC_stacked_bar.png', width = 6, height = 3
#' )
PlotNPCStackedBar <- function(mmo, group_col, output_file, width = 6, height = 3) {
  mmo <- SwitchGroup(mmo, group_col)
  feature_data <- mmo$feature_data
  metadata <- mmo$metadata

  # For each group, get features present in any sample
  group_features <- lapply(unique(metadata$group), function(grp) {
    samples <- metadata |> filter(.data$group == grp) |> pull(.data$sample)
    present <- feature_data |> select(all_of(samples))
    feature_data$feature[base::rowSums(!is.na(present) & present > 0) > 0]
  })
  names(group_features) <- unique(metadata$group)

  # Build long data frame: feature, group, NPC_pathway
  annot <- mmo$sirius_annot[, c("feature", "NPC#pathway")]
  colnames(annot) <- c("feature", "NPC_pathway")

  bar_df <- do.call(rbind, lapply(names(group_features), function(grp) {
    data.frame(
      feature = group_features[[grp]],
      group = grp,
      stringsAsFactors = FALSE
    )
  }))
  bar_df <- merge(bar_df, annot, by = "feature", all.x = TRUE)
  bar_df <- bar_df[!is.na(bar_df$NPC_pathway) & bar_df$NPC_pathway != "", ]

  # Count features per group and NPC_pathway
  plot_df <- bar_df |>
    dplyr::group_by(.data$group, .data$NPC_pathway) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop")

  # Set colors for NPC_pathway
  npc_pathways <- unique(plot_df$NPC_pathway)
  .require_pkg("RColorBrewer")
  bar_colors <- setNames(RColorBrewer::brewer.pal(min(length(npc_pathways), 8), "Set2"), npc_pathways)

  # Plot stacked bar by NPC_pathway
  ggplot(plot_df, aes(x = .data$group, y = .data$count, fill = .data$NPC_pathway, label = .data$count)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(group = .data$NPC_pathway), position = position_stack(vjust = 0.5), size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(values = bar_colors) +
    coord_flip() +
    labs(x = "Group", y = "Feature Count", fill = "NPC Pathway") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(output_file, width = width, height = height)
}

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
#' @export
#' @examplesIf FALSE
#' # Perform enrichment analysis for a list of features using NPC_pathway level
#' sig_terms <- CanopusLevelEnrichmentAnal(
#'  mmo, list_test = c("feature1", "feature2"), pthr = 0.1, 
#'  sig = TRUE, term_level = 'NPC_pathway', representation = 'greater'
#' )
#' # Perform enrichment analysis for a list of features using ClassyFire_class level and return all terms
#' all_terms <- CanopusLevelEnrichmentAnal(
#'  mmo, list_test = c("feature1", "feature2"), pthr = 0.1, 
#'  sig = FALSE, term_level = 'ClassyFire_class', representation = 'greater'
#' )
CanopusLevelEnrichmentAnal <- function(mmo,list_test, pthr = 0.1, sig=TRUE, term_level = 'NPC_pathway', representation = 'greater'){
  all_feature <- mmo$sirius_annot
  subset_feature <- mmo$sirius_annot |> filter(.data$feature %in% list_test)
  # print(paste('total features:', nrow(all_feature), 'list_test features:', nrow(subset_feature)))
  # Select the appropriate term level for enrichment analysis
  if (term_level == "NPC_pathway") {
    all_feature$classifications_split <- all_feature[['NPC#pathway']]
    subset_feature$classifications_split <- subset_feature[['NPC#pathway']]
  } else if (term_level == "NPC_superclass") {
    all_feature$classifications_split <- all_feature[['NPC#superclass']]
    subset_feature$classifications_split <- subset_feature[['NPC#superclass']]
  } else if (term_level == "NPC_class") {
    all_feature$classifications_split <- all_feature[['NPC#class']]
    subset_feature$classifications_split <- subset_feature[['NPC#class']]
  } else if (term_level == "ClassyFire_superclass") {
    all_feature$classifications_split <- all_feature[['ClassyFire#superclass']]
    subset_feature$classifications_split <- subset_feature[['ClassyFire#superclass']]
  } else if (term_level == "ClassyFire_class") {
    all_feature$classifications_split <- all_feature[['ClassyFire#class']]
    subset_feature$classifications_split <- subset_feature[['ClassyFire#class']]
  } else if (term_level == "ClassyFire_subclass") {
    all_feature$classifications_split <- all_feature[['ClassyFire#subclass']]
    subset_feature$classifications_split <- subset_feature[['ClassyFire#subclass']]
  } else if (term_level == "ClassyFire_level5") {
    all_feature$classifications_split <- all_feature[['ClassyFire#level 5']]
    subset_feature$classifications_split <- subset_feature[['ClassyFire#level 5']]
  } else if (term_level == "ClassyFire_most_specific") {
    all_feature$classifications_split <- all_feature[['ClassyFire#most specific class']]
    subset_feature$classifications_split <- subset_feature[['ClassyFire#most specific class']]
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
    filter(.data$fdr < pthr) |>
    arrange(.data$fdr)
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
#' @export 
#' @examplesIf FALSE
#' CanopusListEnrichmentPlot(
#'  mmo, feature_list = DAMs_up$control_vs_treatment1.up, 
#'  pthr = 0.05, outdir = 'canopus_enrichment_plot.pdf', 
#'  height = 5, width = 5
#' )
#' 
CanopusListEnrichmentPlot <- function(mmo, feature_list, pthr = 0.05, outdir, height = 5, width = 5){
  term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  sig.canopus <- data.frame(term = character(),  term_level = character(),subsetcount = double(), totalcount = double(), foldenrichment = double(), pval = double(), fdr = double())
  for (term_level in term_levels){
    sig.canopus <- rbind(sig.canopus, CanopusLevelEnrichmentAnal(mmo, feature_list, pthr = pthr, sig = TRUE, term_level = term_level, representation = 'greater'))
  }
  sig.canopus <- sig.canopus |> arrange(dplyr::desc(.data$foldenrichment))
  ggplot(sig.canopus, aes(x = .data$foldenrichment, y = reorder(.data$term, .data$foldenrichment), color = -log(.data$fdr), size = .data$subsetcount)) +
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
#' @export
#' @examplesIf FALSE
#' CanopusListEnrichmentPlot_2(
#'  mmo, feature_list = DAMs_up$control_vs_treatment1.up, 
#'  pthr = 0.05, outdir = 'canopus_enrichment_plot_topn.pdf', 
#'  height = 5, width = 5, topn = 5
#' )
CanopusListEnrichmentPlot_2 <- function(mmo, feature_list, pthr = 0.05, outdir, height = 5, width = 5, topn = 5){
  term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  sig.canopus <- data.frame(term = character(),  term_level = character(),subsetcount = double(), totalcount = double(), foldenrichment = double(), pval = double(), fdr = double())
  for (term_level in term_levels){
    sig.canopus <- rbind(sig.canopus, CanopusLevelEnrichmentAnal(mmo, feature_list, pthr = pthr, sig = TRUE, term_level = term_level, representation = 'greater'))
  }
  sig.canopus$term <- paste(sig.canopus$term, ';', sig.canopus$term_level)
  sig.canopus <- sig.canopus |> dplyr::slice_max(order_by = -.data$pval, n = topn)
  sig.canopus <- sig.canopus |> dplyr::arrange(dplyr::desc(.data$foldenrichment))
  ggplot(sig.canopus, aes(x = .data$foldenrichment, y = reorder(.data$term, .data$foldenrichment), color = -log(.data$fdr), size = .data$subsetcount)) +
    geom_point() +
    scale_color_gradient(low = 'grey', high = 'red') +
    theme_classic()+
    xlim(0,max(sig.canopus$foldenrichment+1))+
    ylab('Chemical Class')

  ggsave(outdir, height = height, width = width)
}

#' Generate a plot for enrichment analysis of Canopus-predicted terms at a specific level using a list of vectors of features
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
#' @export
#' @examplesIf FALSE
#' # Perform enrichment analysis for multiple comparisons using NPC_pathway level
#' comp.list <- list(
#'   comparison1 = DAMs_up$control_vs_treatment1.up,
#'   comparison2 = DAMs_up$control_vs_treatment2.up
#' )
#' CanopusLevelEnrichmentPlot(
#'  mmo, comp.list = comp.list, term_level = 'NPC_pathway', 
#'  pthr = 0.1, representation = 'greater', prefix = 'enrichment_plot', 
#'  height = 5, width = 5
#' )
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
      .data$fdr,
        breaks = c(0,0.001, 0.01, 0.05, 0.1, 1),
        labels = c("***", "**", "*", ".", "")
    ))
  
  enrichment_plot <- ggplot(data = df.EA.sig, aes(x = .data$comp, y = .data$term, label = .data$label))+
    geom_point(aes(size = .data$subsetcount, color = .data$fdr))+
    geom_text()+
    scale_size_area(name = 'Count', max_size = 10)+
    scale_color_gradient2(low = 'red', high = 'grey', mid = 'grey', midpoint = 0.4)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))+
    xlab('Comparisons')+
    ylab('Chemical classes')
  enrichment_plot
  write.csv(df.EA, paste0(prefix, '.csv'), row.names = FALSE)
  write.csv(df.EA.sig, paste0(prefix, '_sig.csv'), row.names = FALSE)
  ggsave(paste0(prefix, '.pdf'), width = width, height = height)
  return(df.EA)
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
#' @export
#' @examplesIf FALSE
#' comp.list <- list(
#'   comparison1 = DAMs_up$control_vs_treatment1.up,
#'   comparison2 = DAMs_up$control_vs_treatment2.up
#' )
#' CanopusAllLevelEnrichmentPlot(
#'  mmo, comp.list = comp.list, terms = 'all_terms', 
#'  pthr = 0.1, representation = 'greater', prefix = 'enrichment_all_levels', 
#'  height = 10, width = 8
#' )
#' CanopusAllLevelEnrichmentPlot(
#'  mmo, comp.list = comp.list, terms = 'NPC', 
#'  pthr = 0.1, representation = 'greater', prefix = 'enrichment_NPC_levels', 
#'  height = 10, width = 8
#' )
#' CanopusAllLevelEnrichmentPlot(
#'  mmo, comp.list = comp.list, terms = 'ClassyFire', 
#'  pthr = 0.1, representation = 'greater', prefix = 'enrichment_ClassyFire_levels', 
#'  height = 10, width = 8
#' )
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
        .data$fdr,
          breaks = c(0,0.001, 0.01, 0.05, 0.1, 1),
          labels = c("***", "**", "*", ".", "")
      ))
  }
  enrichment_plot <- ggplot(data = df.EA.sig, aes(x = .data$comp, y = .data$term, label = .data$label))+
    geom_point(aes(size = .data$subsetcount, color = .data$fdr))+
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
  return(df.EA)
}


#' Metabolite Set Enrichment Analysis (MSEA)
#' 
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
#' @return A data frame containing the MSEA results, including pathway, NES, p-value, and adjusted p-value (FDR)
#' @export 
#' @examplesIf FALSE
#' # Perform MSEA using NPC_class level
#' MSEA(
#'  mmo, feature_name = rownames(DE_results), feature_score = DE_results$log2FoldChange, 
#'  term_level = 'NPC_class', pthr = 0.05, prefix = 'MSEA_NPC_class', 
#'  width = 8, height = 12, sig = FALSE
#' )
MSEA <- function(mmo, feature_name, feature_score, term_level = 'NPC_class', pthr = 0.05, prefix = 'MSEA', width = 8, height = 12, sig = FALSE){
  # Create a named vector of feature scores
  .require_pkg("fgsea")
  ranked_list <- feature_score
  names(ranked_list) <- feature_name
  ranked_list <- sort(ranked_list, decreasing = TRUE)

  # Retrieve metabolite sets based on the specified term level
  if(term_level == 'NPC_class'){
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[['NPC#class']])
  } else if (term_level == 'NPC_superclass'){
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[['NPC#superclass']])
  } else if (term_level == 'NPC_pathway'){
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[['NPC#pathway']])
  } else if (term_level == "ClassyFire_superclass") {
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[['ClassyFire#superclass']])
  } else if (term_level == "ClassyFire_class") {
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[['ClassyFire#class']])
  } else if (term_level == "ClassyFire_subclass") {
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[['ClassyFire#subclass']])
  } else if (term_level == "ClassyFire_level5") {
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[['ClassyFire#level 5']])
  } else if (term_level == "ClassyFire_most_specific") {
    metabolite_sets <- split(mmo$sirius_annot$feature, mmo$sirius_annot[['ClassyFire#most specific class']])
  } else {
    stop("Invalid term level. Please choose a valid term level.")
  }
  msea_res <- fgsea::fgsea(pathways = metabolite_sets,
                       stats    = ranked_list,
                       minSize  = 5,   # minimum number of features in a class
                       maxSize  = 1500,
                       nPermSimple = 10000)
  msea_res <- msea_res |> arrange(.data$padj)
  readr::write_csv(msea_res, paste0(prefix,'_', term_level,'_results.csv'))
  if (sig) {
    msea_res <- msea_res |> filter(.data$padj < pthr)
  }
  ggplot(msea_res, aes(x = reorder(.data$pathway, .data$NES), y = .data$NES)) +
    geom_point(shape = 21, aes(color = .data$padj < 0.05, size = .data$size, fill = -log(.data$padj)), stroke = 1) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    scale_fill_gradient(low = "grey", high = "red") +
    scale_color_manual(values = c("TRUE" = 'black', "FALSE" = 'white')) +
    guides(shape = "none") +
    labs(x = "Metabolite Class", y = "Normalized Enrichment Score (NES)", title = "MSEA Results", color = "-log10(padj)", size = "Set Size") +
    theme_classic() +
    theme(legend.position = "top", axis.text.y = element_text(size = 6)) 
  ggsave(paste0(prefix,'_', term_level,'.pdf'), width = width, height = height)
  return (msea_res)
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
#' @export

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
#' @export
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
#' @export
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
#' @param cor_method The correlation method to use. Options are 'pearson', 'spearman', or 'kendall' (default: 'pearson')
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @return A data frame containing correlation results for each feature, including effect size (correlation coefficient), p-value, and fold change columns for specified comparisons.
#' @export
GetPerformanceFeatureCorrelation <- function(mmo, phenotype, groups, DAM.list, comparisons, cor_method = 'pearson', normalization = 'None'){
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
    cor <- cor.test(combined_df$performance, combined_df$feature_value, method = cor_method)
    pval <- cor[[3]]
    cor <- cor[[4]]
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
#' @export
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
#' @export
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

#' PlotFoldchangeResistanceQuad
#'
#' This function plots the fold change resistance in a quadrant plot, categorizing points into quadrants based on their effect size and fold change.
#' It also performs a binomial test to assess the distribution of points across quadrants.
#' @param performance_regression The regression results data frame containing effect size, fold change, and tag. The output from GetPerformanceFeatureRegression, GetPerformanceFeatureLMM, or GetPerformanceFeatureCorrelation.
#' @param fold_change The name of the fold change column in the performance_regression dataframe
#' @param color A vector of colors for the points in the plot
#' @param output_dir The output file path for the quadrant plot
#' @export
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

#' Generate barplots for each feature and perform ANOVA
#'
#' This function generates bar plots for a specified feature across different groups in the metadata, performing ANOVA and Tukey's HSD test for post-hoc analysis.
#'
#' @param mmo The mmo object containing metadata and feature data
#' @param ID_list A list of feature IDs to analyze. Use FeatureToID() to convert feature names to IDs.
#' @param outdir The output directory to save the bar plots and ANOVA results
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param filter_group A boolean indicating whether to filter the feature values by a specific group list (default: FALSE)
#' @param group_list A list of groups to filter the feature values by, if filter_group is TRUE (default: NULL)
#' @export
#' @examplesIf FALSE
#' AnovaBarPlot(mmo, ID_list = c("ID_1", "ID_2"), outdir = "output_directory", normalization = 'Z')
#' AnovaBarPlot(
#'  mmo, ID_list = c("ID_1", "ID_2"), outdir = "output_directory", normalization = 'Z', 
#'  filter_group = TRUE, group_list = c("Group1", "Group2")
#' )
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
#' @export
#' @examplesIf FALSE
#' ExportFeaturesToCSV(mmo, feature_list = Glucosinolates, normalization = 'Z', output_dir = 'output.csv')
#' ExportFeaturesToCSV(mmo, feature_list = DAMs_up$control_vs_treatment1.up, normalization = 'None', output_dir = 'output.csv')
#' 
ExportFeaturesToCSV <- function(mmo, feature_list, normalization = 'None', output_dir){
  feature <- GetNormFeature(mmo, normalization = normalization) # Get normalized feature data
  # Filter the feature data, annotation, and DA analysis for the list provided
  selected_feature <- feature |> filter(feature %in% feature_list)
  selected_pairwise <- mmo$pairwise |> filter(feature %in% feature_list)
  # Merge all
  merged_df <- merge(mmo$sirius_annot, selected_feature, by = 'feature')
  merged_df <- merge(merged_df, selected_pairwise, by = 'feature')

  write.csv(merged_df, output_dir)
}


#' GetRichness
#' 
#' This function calculates the richness of features for each sample in the mmo object.
#' Richness is defined as the number of non-missing features observed in each sample.
#' @param mmo The mmo object containing feature data and metadata
#' @param filter_feature A boolean indicating whether to filter features based on a provided list (default: FALSE)
#' @param feature_list A list of features to include in the richness calculation if filter_feature is TRUE (default: NULL)
#' @return A data frame containing the richness for each sample, with columns for sample, richness, and group.
#' @export
GetRichness <- function(mmo, filter_feature = FALSE, feature_list = NULL) {
  feature_data <- mmo$feature_data
  if (filter_feature) {
    feature_data <- feature_data |> filter(.data$feature %in% feature_list)
  }
  richness <- apply(feature_data[, -(1:2)], 2, function(x) sum(!is.na(x)))
  
  metadata <- mmo$metadata
  groups <- c()
  for (col in colnames(feature_data)[-c(1, 2)]) {
    groups <- append(groups, metadata[metadata$sample == col, ]$group)
  }
  
  richness_df <- data.frame(sample = colnames(feature_data)[-c(1, 2)], richness = richness, group = groups)
  return(richness_df)
}

#' CalculateCumulativeRichness
#' 
#' This function calculates the cumulative richness of features across groups in the metadata.
#' Cumulative richness is defined as the total number of unique features observed as groups are added sequentially.
#' @param mmo The mmo object containing feature data and metadata
#' @param groups A vector specifying the order of groups to consider for cumulative richness calculation
#' @return A data frame containing the cumulative richness for each group in the specified order, with columns for group and cumulative richness.
#' @export
#' @examplesIf FALSE
#' groups <- c("Control", "Treatment1", "Treatment2")
#' cumulative_richness <- CalculateCumulativeRichness(mmo, groups)
CalculateCumulativeRichness <- function(mmo, groups) {
  feature_data <- mmo$feature_data
  metadata <- mmo$metadata
  cumulative_richness <- numeric(length(groups))
  selected_features <- rep(FALSE, nrow(feature_data))
  for (i in seq_along(groups)) {
    selected_groups <- groups[1:i]
    selected_samples <- metadata |> filter(.data$group %in% selected_groups) |> pull(.data$sample)
    selected_data <- feature_data |> select(all_of(selected_samples))
    selected_features <- selected_features | (rowSums(!is.na(selected_data)) > 0)
    cumulative_richness[i] <- sum(selected_features)
  }
  data.frame(group = groups, cumulative_richness = cumulative_richness)
}

#' BootstrapCumulativeRichness
#' 
#' This function bootstraps the cumulative richness of features across groups in the metadata by randomizing the order of groups.
#' It performs multiple bootstrap iterations to estimate the mean and confidence intervals of cumulative richness at each step.
#' @param mmo The mmo object containing feature data and metadata
#' @param groups A vector of group names from the metadata to consider for cumulative richness calculation
#' @param n_boot The number of bootstrap iterations to perform (default: 1000)
#' @param ci The confidence interval width (e.g., 0.95 for 95% CI) (default: 0.95)
#' @return A data frame containing the mean cumulative richness and confidence intervals for each group index, with columns for group index, mean, lower CI, and upper CI.
#' @export
#' @examplesIf FALSE
#' groups <- c("Control", "Treatment1", "Treatment2")
#' bootstrapped_richness <- BootstrapCumulativeRichness(mmo, groups, n_boot = 1000, ci = 0.95)
BootstrapCumulativeRichness <- function(mmo, groups, n_boot = 1000, ci = 0.95) {
  # Bootstraps cumulative richness by randomizing group order within a direction
  # ci: confidence interval width (e.g., 0.5 for 25%-75%)
  lower_q <- (1 - ci) / 2
  upper_q <- 1 - lower_q
  n_groups <- length(groups)
  boot_mat <- matrix(NA, nrow = n_boot, ncol = n_groups)
  for (i in seq_len(n_boot)) {
    rand_order <- sample(groups)
    boot_mat[i, ] <- CalculateCumulativeRichness(mmo, rand_order)$cumulative_richness
  }
  # Each row is a bootstrap, each column is the cumulative richness after adding that many groups
  boot_df <- data.frame(
    group_index = rep(seq_len(n_groups), times = n_boot),
    bootstrap = rep(seq_len(n_boot), each = n_groups),
    richness = as.vector(t(boot_mat))
  )
  boot_summary <- boot_df |>
    dplyr::group_by(.data$group_index) |>
    dplyr::summarise(
      mean = mean(.data$richness),
      lower = stats::quantile(.data$richness, lower_q),
      upper = stats::quantile(.data$richness, upper_q)
    ) |>
    dplyr::ungroup()
  as.data.frame(boot_summary)
}

#' CalculateNullCumulativeRichness
#' 
#' This function calculates the null model of cumulative richness by randomizing samples regardless of group.
#' It performs multiple bootstrap iterations to estimate the mean and confidence intervals of cumulative richness at each step.
#' @param mmo The mmo object containing feature data and metadata
#' @param n_boot The number of bootstrap iterations to perform (default: 1000)
#' @param n_groups The number of groups to simulate for cumulative richness calculation 
#' @param ci The confidence interval width (e.g., 0.95 for 95% CI) (default: 0.95)
#' @return A data frame containing the mean cumulative richness and confidence intervals for each group index, with columns for group index, mean, lower CI, and upper CI.
#' @export
#' @examplesIf FALSE
#' null_richness <- CalculateNullCumulativeRichness(mmo, n_boot = 1000, n_groups = 5, ci = 0.95)
CalculateNullCumulativeRichness <- function(mmo, n_boot = 1000, n_groups, ci = 0.95) {
  # Null model: randomize samples regardless of group, then add samples one by one
  feature_data <- mmo$feature_data
  metadata <- mmo$metadata
  all_samples <- metadata$sample
  n_features <- nrow(feature_data)
  samples_per_group <- ceiling(length(all_samples) / n_groups)
  boot_mat <- matrix(NA, nrow = n_boot, ncol = n_groups)
  lower_q <- (1 - ci) / 2
  upper_q <- 1 - lower_q
  for (i in seq_len(n_boot)) {
    rand_samples <- sample(all_samples)
    selected_features <- rep(FALSE, n_features)
    for (j in seq_len(n_groups)) {
      end_idx <- min(j * samples_per_group, length(rand_samples))
      selected_data <- feature_data |> select(all_of(rand_samples[1:end_idx]))
      selected_features <- selected_features | (rowSums(!is.na(selected_data)) > 0)
      boot_mat[i, j] <- sum(selected_features)
    }
  }
  boot_df <- data.frame(
    group_index = rep(seq_len(n_groups), times = n_boot),
    bootstrap = rep(seq_len(n_boot), each = n_groups),
    richness = as.vector(t(boot_mat))
  )
  boot_summary <- boot_df |>
    dplyr::group_by(.data$group_index) |>
    dplyr::summarise(
      mean = mean(.data$richness),
      lower = stats::quantile(.data$richness, lower_q),
      upper = stats::quantile(.data$richness, upper_q)
    ) |>
    dplyr::ungroup()
  as.data.frame(boot_summary)
}

#' CalcNormalizedAUC
#' 
#' This function calculates the normalized area under the curve (AUC) for a cumulative richness curve.
#' The normalized AUC is computed by dividing the AUC by the maximum possible area, which is the product of the maximum group index and maximum cumulative richness.
#' @param curve A data frame containing the cumulative richness curve with columns for group index and cumulative richness
#' @return The normalized AUC value
#' @export
#' @examplesIf FALSE
#' curve <- CalculateCumulativeRichness(mmo, group =c("Control", "Treatment1", "Treatment2"))
#' norm_auc <- CalcNormalizedAUC(curve)
CalcNormalizedAUC <- function(curve) {
  curve$group_index <- seq_len(nrow(curve))
  x <- curve$group_index
  y <- curve$cumulative_richness
  auc <- sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
  norm_auc <- auc / (max(x) * max(y))
  norm_auc
}

#' BootCumulRichnessAUC
#' 
#' This function bootstraps the normalized area under the curve (AUC) for cumulative richness by randomizing the order of groups.
#' It performs multiple bootstrap iterations to estimate the distribution of normalized AUC values.
#' @param mmo The mmo object containing feature data and metadata
#' @param groups A vector of group names from the metadata to consider for cumulative richness calculation
#' @param n_boot The number of bootstrap iterations to perform (default: 500)
#' @return A numeric vector containing the normalized AUC values from each bootstrap iteration
#' @export
#' @examplesIf FALSE
#' groups <- c("Control", "Treatment1", "Treatment2")
#' bootstrapped_aucs <- BootCumulRichnessAUC(mmo, groups, n_boot = 500)
BootCumulRichnessAUC <- function(mmo, groups, n_boot = 500) {
  aucs <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    rand_order <- sample(groups)
    curve <- CalculateCumulativeRichness(mmo, rand_order)
    aucs[i] <- CalcNormalizedAUC(curve)
  }
  aucs
}

#' GetFunctionalHillNumber
#'
#' This function calculates the functional Hill number for a given mmo object, normalization method, and distance metric.
#' See https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.18685 for details of the functional Hill number calculation.
#'
#' @param mmo The mmo object containing feature data and metadata
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param q The order of the Hill number to calculate (default: 1). Larger q values give more weight to evenness portion of the hill number over richness.
#' @param distance The distance metric to use for calculating dissimilarity. Options are 'dreams', 'm2ds', or 'cosine' (default: 'dreams')
#' @param filter_feature A boolean indicating whether to filter the feature data by a specific list (default: FALSE)
#' @param feature_list A list of feature names to filter the feature data by, if filter_feature is TRUE (default: NULL)
#' @return A data frame containing the functional Hill number for each group in the metadata, with columns for group and hill number.
#' @export 
#' @examplesIf FALSE
#' hill_number <- GetFunctionalHillNumber(mmo, normalization = 'Z', q = 1, distance = 'dreams', filter_feature = FALSE)
#' hill_number <- GetFunctionalHillNumber(mmo, normalization = 'Z', q = 3, distance = 'dreams', filter_feature = TRUE, feature_list = Glucosinolates)
#
GetFunctionalHillNumber <- function(mmo, normalization = 'None',q = 1, distance = 'dreams', filter_feature = FALSE, feature_list = NULL){
  feature <- GetNormFeature(mmo, normalization = normalization)
  metadata <- mmo$metadata
  distance_matrix <- GetDistanceMat(mmo, distance = distance)
  # Scale the  distance matrix to be between 0 and 1

  if (filter_feature == TRUE){
    id_list <- FeatureToID(mmo, feature_list)
    id_list <- rownames(distance_matrix)[rownames(distance_matrix) %in% id_list]
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
#' This function calculates the Hill numbers for a given mmo object, normalization method, and order of the Hill number without considering feature dissimilarity.
#' 
#'
#' @param mmo The mmo object containing feature data and metadata
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param q The order of the Hill number to calculate (default: 0)
#' @param filter_feature A boolean indicating whether to filter the feature data by a specific list (default: FALSE)
#' @param feature_list A list of feature names to filter the feature data by, if filter_feature is TRUE (default: NULL)
#' @return A data frame containing the Hill number for each group in the metadata, with columns for group and hill number.
#' @export
#' @examplesIf FALSE
#' hill_number <- GetHillNumbers(mmo, normalization = 'Z', q = 1, filter_feature = FALSE)
#' hill_number <- GetHillNumbers(mmo, normalization = 'Z', q = 2, filter_feature = TRUE, feature_list = Glucosinolates)
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
#' Unweighted mode uses Hill numbers without considering feature dissimilarity, while weighted mode uses functional Hill numbers that account for feature dissimilarity.
#' 
#' @param mmo The mmo object containing feature data and metadata
#' @param q The order of the Hill number to calculate (default: 1)
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param mode The mode of diversity calculation. Options are 'weighted' or 'unweighted' for chemical distance(default: 'weighted')
#' @param distance The distance metric to use for calculating dissimilarity. Options are 'dreams', 'm2ds', or 'cosine' (default: 'dreams')
#' @param filter_feature A boolean indicating whether to filter the feature data by a specific list (default: FALSE)
#' @param feature_list A list of feature names to filter the feature data by, if filter_feature is TRUE (default: NULL)
#' @return A data frame containing the alpha diversity for each group in the metadata, with columns for group and alpha diversity value.
#' @export
#' @examplesIf FALSE
#' alpha_diversity <- GetAlphaDiversity(mmo, q = 1, normalization = 'None', 
#'  mode = 'weighted', distance = 'dreams', filter_feature = FALSE)
#' alpha_diversity <- GetAlphaDiversity(mmo, q = 2, normalization = 'Z', 
#'  mode = 'unweighted', filter_feature = TRUE, feature_list = Glucosinolates)
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
#' @export
#' @return A data frame containing the specialization index for each group in the metadata, with columns for group and specialization index.
#' @examplesIf FALSE
#' specialization_index <- GetSpecializationIndex(mmo, normalization = 'None', filter_group = FALSE)
#' specialization_index <- GetSpecializationIndex(mmo, normalization = 'Z', filter_group = TRUE, group_list = c('Control', 'Treatment1'), filter_feature = TRUE)
GetSpecializationIndex <- function(mmo, normalization = 'None', filter_group = FALSE, group_list = NULL, filter_feature = FALSE, feature_list = NULL){
  metadata <- mmo$metadata
  feature <- GetNormFeature(mmo, normalization)

  # All feature or filtered feature
  if (filter_feature == TRUE){
    feature <- feature |> filter(feature %in% feature_list)
  }
  if (filter_group == TRUE){
    samples <- c()
    for (group in group_list){
      samples <- append(samples, metadata |> filter(.data$group == !!group) |> pull(.data$sample))
    }
    feature <- feature |> select(.data$id, .data$feature, all_of(samples))
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
#' The Generalized UniFrac and CSCS method requires a distance matrix of feature dissimilarity, which is calculated using the specified distance metric.
#' Bray-Curtis and Jaccard methods are calculated using the vegan package, not considering feature dissimilarity.
#'
#' @param mmo The mmo object containing feature data and metadata
#' @param method The method of beta diversity calculation. Options are 'Gen.Uni' for Generalized UniFrac, 'bray' for Bray-Curtis, 'jaccard' for Jaccard, or 'CSCS' for Compound Similarity and Chemical structural and compositional similarity (default: 'Gen.Uni')
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param distance The distance metric to use for calculating dissimilarity. Options are 'dreams', 'm2ds', or 'cosine' (default: 'dreams')
#' @param filter_feature A boolean indicating whether to filter the feature data by a specific list (default: FALSE)
#' @param feature_list A list of feature names to filter the feature data by, if filter_feature is TRUE (default: NULL)
#' @param filter_group A boolean indicating whether to filter the feature data by a specific group list (default: FALSE)
#' @param group_list A list of groups to filter the feature data by, if filter_group is TRUE (default: NULL)
#' @return A distance matrix of beta diversity values between samples.
#' @export
#' @examplesIf FALSE
#' beta_diversity <- GetBetaDiversity(mmo, method = 'Gen.Uni', 
#'  normalization = 'None', distance = 'dreams', filter_feature = FALSE)
#' beta_diversity <- GetBetaDiversity(mmo, method = 'bray', 
#'  normalization = 'Z', filter_feature = TRUE, feature_list = Glucosinolates, 
#'  filter_group = TRUE, group_list = c('Control', 'Treatment1'))
GetBetaDiversity <- function(mmo, method = 'Gen.Uni', normalization = 'None', distance = 'dreams', filter_feature = FALSE, feature_list = NULL, filter_group = FALSE, group_list = NULL){
  # Get compound distance and build tree for UniFrac
  scaled_dissimilarity <- GetDistanceMat(mmo, distance = distance) / max(GetDistanceMat(mmo, distance = distance))
  if (filter_feature == TRUE) {
    id_list <- FeatureToID(mmo, feature_list)
    scaled_dissimilarity <- scaled_dissimilarity[id_list, id_list]
  }
  compound_tree <- ape::as.phylo(hclust(as.dist(scaled_dissimilarity), method = "average"))

  # Get feature matrix of relative proportion
  metadata <- mmo$metadata
  feature <- GetNormFeature(mmo, normalization)
  if (filter_group == TRUE){
    samples <- c()
    for (group in group_list){
      samples <- append(samples, metadata |> filter(.data$group == !!group) |> pull(sample))
    }
    feature <- feature |> dplyr::select(.data$id, .data$feature, all_of(samples))
  }
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


#' NMDSplot
#'
#' This function generates a Non-metric Multidimensional Scaling (NMDS) plot based on the provided beta diversity distance matrix.
#' It also performs PERMANOVA analysis to assess the significance of group differences and saves the results to CSV files.
#' @param mmo The mmo object containing metadata
#' @param betadiv The beta diversity distance matrix, output of GetBetaDiversity()
#' @param prefix The prefix for the output files
#' @param width The width of the output NMDS plot (default: 6)
#' @param height The height of the output NMDS plot (default: 6)
#' @param color A vector of colors for the groups in the plot
#' @export
#' @examplesIf FALSE
#' beta_diversity <- GetBetaDiversity(mmo, method = 'Gen.Uni', 
#'  normalization = 'None', distance = 'dreams', filter_feature = FALSE)
#' # Use method = 'bray' or 'jaccard' if you want to use just feature abundance 
#' # without considering feature spectral dissimilarity
#' NMDSplot(mmo, betadiv = beta_diversity, prefix = 'output/NMDS', width = 6, height = 6)
NMDSplot <- function(mmo, betadiv, prefix, width = 6, height = 6, color){
  .require_pkg("vegan")
  .require_pkg("ggrepel")
  metadata <- mmo$metadata
  nmds <- vegan::metaMDS(betadiv, k = 2, try = 50, trymax = 100)
  
  # Extract NMDS coordinates
  nmds_coords <- as.data.frame(vegan::scores(nmds, display = "sites"))
  groups <- c()
  for (row in rownames(nmds_coords)) {
    groups <- append(groups, metadata[metadata$sample == row, ]$group)
  }
  nmds_coords$group <- groups
  
  # Plot NMDS
  ggplot(nmds_coords, aes(x = .data$NMDS1, y = .data$NMDS2, color = .data$group)) +
    geom_point(size = 3) +
    #geom_text_repel(aes(label = group), size = 3) +
    theme_classic() +
    stat_ellipse(level = 0.90) +
    labs(x = "NMDS1", y = "NMDS2") +
    scale_color_manual(values = color) +
    theme(legend.position = "right")
  ggsave(paste0(prefix, '_NMDS.pdf'), height = height, width = width)

  permanova <- permanova_stat(betadiv, mmo$metadata, mode = 'distance')
  write.csv(permanova$permanova_res, paste0(prefix, '_permanova_results.csv'))
  write.csv(as.data.frame(permanova$pairwise_raw), paste0(prefix, '_pairwise_permanova_results.csv'))
  write.csv(as.data.frame(permanova$pairwise_p_matrix), paste0(prefix, '_pairwise_permanova_pvalue_matrix.csv'))
  write.csv(as.data.frame(permanova$pairwise_F_matrix), paste0(prefix, '_pairwise_permanova_F_matrix.csv'))
  write.csv(as.data.frame(permanova$pairwise_R2_matrix), paste0(prefix, '_pairwise_permanova_R2_matrix.csv'))
}

#' PCoAplot
#' 
#' This function generates a Principal Coordinates Analysis (PCoA) plot based on the provided beta diversity distance matrix.
#' It also performs PERMANOVA analysis to assess the significance of group differences and saves the
#' results to CSV files.
#' @param mmo The mmo object containing metadata
#' @param betadiv The beta diversity distance matrix, output of GetBetaDiversity
#' @param prefix The prefix for the output files
#' @param width The width of the output PCoA plot (default: 6
#' @param height The height of the output PCoA plot (default: 6)
#' @param color A vector of colors for the points in the plot
#' @export
#' @examplesIf FALSE
#' beta_diversity <- GetBetaDiversity(mmo, method = 'Gen.Uni', 
#'  normalization = 'None', distance = 'dreams', filter_feature = FALSE)
#' PCoAplot(mmo, betadiv = beta_diversity, prefix = 'output/PCoA', width = 6, height = 6)
PCoAplot <- function(mmo, betadiv, prefix, width = 6, height = 6, color){
  .require_pkg('ape')
  metadata <- mmo$metadata
  pcoa_res <- ape::pcoa(betadiv)
  pcoa_coords <- as.data.frame(pcoa_res$vectors[, 1:2])
  colnames(pcoa_coords) <- c("PCoA1", "PCoA2")
  pcoa_coords$group <- metadata$group[match(rownames(pcoa_coords), metadata$sample)]
  
  ggplot(pcoa_coords, aes(x = .data$PCoA1, y = .data$PCoA2, color = .data$group)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.90) +
    theme_classic() +
    labs(x = "PCoA1", y = "PCoA2") +
    scale_color_manual(values = color) +
    theme(legend.position = "right")
  ggsave(paste0(prefix, '_PCoA.pdf'), height = height, width = width)
  
  permanova <- permanova_stat(betadiv, mmo$metadata, mode = 'distance')
  write.csv(permanova$permanova_res, paste0(prefix, '_permanova_results.csv'))
  write.csv(as.data.frame(permanova$pairwise_raw), paste0(prefix, '_pairwise_permanova_results.csv'))
  write.csv(as.data.frame(permanova$pairwise_p_matrix), paste0(prefix, '_pairwise_permanova_pvalue_matrix.csv'))
  write.csv(as.data.frame(permanova$pairwise_F_matrix), paste0(prefix, '_pairwise_permanova_F_matrix.csv'))
  write.csv(as.data.frame(permanova$pairwise_R2_matrix), paste0(prefix, '_pairwise_permanova_R2_matrix.csv'))
}

#' CalculateGroupBetaDistance
#'
#' This function calculates the beta diversity distance between a reference group and other groups in the metadata.
#'
#' @param mmo The mmo object containing feature data and metadata
#' @param beta_div The beta diversity distance matrix, output of GetBetaDiversity()
#' @param reference_group The name of the reference group to compare against
#' @param groups A vector of group names from the metadata to calculate distances for
#' @return A data frame containing the group names, sample names, and their corresponding beta diversity distances from the reference group.
#' @export
#' @examplesIf FALSE
#' beta_diversity <- GetBetaDiversity(mmo, method = 'Gen.Uni', 
#'  normalization = 'None', distance = 'dreams', filter_feature = FALSE)
#' group_distances <- CalculateGroupBetaDistance(mmo, beta_div = beta_diversity, 
#'  reference_group = 'Control', groups = c('Control', 'Treatment1', 'Treatment2'))
CalculateGroupBetaDistance <- function(mmo, beta_div, reference_group, groups) {
  metadata <- mmo$metadata
  distances <- data.frame(group = character(), distance = numeric())

  for (group in groups) {
      group_samples <- metadata |> filter(group == !!group) |> pull(sample)
      reference_samples <- metadata |> filter(group == !!reference_group) |> pull(sample)

      for (sample in group_samples) {
        for (ref_sample in reference_samples) {
          distance <- beta_div[sample, ref_sample]
          distances <- rbind(distances, data.frame(group = group,sample = sample, distance = distance))
        }
      
    }
  }

  return(distances)
}
