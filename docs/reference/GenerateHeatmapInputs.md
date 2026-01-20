# Generate input files to be used for pheatmap from the mmo object

This function generates heatmap inputs from the mmo object, including
fold change or mean values, distance matrix, and row labels for
custom-annotated features.

## Usage

``` r
GenerateHeatmapInputs(
  mmo,
  filter_feature = FALSE,
  feature_list = NULL,
  filter_group = FALSE,
  group_list = NULL,
  summarize = "mean",
  control_group = "ctrl",
  normalization = "None",
  distance = "dreams"
)
```

## Arguments

- mmo:

  The mmo object with sirius annotation and normalized data

- filter_feature:

  Boolean to filter features by feature_list (default: FALSE)

- feature_list:

  A vector of feature names to filter (default: NULL)

- filter_group:

  Boolean to filter groups by group_list (default: FALSE)

- group_list:

  A vector of group names to filter (default: NULL)

- summarize:

  The summarization method to use. Options are 'fold_change' or 'mean'
  (default: 'mean')

- control_group:

  The group to use as control for fold change calculation (default:
  'ctrl')

- normalization:

  The normalization method to use. Options are 'None', 'Log',
  'Meancentered', or 'Z'

- distance:

  The distance metric to use. Options are 'dreams', 'cosine', or 'm2ds'
  (default: 'dreams')

## Value

A list containing the following elements:

- FC_matrix: A matrix of fold change or mean values

- dist_matrix: A distance matrix based on the specified distance metric

- row_label: A vector of row labels for custom-annotated features (See
  AddCustomAnnot()). If no custom annotation is available, feature IDs
  are used.

- heatmap_data: A data frame containing the heatmap data with feature
  IDs and values

## Examples

``` r
if (FALSE) {
# Generate heatmap inputs to visualize fold change values with log normalization and dreams distance
heatmap_inputs <- GenerateHeatmapInputs(
 mmo, summarize = 'fold_change', control_group = 'Control',
 normalization = 'None', distance = 'dreams'
)
# Generate heatmap inputs to visualize mean values
heatmap_inputs <- GenerateHeatmapInputs(
 mmo, summarize = 'mean', normalization = 'None', distance = 'dreams'
)
# The resulting list contains FC_matrix, dist_matrix, row_label, and heatmap_data
# A heatmap can be generated using pheatmap
# 'clustering_distance_rows' option make the dendrogram follows chemical distances of features.
#  -Delete this option to visualize the heatmap following cannonical clustering
pheatmap(mat = heatmap_inputs$FC_matrix,
    cluster_rows = TRUE, #do not change
    clustering_distance_rows = heatmap_inputs$dist_matrix,
    cluster_cols = TRUE,
    clustering_method = "average", #UPGMA
    show_rownames = TRUE,
    show_colnames = TRUE,
    cellwidth = 25,
    cellheight = 0.05,
    treeheight_row = 100,
    fontsize_row = 3,
    fontsize_col = 15,
    scale = 'none',
    annotation_names_row = TRUE,
    labels_row = heatmap_inputs$row_label,
    )
}
```
