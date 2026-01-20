# Plots PCA and performs PERMANOVA

This function performs PCA analysis and generates a PCA plot with
optional filtering of features and groups. It also conducts PERMANOVA
and saves the results to CSV files.

## Usage

``` r
PCAplot(
  mmo,
  color,
  outdir = "PCA",
  normalization = "Z",
  filter_feature = FALSE,
  feature_list = NULL,
  filter_group = FALSE,
  group_list = NULL,
  label = TRUE,
  save_output = TRUE
)
```

## Arguments

- mmo:

  The mmo object with feature data and metadata

- color:

  A vector of colors for the groups in the plot. Make sure the names
  correspond to the group names in metadata

- outdir:

  The output file path for the PCA plot (default: 'PCA')

- normalization:

  The normalization method to use for feature data. Options are 'None',
  'Log', 'Meancentered', or 'Z' (default: 'Z')

- filter_feature:

  Boolean to filter features by feature_list (default: FALSE)

- feature_list:

  A vector of feature names to filter (default: NULL)

- filter_group:

  Boolean to filter groups by group_list (default: FALSE)

- group_list:

  A vector of group names to filter (default: NULL)

- label:

  Boolean to indicate whether to label points with sample names
  (default: TRUE)

- save_output:

  Boolean; if TRUE (default) write plot (.pdf) and PERMANOVA tables
  using `outdir` as prefix. If FALSE, nothing is written.

## Value

A list with elements `plot` (ggplot), `df` (raw data to generate plots),
and `permanova` (results from `permanova_stat`).

## Examples

``` r
if (FALSE) {
PCAplot(
 mmo, color = c("Control" = "blue", "Treatment1" = "red", "Treatment2" = "green"),
 outdir = 'PCA_plot', normalization = 'None',
 filter_feature = FALSE, filter_group = FALSE, label = FALSE
)
PCAplot(
 mmo, color = c("Control" = "blue", "Treatment1" = "red"),
 outdir = 'PCA_plot', normalization = 'Z',
 filter_feature = TRUE, feature_list = Glucosinolates,
 filter_group = TRUE, group_list = c("Control", "Treatment1"), label = TRUE
)
}
```
