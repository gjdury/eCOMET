# PLS-DA plot with feature loadings

PLS-DA plot with feature loadings

## Usage

``` r
PLSDAplot(
  mmo,
  color,
  topk = 10,
  outdir,
  normalization = "Z",
  filter_feature = FALSE,
  feature_list = NULL,
  filter_group = FALSE,
  group_list = NULL,
  save_output = TRUE
)
```

## Arguments

- mmo:

  The mmo object with feature data and metadata

- color:

  A vector of colors for the groups in the plot. Make sure the names
  correspond to the group names in metadata

- topk:

  The number of top features to display in the plot (default: 10)

- outdir:

  The output file path for the PLS-DA plot

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

- save_output:

  Boolean; if TRUE (default) write plot (.pdf) and loadings tables using
  `outdir` as prefix. If FALSE, nothing is written.

## Value

A list with elements `plot` (ggplot), `df` (raw data to generate plots),
and `loadings` (loadings for PLSDA).

## Examples

``` r
if (FALSE) {
PLSDAplot(
 mmo, color = c("Control" = "blue", "Treatment1" = "red", "Treatment2" = "green"),
 topk = 10, outdir = 'PLSDA_plot.pdf', normalization = 'Z',
 filter_feature = FALSE, filter_group = FALSE
)
PLSDAplot(
 mmo, color = c("Control" = "blue", "Treatment1" = "red"),
 topk = 5, outdir = 'PLSDA_plot.pdf', normalization = 'Log',
 filter_feature = TRUE, feature_list = Glucosinolates,
 filter_group = TRUE, group_list = c("Control", "Treatment1")
)
}
```
