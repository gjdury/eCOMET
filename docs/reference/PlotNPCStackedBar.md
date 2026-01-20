# PlotNPCStackedBar

This function generates a stacked bar plot showing the count of features
in each group categorized by NPC_pathway. It uses the mmo object with
sirius annotation and normalized data. Make sure you don't run
ReplaceZero() before using this function, as it may remove
presence/absence information.

## Usage

``` r
PlotNPCStackedBar(
  mmo,
  group_col,
  outdir,
  width = 6,
  height = 3,
  save_output = TRUE
)
```

## Arguments

- mmo:

  The mmo object with sirius annotation and normalized data

- group_col:

  The column name in metadata to use for grouping samples

- outdir:

  The output file path for the stacked bar plot (e.g.,
  'NPC_stacked_bar.png')

- width:

  The width of the output plot

- height:

  The height of the output plot

- save_output:

  boolean, whether to save the output plot

## Value

A list containing the stacked bar plot and the data used to generate it

## Examples

``` r
if (FALSE) {
PlotNPCStackedBar(
 mmo, group_col = 'treatment',
 outdir = 'NPC_stacked_bar.png', width = 6, height = 3
)
}
```
