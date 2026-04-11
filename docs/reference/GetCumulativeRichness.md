# GetCumulativeRichness

This function calculates the cumulative richness of features across
groups in the metadata. Cumulative richness is defined as the total
number of unique features observed as groups are added sequentially.

## Usage

``` r
GetCumulativeRichness(mmo, groups)
```

## Arguments

- mmo:

  The mmo object containing feature data and metadata

- groups:

  A vector specifying the order of groups to consider for cumulative
  richness calculation

## Value

A data frame containing the cumulative richness for each group in the
specified order, with columns for group and cumulative richness.

## Examples

``` r
if (FALSE) {
groups <- c("Control", "Treatment1", "Treatment2")
cumulative_richness <- GetCumulativeRichness(mmo, groups)
}
```
