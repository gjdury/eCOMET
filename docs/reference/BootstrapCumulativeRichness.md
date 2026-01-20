# BootstrapCumulativeRichness

This function bootstraps the cumulative richness of features across
groups in the metadata by randomizing the order of groups. It performs
multiple bootstrap iterations to estimate the mean and confidence
intervals of cumulative richness at each step.

## Usage

``` r
BootstrapCumulativeRichness(mmo, groups, n_boot = 1000, ci = 0.95)
```

## Arguments

- mmo:

  The mmo object containing feature data and metadata

- groups:

  A vector of group names from the metadata to consider for cumulative
  richness calculation

- n_boot:

  The number of bootstrap iterations to perform (default: 1000)

- ci:

  The confidence interval width (e.g., 0.95 for 95% CI) (default: 0.95)

## Value

A data frame containing the mean cumulative richness and confidence
intervals for each group index, with columns for group index, mean,
lower CI, and upper CI.

## Examples

``` r
if (FALSE) {
groups <- c("Control", "Treatment1", "Treatment2")
bootstrapped_richness <- BootstrapCumulativeRichness(mmo, groups, n_boot = 1000, ci = 0.95)
}
```
