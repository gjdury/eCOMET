# CalculateNullCumulativeRichness

This function calculates the null model of cumulative richness by
randomizing samples regardless of group. It performs multiple bootstrap
iterations to estimate the mean and confidence intervals of cumulative
richness at each step.

## Usage

``` r
CalculateNullCumulativeRichness(mmo, n_boot = 1000, n_groups, ci = 0.95)
```

## Arguments

- mmo:

  The mmo object containing feature data and metadata

- n_boot:

  The number of bootstrap iterations to perform (default: 1000)

- n_groups:

  The number of groups to simulate for cumulative richness calculation

- ci:

  The confidence interval width (e.g., 0.95 for 95% CI) (default: 0.95)

## Value

A data frame containing the mean cumulative richness and confidence
intervals for each group index, with columns for group index, mean,
lower CI, and upper CI.

## Examples

``` r
if (FALSE) {
null_richness <- CalculateNullCumulativeRichness(mmo, n_boot = 1000, n_groups = 5, ci = 0.95)
}
```
