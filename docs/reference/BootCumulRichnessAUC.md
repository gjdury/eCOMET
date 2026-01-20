# BootCumulRichnessAUC

This function bootstraps the normalized area under the curve (AUC) for
cumulative richness by randomizing the order of groups. It performs
multiple bootstrap iterations to estimate the distribution of normalized
AUC values.

## Usage

``` r
BootCumulRichnessAUC(mmo, groups, n_boot = 500)
```

## Arguments

- mmo:

  The mmo object containing feature data and metadata

- groups:

  A vector of group names from the metadata to consider for cumulative
  richness calculation

- n_boot:

  The number of bootstrap iterations to perform (default: 500)

## Value

A numeric vector containing the normalized AUC values from each
bootstrap iteration

## Examples

``` r
if (FALSE) {
groups <- c("Control", "Treatment1", "Treatment2")
bootstrapped_aucs <- BootCumulRichnessAUC(mmo, groups, n_boot = 500)
}
```
