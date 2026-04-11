# CalcNormalizedAUC

This function calculates the normalized area under the curve (AUC) for a
cumulative richness curve. The normalized AUC is computed by dividing
the AUC by the maximum possible area, which is the product of the
maximum group index and maximum cumulative richness.

## Usage

``` r
CalcNormalizedAUC(curve)
```

## Arguments

- curve:

  A data frame containing the cumulative richness curve with columns for
  group index and cumulative richness

## Value

The normalized AUC value

## Examples

``` r
if (FALSE) {
curve <- GetCumulativeRichness(mmo, group =c("Control", "Treatment1", "Treatment2"))
norm_auc <- CalcNormalizedAUC(curve)
}
```
