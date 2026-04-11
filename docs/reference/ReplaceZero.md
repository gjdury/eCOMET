# \#' Replace zero and NA values in the mmo object

This function replaces zero values in the feature table of an mmo
object. Imputed data are stored in mmo\$imputed_feature_data, to be used
for downsteam analyses including PairwiseComp(). Note that imputation
affects normalizations (Log-transformation, etc.), as well as chemical
diversity calculations that uses presence/absence.

## Usage

``` r
ReplaceZero(mmo, method = c("one", "half_min"))
```

## Arguments

- mmo:

  An mmo object containing `feature_data`

- method:

  Replacement method:

  - "one": replace zeros and NA values with 1

  - "half_min": replace zeros and NA values with half of the smallest
    non-zero value in each feature (row)

## Value

The updated mmo object
