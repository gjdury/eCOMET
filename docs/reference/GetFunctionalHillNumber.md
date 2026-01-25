# GetFunctionalHillNumber

This function calculates the functional Hill number for a given mmo
object, normalization method, and distance metric. See
https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.18685 for
details of the functional Hill number calculation.

## Usage

``` r
GetFunctionalHillNumber(
  mmo,
  normalization = "None",
  q = 1,
  distance = "dreams",
  filter_feature = FALSE,
  feature_list = NULL
)
```

## Arguments

- mmo:

  The mmo object containing feature data and metadata

- normalization:

  The normalization method to use for feature data. Options are 'None',
  'Log', 'Meancentered', or 'Z' (default: 'None')

- q:

  The order of the Hill number to calculate (default: 1). Larger q
  values give more weight to evenness portion of the hill number over
  richness.

- distance:

  The distance metric to use for calculating dissimilarity. Options are
  'dreams', 'm2ds', or 'cosine' (default: 'dreams')

- filter_feature:

  A boolean indicating whether to filter the feature data by a specific
  list (default: FALSE)

- feature_list:

  A list of feature names to filter the feature data by, if
  filter_feature is TRUE (default: NULL)

## Value

A data frame containing the functional Hill number for each group in the
metadata, with columns for group and hill number.

## Examples

``` r
if (FALSE) {
hill_number <- GetFunctionalHillNumber(mmo,
                                       normalization = 'Z',
                                       q = 1, distance = 'dreams',
                                       filter_feature = FALSE)
hill_number <- GetFunctionalHillNumber(mmo, normalization = 'Z',
                                       q = 3, distance = 'dreams',
                                       filter_feature = TRUE,
                                       feature_list = Glucosinolates)
}
```
