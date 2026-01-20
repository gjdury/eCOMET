# GetHillNumbers

This function calculates the Hill numbers for a given mmo object,
normalization method, and order of the Hill number without considering
feature dissimilarity.

## Usage

``` r
GetHillNumbers(
  mmo,
  normalization = "None",
  q = 0,
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

  The order of the Hill number to calculate (default: 0)

- filter_feature:

  A boolean indicating whether to filter the feature data by a specific
  list (default: FALSE)

- feature_list:

  A list of feature names to filter the feature data by, if
  filter_feature is TRUE (default: NULL)

## Value

A data frame containing the Hill number for each group in the metadata,
with columns for group and hill number.

## Examples

``` r
if (FALSE) {
hill_number <- GetHillNumbers(mmo, normalization = 'Z', q = 1, filter_feature = FALSE)
hill_number <- GetHillNumbers(mmo, normalization = 'Z', q = 2, filter_feature = TRUE, feature_list = Glucosinolates)
}
```
