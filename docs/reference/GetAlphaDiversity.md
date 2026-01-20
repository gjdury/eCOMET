# GetAlphaDiversity

This function calculates the alpha diversity for a given mmo object,
order of the Hill number, normalization method, mode (weighted or
unweighted), distance metric, and optional feature filtering. Unweighted
mode uses Hill numbers without considering feature dissimilarity, while
weighted mode uses functional Hill numbers that account for feature
dissimilarity.

## Usage

``` r
GetAlphaDiversity(
  mmo,
  q = 1,
  normalization = "None",
  mode = "weighted",
  distance = "dreams",
  filter_feature = FALSE,
  feature_list = NULL
)
```

## Arguments

- mmo:

  The mmo object containing feature data and metadata

- q:

  The order of the Hill number to calculate (default: 1)

- normalization:

  The normalization method to use for feature data. Options are 'None',
  'Log', 'Meancentered', or 'Z' (default: 'None')

- mode:

  The mode of diversity calculation. Options are 'weighted' or
  'unweighted' for chemical distance(default: 'weighted')

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

A data frame containing the alpha diversity for each group in the
metadata, with columns for group and alpha diversity value.

## Examples

``` r
if (FALSE) {
alpha_diversity <- GetAlphaDiversity(mmo, q = 1, normalization = 'None',
 mode = 'weighted', distance = 'dreams', filter_feature = FALSE)
alpha_diversity <- GetAlphaDiversity(mmo, q = 2, normalization = 'Z',
 mode = 'unweighted', filter_feature = TRUE, feature_list = Glucosinolates)
}
```
