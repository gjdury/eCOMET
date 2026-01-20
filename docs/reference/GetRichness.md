# GetRichness

This function calculates the richness of features for each sample in the
mmo object. Richness is defined as the number of non-missing features
observed in each sample.

## Usage

``` r
GetRichness(mmo, filter_feature = FALSE, feature_list = NULL)
```

## Arguments

- mmo:

  The mmo object containing feature data and metadata

- filter_feature:

  A boolean indicating whether to filter features based on a provided
  list (default: FALSE)

- feature_list:

  A list of features to include in the richness calculation if
  filter_feature is TRUE (default: NULL)

## Value

A data frame containing the richness for each sample, with columns for
sample, richness, and group.
