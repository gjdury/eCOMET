# GetBetaDiversity

This function calculates the beta diversity for a given mmo object,
method (Generalized Unifrac, bray, jaccard, or CSCS), normalization
method, distance metric, and optional feature filtering. Then it returns
a distance matrix of beta diversity values. The Generalized UniFrac and
CSCS method requires a distance matrix of feature dissimilarity, which
is calculated using the specified distance metric. Bray-Curtis and
Jaccard methods are calculated using the vegan package, not considering
feature dissimilarity.

## Usage

``` r
GetBetaDiversity(
  mmo,
  method = "Gen.Uni",
  normalization = "None",
  distance = "dreams",
  filter_feature = FALSE,
  feature_list = NULL,
  filter_group = FALSE,
  group_list = NULL
)
```

## Arguments

- mmo:

  The mmo object containing feature data and metadata

- method:

  The method of beta diversity calculation. Options are 'Gen.Uni' for
  Generalized UniFrac, 'bray' for Bray-Curtis, 'jaccard' for Jaccard, or
  'CSCS' for Compound Similarity and Chemical structural and
  compositional similarity (default: 'Gen.Uni')

- normalization:

  The normalization method to use for feature data. Options are 'None',
  'Log', 'Meancentered', or 'Z' (default: 'None')

- distance:

  The distance metric to use for calculating dissimilarity. Options are
  'dreams', 'm2ds', or 'cosine' (default: 'dreams')

- filter_feature:

  A boolean indicating whether to filter the feature data by a specific
  list (default: FALSE)

- feature_list:

  A list of feature names to filter the feature data by, if
  filter_feature is TRUE (default: NULL)

- filter_group:

  A boolean indicating whether to filter the feature data by a specific
  group list (default: FALSE)

- group_list:

  A list of groups to filter the feature data by, if filter_group is
  TRUE (default: NULL)

## Value

A distance matrix of beta diversity values between samples.

## Examples

``` r
if (FALSE) {
beta_diversity <- GetBetaDiversity(mmo, method = 'Gen.Uni',
 normalization = 'None', distance = 'dreams', filter_feature = FALSE)
beta_diversity <- GetBetaDiversity(mmo, method = 'bray',
 normalization = 'Z', filter_feature = TRUE, feature_list = Glucosinolates,
 filter_group = TRUE, group_list = c('Control', 'Treatment1'))
}
```
