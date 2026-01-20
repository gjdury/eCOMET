# FeaturePerformanceRegression

This function performs regression analysis of a specific feature against
a phenotype performance in the metadata. It can use linear mixed models
(LMM), simple linear regression (LM), or Pearson correlation.

## Usage

``` r
FeaturePerformanceRegression(
  mmo,
  target,
  phenotype,
  groups,
  model = "lmm",
  normalization = "Z",
  outdir = "FeaturePerformanceRegression",
  width = 6,
  height = 6,
  save_output = TRUE
)
```

## Arguments

- mmo:

  The mmo object with feature data and metadata

- target:

  The name of the feature to analyze

- phenotype:

  The name of the phenotype performance in the metadata

- groups:

  A vector of group names from the metadata containing performance data

- model:

  The type of regression model to use. Options are 'lmm' for linear
  mixed model, 'lm' for simple linear regression, or 'pearson' for
  Pearson correlation (default: 'lmm')

- normalization:

  The normalization method to use for feature data. Options are 'None',
  'Log', 'Meancentered', or 'Z' (default: 'Z')

- outdir:

  The directory to save the output files

- width:

  The width of the output plot in inches (defa ult: 8)

- height:

  The height of the output plot in inches (default: 12)

- save_output:

  A logical value indicating whether to save the output plot (default:
  TRUE)

## Value

A list of the plot and the raw data
