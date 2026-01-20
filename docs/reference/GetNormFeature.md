# Retrieve feature data from the mmo object, with normalization options

This function retrieves the feature data from the mmo object based on
the specified normalization method.

## Usage

``` r
GetNormFeature(mmo, normalization)
```

## Arguments

- mmo:

  The mmo object

- normalization:

  The normalization method to use. Options are 'None', 'Log',
  'Meancentered', or 'Z'

## Value

The feature data corresponding to the specified normalization method

## Examples

``` r
if (FALSE) {
feature_data <- GetNormFeature(mmo, normalization = 'Log')
}
```
