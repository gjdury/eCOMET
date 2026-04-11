# Mean-center the peak area in the mmo object

This function applies mean-centering to the peak area in the feature
data of the mmo object. Mean-centering is performed per feature (row)
across samples. Features with zero variance are returned as all zeros
and are reported in a warning.

## Usage

``` r
MeancenterNormalization(mmo, imputed_data = FALSE)
```

## Arguments

- mmo:

  The mmo object

- imputed_data:

  Whether to use imputed feature data (default = FALSE)

## Value

The mmo object with mean-centered feature data stored in
`mmo$meancentered`
