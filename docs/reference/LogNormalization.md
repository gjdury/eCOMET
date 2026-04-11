# Log-normalize the peak area in the mmo object

This function applies log2 transformation to the peak area in the
feature data of the mmo object. Run ReplaceZero() before this function
to avoid -Inf values.

## Usage

``` r
LogNormalization(mmo, imputed_data = FALSE)
```

## Arguments

- mmo:

  The mmo object

- imputed_data:

  Whether to use imputed feature data (default = FALSE)

## Value

The mmo object with log-normalized feature data (mmo\$log)

## Examples

``` r
if (FALSE) {
mmo <- LogNormalization(mmo)
}
```
