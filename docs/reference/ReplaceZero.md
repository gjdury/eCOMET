# Replace zero and NA values in the mmo object

This function replaces zero and NA values in the feature data of the mmo
object. Run this function before MassNormalization(),
LogNormalization(), MeancenterNormalization(), or ZNormalization().

## Usage

``` r
ReplaceZero(mmo, method = "one")
```

## Arguments

- mmo:

  The mmo object

- method:

  The method to use for replacement. Options are 'one' (replace with 1)
  or 'half_mean' (replace with half of the smallest non-zero value in
  the row)

## Value

The mmo object with replaced values in the feature data
(mmo\$feature_data)

## Examples

``` r
if (FALSE) {
mmo <- ReplaceZero(mmo, method = 'one')
}
```
