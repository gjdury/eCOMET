# Get the distance matrix from the mmo object based on the specified distance metric

This function retrieves the distance matrix from the mmo object based on
the specified distance metric.

## Usage

``` r
GetDistanceMat(mmo, distance = "dreams")
```

## Arguments

- mmo:

  The mmo object

- distance:

  The distance metric to use. Options are 'dreams', 'cosine', or 'm2ds'

## Value

The distance matrix corresponding to the specified distance metric

## Examples

``` r
if (FALSE) {
distance_matrix <- GetDistanceMat(mmo, distance = 'dreams')
}
```
