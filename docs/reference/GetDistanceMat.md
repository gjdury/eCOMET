# Get the distance matrix from the mmo object based on the specified distance metric

Retrieve a feature distance matrix from the mmo object

## Usage

``` r
GetDistanceMat(mmo, distance = "dreams")
```

## Arguments

- mmo:

  The mmo object

- distance:

  Name of the distance matrix to retrieve. Built-in options are
  `'dreams'`, `'cosine'`, and `'m2ds'`. Any name passed to
  `AddCustomDist(mmo, name = ...)` is also valid.

## Value

The distance matrix (numeric matrix with feature ID row/col names).

## Details

Looks up a dissimilarity matrix stored in the mmo object by name. Works
with the three built-in matrices added by
[`AddChemDist()`](https://phytoecia.github.io/eCOMET/reference/AddChemDist.md)
as well as any custom matrix added via
[`AddCustomDist()`](https://phytoecia.github.io/eCOMET/reference/AddCustomDist.md).

## Examples

``` r
if (FALSE) {
distance_matrix <- GetDistanceMat(mmo, distance = 'dreams')
distance_matrix <- GetDistanceMat(mmo, distance = 'tanimoto')  # custom
}
```
