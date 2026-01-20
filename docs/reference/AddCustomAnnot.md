# Add custom annotations to an mmo object

Match features to a custom DB by m/z (ppm) and RT (minutes) tolerances
and attach a list-column of candidate compound IDs per feature.

## Usage

``` r
AddCustomAnnot(mmo, DB_file, mztol = 5, rttol = 0.5)
```

## Arguments

- mmo:

  An `mmo` object created by
  [`GetMZmineFeature()`](https://phytoecia.github.io/eCOMET/reference/GetMZmineFeature.md).

- DB_file:

  CSV path with at least columns `compound`, `mz`, `rt`.

- mztol:

  m/z tolerance in ppm (default 5).

- rttol:

  RT tolerance in minutes (default 0.5).

## Value

The same `mmo` object with `mmo$custom_annot` (id, feature,
custom_annot).

## Examples

``` r
if (FALSE) {
mmo <- AddCustomAnnot(mmo, DB_file = "path/to/custom_db.csv", mztol = 5, rttol = 0.5)
}
```
