# pool_mmo_by_group

Pools sample columns within each group into one pseudo-sample per group.
Feature rows are preserved (no filtering of features). Keeps all other
slots of the mmo object unchanged by copying mmo first.

## Usage

``` r
pool_mmo_by_group(mmo, group_col = "group")
```

## Arguments

- mmo:

  mmo object

- group_col:

  column in mmo\$metadata used for grouping (default: "group")

## Value

mmo object with feature_data containing one column per group and
metadata updated accordingly
