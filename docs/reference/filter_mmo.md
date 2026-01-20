# Filter an mmo object by samples, groups, and/or features

Subset all components of an mmo object (feature tables, metadata,
distance matrices, and annotations) to a given set of samples, groups,
and/or feature IDs. Filtering is applied consistently across all slots
present in `mmo`.

## Usage

``` r
filter_mmo(
  mmo,
  sample_list = NULL,
  group_list = NULL,
  feature_list = NULL,
  sample_col = "sample",
  group_col = "group",
  drop_empty_feat = TRUE,
  empty_threshold = NULL
)
```

## Arguments

- mmo:

  A list-like mmo object as returned by
  [`GetMZmineFeature()`](https://phytoecia.github.io/eCOMET/reference/GetMZmineFeature.md).

- sample_list:

  Optional character vector of sample IDs (matching `sample_col` in
  `mmo$metadata`) to retain.

- group_list:

  Optional character vector of group labels (matching `group_col` in
  `mmo$metadata`) to retain. Mutually exclusive with `sample_list`.

- feature_list:

  Optional character vector of feature IDs to retain. If `NULL`,
  features are determined from `feature_data` and optionally filtered by
  `drop_empty_feat`.

- sample_col:

  Column name in `mmo$metadata` containing sample IDs. Default is
  `"sample"`.

- group_col:

  Column name in `mmo$metadata` containing group labels. Default is
  `"group"`.

- drop_empty_feat:

  Logical; if `TRUE` (default) drop features with no non-zero values in
  the retained samples.

- empty_threshold:

  Optional numeric threshold used to define “empty” features. If `NULL`
  (default), the smallest positive, non-NA intensity in the retained
  samples (`min_val`) is used. Features are kept if they have at least
  one value \> threshold across retained samples.

## Value

A filtered mmo object with the same structure as `mmo`, but restricted
to the requested samples / groups / features.

## Examples

``` r
if (FALSE) { # \dontrun{
mmo_sub <- filter_mmo(mmo, group_list = c("Species1", "Species2"))
} # }
```
