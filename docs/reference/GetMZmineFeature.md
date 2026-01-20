# Import mzmine feature data and metadata to create a mmo object

This function reads mzmine feature data and metadata from specified
directories, to initiate a mmo object containing feature data and
metadata

Import mzmine feature data and metadata to create a mmo object

## Usage

``` r
GetMZmineFeature(
  mzmine_dir,
  metadata_dir,
  group_col,
  sample_col,
  mz_col = NULL,
  rt_col = NULL
)
```

## Arguments

- mzmine_dir:

  Path to the mzmine feature data CSV file

- metadata_dir:

  Path to the metadata CSV file (must include sample_col and group_col)

- group_col:

  Column name in the metadata file used for grouping samples together
  i.e into treatments or species.

- sample_col:

  Column in metadata file used to identify and match individual samples

- mz_col:

  Optional m/z column name (defaults to "mz" or "row m/z")

- rt_col:

  Optional RT column name (defaults to "rt" or "row retention time")

## Value

A mmo object containing the feature data and metadata
