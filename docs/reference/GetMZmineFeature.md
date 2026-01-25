# Import MZmine feature table and metadata to create a mmo object

Reads an MZmine exported feature table (typically the full feature
table) and a sample metadata table to initialize an mmo object
containing:

- `mmo$feature_data`: feature-by-sample abundance matrix (peak areas)

- `mmo$feature_info`: feature-level annotations (e.g., mz, rt, ranges,
  IDs)

- `mmo$metadata`: sample metadata with standardized `sample` and `group`
  columns

Sample columns in the MZmine table are matched to metadata using fuzzy
string matching on area column names.

## Usage

``` r
GetMZmineFeature(
  mzmine_dir,
  metadata_dir,
  group_col,
  sample_col,
  drop_missing_samples = FALSE,
  mz_col = NULL,
  rt_col = NULL,
  max_distance = 5,
  feature_info_cols = c("id", "rt", "rt_range:min", "rt_range:max", "mz", "mz_range:min",
    "mz_range:max", "feature_group", "ion_identities:iin_id",
    "ion_identities:ion_identities")
)
```

## Arguments

- mzmine_dir:

  Path to the MZmine feature table CSV (should include feature-level
  columns and per-sample area columns)

- metadata_dir:

  Path to the metadata CSV file (must include sample_col and group_col)

- group_col:

  Column name in metadata used for grouping samples (e.g.,
  treatment/species)

- sample_col:

  Column name in metadata used to identify and match samples to MZmine
  area columns

- drop_missing_samples:

  Logical. If FALSE (default), error when metadata samples are missing
  from the MZmine table area columns. If TRUE, drop those samples from
  metadata (with a warning) and continue.

- mz_col:

  Optional m/z column name in the MZmine table (defaults to "mz" or "row
  m/z")

- rt_col:

  Optional RT column name in the MZmine table (defaults to "rt" or "row
  retention time")

- max_distance:

  Maximum edit distance used when fuzzy-matching metadata sample names
  to MZmine area column names (default 5). Lower this for stricter
  matching.

- feature_info_cols:

  Character vector of feature-level columns to retain in
  `mmo$feature_info`. Columns not present in the MZmine table are
  skipped with a warning.

## Value

A mmo object
