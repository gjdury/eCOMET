# Adding annotation from SIRIUS to the mmo object

This function reads SIRIUS structure identification and formula summary
files, and adds the annotations to the mmo object.

## Usage

``` r
AddSiriusAnnot(
  mmo,
  canopus_structuredir,
  canopus_formuladir,
  filter_annot = FALSE,
  filter_threshold = 0.5
)
```

## Arguments

- mmo:

  The mmo object

- canopus_structuredir:

  Path to the SIRIUS structure_identification.tsv file

- canopus_formuladir:

  Path to the SIRIUS canopus_formula_summary.tsv file

- filter_annot:

  Logical. If TRUE, filter the annotations by probability threshold in
  CANOPUS.

- filter_threshold:

  Numeric between 0 and 1. The probability threshold for filtering
  annotations.

## Value

The mmo object with SIRIUS annotations added

## Examples

``` r
if (FALSE) {
mmo <- AddSiriusAnnot(mmo,
 canopus_structuredir = "path/to/structure_identification.tsv",
 canopus_formuladir = "path/to/canopus_formula_summary.tsv"
)
}
```
