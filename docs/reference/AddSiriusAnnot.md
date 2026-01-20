# Adding annotation from SIRIUS to the mmo object

This function reads SIRIUS structure identification and formula summary
files, and adds the annotations to the mmo object.

## Usage

``` r
AddSiriusAnnot(mmo, canopus_structuredir, canopus_formuladir)
```

## Arguments

- mmo:

  The mmo object

- canopus_structuredir:

  Path to the SIRIUS structure_identification.tsv file

- canopus_formuladir:

  Path to the SIRIUS canopus_formula_summary.tsv file

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
