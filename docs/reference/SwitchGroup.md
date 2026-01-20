# Switch the group column in the mmo object

This function switches the group column in the metadata of the mmo
object to a new specified column. The new group column must exist in the
metadata file.

## Usage

``` r
SwitchGroup(mmo, new_group_col)
```

## Arguments

- mmo:

  The mmo object

- new_group_col:

  The name of the new group column in the metadata file

## Value

The mmo object with the updated group column

## Examples

``` r
if (FALSE) {
mmo <- SwitchGroup(mmo, new_group_col = "genotype")
}
```
