# Add chemical distance matrices to the mmo object

This function reads cosine, DREAMS, and MS2DeepScore molecular
networking outputs from MZmine, then transform the similarity to
distance and adds the dissimilarity matrices to the mmo object.

## Usage

``` r
AddChemDist(mmo, cos_dir = NULL, dreams_dir = NULL, m2ds_dir = NULL)
```

## Arguments

- mmo:

  The mmo object

- cos_dir:

  Path to the cosine similarity CSV file from MZMine (molecular
  networking)

- dreams_dir:

  Path to the DREAMS similarity CSV file from MZMine (molecular
  networking)

- m2ds_dir:

  Path to the MS2DeepScore similarity CSV file from MZMine (molecular
  networking)

## Value

The mmo object with dissimilarity matrices added (mmo\$cos.dissim,
mmo\$dreams.dissim, mmo\$m2ds.dissim)

## Examples

``` r
if (FALSE) {
mmo <- AddChemDist(mmo,
 cos_dir = "path/to/cosine_similarity.csv",
 dreams_dir = "path/to/dreams_similarity.csv",
 m2ds_dir = "path/to/ms2deepscore_similarity.csv"
)
}
```
