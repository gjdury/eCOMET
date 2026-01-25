# HCplot

Hierarchical clustering of samples from a precomputed beta-diversity
matrix and plotting as a phylogram with tip labels colored by species
(or any grouping column).

## Usage

``` r
HCplot(
  mmo,
  betadiv,
  outdir,
  group_col = "Species_binomial",
  sample_col = "sample",
  hclust_method = "average",
  palette = "Dark 3",
  cex = 0.6,
  width = 10,
  height = 7,
  save_output = TRUE
)
```

## Arguments

- mmo:

  The mmo object containing metadata in mmo\$metadata

- betadiv:

  The beta diversity distance matrix, output of GetBetaDiversity()

- outdir:

  Output prefix for files (e.g., "output/HC"). If save_output=TRUE a PDF
  is saved.

- group_col:

  Metadata column name used to color tips (default: "Species_binomial")

- sample_col:

  Metadata column name containing sample IDs (default: "sample")

- hclust_method:

  hclust linkage method (default: "average"; alternatives:
  "complete","ward.D2")

- palette:

  Qualitative palette name for colorspace::qualitative_hcl (default:
  "Dark 3")

- cex:

  Tip label size (default: 0.6)

- width:

  PDF width (default: 10)

- height:

  PDF height (default: 7)

- save_output:

  Whether to save the plot to PDF (default: TRUE)

## Value

A list containing: hc (hclust), phy (phylo), tip_df (mapping), colors
(named palette)

## Details

This function is intended for visualization (no cluster significance is
implied).

## Examples

``` r
if (FALSE) {
bet <- GetBetaDiversity(mmo, method='bray',
        normalization='Log', distance='dreams',
         filter_feature=FALSE)
HCplot(mmo,
       betadiv = bet,
       outdir = "output/HC_dreams_bray")
}
```
