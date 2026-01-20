# Build a feature dendrogram from a dissimilarity matrix

Given a feature-to-feature dissimilarity matrix (e.g., DREAMS, cosine,
MS2DeepScore), cluster the features and optionally write a Newick tree
for visualization in tools like iTOL.

## Usage

``` r
FeatureDendrogram(
  mmo,
  distance = "dreams",
  features = NULL,
  method = "average",
  save_output = FALSE,
  outprefix = "feature_tree",
  layout = c("base", "circular"),
  npl_col = "NPC#class",
  group_col = "group",
  add_group_bars = TRUE,
  abundance_fun = c("mean", "sum")
)
```

## Arguments

- mmo:

  The mmo object containing one or more dissimilarity matrices (e.g.,
  `dreams.dissim`)

- distance:

  Which dissimilarity matrix to use; one of "dreams", "cos", or "m2ds"
  (default: "dreams")

- features:

  Optional character vector of feature IDs to include; if NULL, all
  features in the matrix are used

- method:

  Clustering linkage method passed to
  [`stats::hclust`](https://rdrr.io/r/stats/hclust.html) (default:
  "average")

- save_output:

  Logical; if TRUE, write a Newick tree (`<outprefix>.nwk`) and a PDF
  dendrogram (`<outprefix>.pdf`)

- outprefix:

  File prefix for outputs when `save_output` is TRUE (default:
  "feature_tree")

- layout:

  Choose "base" (original rectangular tree) or "circular" (ggtree-based
  circular layout with annotation rings)

- npl_col:

  Column in `mmo$sirius_annot` used to color tips (e.g., "NPC#class");
  ignored if missing

- group_col:

  Column in `mmo$metadata` that defines sample groups for tip bars

- add_group_bars:

  Logical; when TRUE, draws stacked bars at the tips summarizing feature
  abundance by `group_col`

- abundance_fun:

  Aggregation used for group bars: "mean" (default) or "sum" over
  samples within a group

## Value

A list with `hclust` (the clustering), `dendrogram`, `phylo` (ape::phylo
object), and when `layout = "circular"` a `ggplot` named `plot`

## Examples

``` r
if (FALSE) {
tree <- FeatureDendrogram(mmo, distance = "dreams", features = NULL, method = "average",
                          save_output = TRUE, outprefix = "dreams_tree")
}
```
