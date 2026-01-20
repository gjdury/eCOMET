# Treatment-based_study_tutorial

``` r
# pak::pak("phytoecia/eCOMET")
library(ecomet)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(stringr)
library(here)
#> here() starts at /Users/dlforrister/Library/CloudStorage/OneDrive-SmithsonianInstitution/One_Drive_BackUps_Local_Mac_Files/CODE_GIT_HUB_2017_Aug_31/eCOMET
here()
#> [1] "/Users/dlforrister/Library/CloudStorage/OneDrive-SmithsonianInstitution/One_Drive_BackUps_Local_Mac_Files/CODE_GIT_HUB_2017_Aug_31/eCOMET"
```

## Background

In this tutorial, we will demonstrate how to use the **eCOMET** package
for analyzing metabolomics data from a treatment-based study.

- Treatment-based studies involve comparing metabolite profiles between
  different treatment groups (e.g., control vs treated).

- The tutorial covers data preprocessing, normalization, statistical
  analysis, and visualization techniques commonly used in
  treatment-based metabolomics studies.

- The tutorial files are metabolomics analysis files of *Arabidopsis
  thaliana* (Col-0) attacked by two different herbivores (*Spodoptera
  litura*; sl, and *Lipaphis erysimi*; le).

- Eight replicates were sampled and analyzed for each group. See
  metadata.csv for more details.

- **Tutorial data**

  eCOMET is distributed with a small set of example data files that are
  installed automatically with the package. These files are included
  specifically so that the tutorials can be run immediately, without
  downloading additional data or setting file paths by hand.

  The example datasets represent simplified versions of the feature
  tables, metadata, and annotation outputs used in a typical eCOMET
  analysis workflow. They are meant for demonstration and testing, not
  as complete research datasets.

``` r

# Locate tutorial data shipped with the eCOMET package
data_dir <- system.file(
  "extdata/tutorials/treatment_based",
  package = "ecomet"
)

stopifnot(nzchar(data_dir))  # fail loudly if package/data not installed

# Define file paths
demo_feature <- file.path(data_dir, "feature_table_demo.csv")
demo_metadata <- file.path(data_dir, "metadata_demo.csv")
demo_sirius_formula <- file.path(data_dir, "canopus_formula_summary.tsv")
demo_sirius_structure <- file.path(data_dir, "structure_identifications.tsv")
demo_dreams <- file.path(data_dir, "dreams_sim_demo.csv")
gls_db <- file.path(data_dir, "custom_DB_glucosinolates.csv")
```

## 1. Create an `mmo` object

Following steps are performed to create an `mmo` object and add various
normalizations and annotations. Inspect the structure of the ‘mmo’
object using ‘summary(mmo)’’ after each step to see how the object is
updated.

``` r
# Initialize eCOMET object
mmo <- GetMZmineFeature(mzmine_dir=demo_feature, metadata_dir = demo_metadata, group_col = 'group', sample_col = 'sample')
# Add normalizations
mmo <- ReplaceZero(mmo, method = 'one') # Replace 0 and NA values by 1
mmo <- MassNormalization(mmo) # Normalize peak area by sample mass in metadata
mmo <- MeancenterNormalization(mmo) # Add mean-centered area
mmo <- LogNormalization(mmo) # Add log-transformed area
mmo <- ZNormalization(mmo) # Add Zscore

# Add SIRIUS annotation
mmo <- AddSiriusAnnot(mmo, canopus_structuredir = demo_sirius_structure, canopus_formuladir = demo_sirius_formula)
# Make a vector of flavonoids using SIRIUS annotation
FLVs <- mmo$sirius_annot %>% filter(str_detect(mmo$sirius_annot[['ClassyFire#most specific class']], "Flavonoid")) %>% pull(feature)

# Add custom annotation using inhouse glucosinolate library
mmo <- AddCustomAnnot(mmo, DB = gls_db, mztol = 5, rttol = 0.2)
# Make a vector of glucosinolates using custom annotation
GLSs <- mmo$custom_annot %>% filter(lengths(custom_annot) > 0) %>% pull(feature)

# Add Dreams distance
mmo <- AddChemDist(mmo, dreams_dir = demo_dreams)
```

## 2. Plot dimensionality reduction plots

### PCA plot and PERMANOVA

A PCA plot is commonly used to visualize the overall distribution of
groups. Following PERMANOVA test can be performed to check if the groups
are significantly different.
[`PCAplot()`](https://phytoecia.github.io/eCOMET/reference/PCAplot.md)
function use mmo to perform PCA, plot them, and perform PERMANOVA test.

``` r
# Set colors for groups
colors <- c("ctrl" = "grey", "sl1" = "#fdcdac", "le1" = "#b3e2cd")
# PCA plot
# A plot and permanova output files will be generated in the specified outdir
PCAplot(mmo, color = colors, outdir = 'plot/PCA')
PCAplot(mmo, color = colors, outdir = 'plot/PCA_log', label = FALSE, normalization = 'Log') # Data normalization option
PCAplot(mmo, color = colors, outdir = 'plot/PCA_FLV', label = FALSE, filter_feature = TRUE, feature_list = FLVs) # Feature filtering option
```

Alternatively, you can use the output of PCAplot to generate plot with
your own aesthetics.

``` r
# Get PCA plot inputs
pca_res <- PCAplot(mmo, color = colors, outdir = 'plot/PCA')
# Generate plot with your own aesthetics
pdf("plot/PCA.pdf", width = 10, height = 10)
ggplot(pca_res$df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "PCA plot", x = "PC1", y = "PC2")
dev.off()
```

Likewise, all plot-generating functions in eCOMET return a list of
dataframes and other objects that can be used to generate plots with
your own aesthetics. Use `summary(pca_res)` to see the structure of the
output and try this for other plot-generating functions. \### PLS-DA
plot

``` r
# Plot PLS-DA
PLSDAplot(mmo, color = colors, outdir = 'plot/PLSDA.pdf')
PLSDAplot(mmo, color = colors, outdir = 'plot/PLSDA_meancentered.pdf', normalization = 'Meancentered') # Data normalization option
```

## 3. Identify differentially accumulated metabolites (DAMs)

Many analyses targets to find **Differentially Accumulated Metabolites
(DAMs)**. DAMs can be defined by thresholds of log2-fold change and
adjusted p-value. Those two metrics are calculated by following code.
Note that the divisor group is at the left. -
[`PairwiseComp()`](https://phytoecia.github.io/eCOMET/reference/PairwiseComp.md)
function performs pairwise comparison between two groups using t-test
and fold change. The results are stored in mmo\$pairwise_comp. -
[`GetDAMs()`](https://phytoecia.github.io/eCOMET/reference/GetDAMs.md)
function extracts DAMs based on user-defined cutoffs for p-value and
fold change. - Here, we compare each herbivore treatment group (sl1 and
le1) to the control group (ctrl). - Log2(1.5) = 0.5849625 is used as
fold change cutoff and 0.1 as p-value cutoff to extract more DAMs

``` r
# Run pairwise comparison
mmo <- PairwiseComp(mmo, group1 = 'ctrl', group2 = 'sl1')
mmo <- PairwiseComp(mmo, group1 = 'ctrl', group2 = 'le1')
# Extract DAMs
DAMs <- GetDAMs(mmo, fc_cutoff = 0.5849625, pval_cutoff = 0.1) # log2(1.5) = 0.5849625
DAMs_up <- DAMs$DAMs_up
DAMs_down <- DAMs$DAMs_down
head(DAMs_up)
```

### 3.1. Volcano plot

A volcano plot can be used to visualize the overall distribution of
features based on fold change and p-value.

``` r
# Volcano plot
VolcanoPlot(mmo, comp = 'ctrl_vs_sl1', outdir = 'plot/Volcano_ctrl_vs_sl1.pdf')
VolcanoPlot(mmo, comp = 'ctrl_vs_le1', outdir = 'plot/Volcano_ctrl_vs_le1.pdf', topk = 0) # Remove label by setting topk = 0
```

### 3.2. Venn diagram and upset plot

A Venn diagram or upset plot can be used to visualize the overlap of
DAMs between different comparisons.

``` r
# Define input list
VennInput <- list(
  sl1.up = DAMs_up$ctrl_vs_sl1.up,
  le1.up = DAMs_up$ctrl_vs_le1.up
)
# Venn diagram
require(ggvenn)
library(ggvenn)
ggvenn(VennInput, stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE) +
  theme(legend.position = "none") 
ggsave("plot/Venn_Upreg.pdf", height = 5, width = 5)
# UpSet plot
require(UpSetR)
library(UpSetR)
pdf("plot/Upset_Upreg.pdf", 7, 5)
upset(fromList(VennInput), nsets=10, nintersects=20,order.by='freq', mainbar.y.label='Features in Set', line.size=1, point.size=4, shade.color='white', text.scale=1, show.numbers=FALSE)
dev.off()
```

## 4. Heatmap

To visualize the relative abundance of features across samples, a
heatmap can be generated. The features can be either clusterd using
hierarchical clustering or by chemical distances (e.g., Dreams
distance), following idea of Qemistree (Tripathi et al., 2021, Nat Chem.
Biol.).

- The input matrix for heatmap can be either fold change between each
  treatment and control or mean normalized values across samples.
- Various normalizations can be used (None, Log, Meancentered, and Z).
- Groups or features can be filtered by setting ‘filter_group’ or
  ‘filter_feature’ to TRUE and providing a list of groups or features.
- The input data for the heatmap is generated by
  ‘GenerateHeatmapInputs()’ function, then the heatmap is plotted using
  ‘pheatmap’ package.

``` r
library(pheatmap)
# Generate input for heatmap
# distance is one of the chemical distance (cosine, m2ds, and dreams) for clustering rows
# The values can be either fold_change or mean (use option 'summarize')
heatmap_inputs <- GenerateHeatmapInputs(
  mmo, summarize = 'fold_change', control_group = 'ctrl', 
  normalization = 'None', distance = 'dreams'
) # This generates fold change matrix between each treatment and control

# The resulting list contains FC_matrix, dist_matrix, row_label, and heatmap_data
# A heatmap can be generated using pheatmap
pdf("plot/Heatmap_FC_dreams.pdf", width = 10, height = 10)
pheatmap(mat = heatmap_inputs$FC_matrix, 
     clustering_distance_rows = heatmap_inputs$dist_matrix,  
     clustering_method = "average", 
     cellwidth = 100,
     cellheight = 0.3,
     treeheight_row = 100,
     fontsize_row = 3,
     fontsize_col = 15,
     scale = 'none'
)
dev.off()

# Either, you can visualize mean normalized values across samples
heatmap_inputs <- GenerateHeatmapInputs(
  mmo, summarize = 'mean', normalization = 'Z', distance = 'dreams'
) 
# 'clustering_distance_rows' option make the dendrogram follows chemical distances of features. 
#  -Delete this option to visualize the heatmap following cannonical clustering
pdf("plot/Heatmap_Mean_Z_clustering.pdf", width = 10, height = 10)
pheatmap(mat = heatmap_inputs$FC_matrix, 
     #clustering_distance_rows = heatmap_inputs$dist_matrix,  # Delete this option to visualize the heatmap following cannonical clustering 
     clustering_method = "average", #UPGMA
     cellwidth = 100,
     cellheight = 0.3,
     treeheight_row = 100,
     fontsize_row = 3,
     fontsize_col = 15,
     scale = 'none'
)
dev.off()

# Feature filtering option
# Visualize only glucosinolates
# As glucosinolate annotations are stored in mmo$custom_annot (by AddCustomAnnot() function),
# the compound names are stored in heatmap_inputs$row_label and can be used as row names
heatmap_inputs_GLS <- GenerateHeatmapInputs(
  mmo, summarize = 'mean', normalization = 'Z', distance = 'dreams',
  filter_feature = TRUE, feature_list = GLSs
) 
pdf("plot/Heatmap_Mean_Z_GLS.pdf", width = 10, height = 10)
pheatmap(mat = heatmap_inputs_GLS$FC_matrix, 
     #clustering_distance_rows = heatmap_inputs_GLS$dist_matrix,  # Delete this option to visualize the heatmap following cannonical clustering 
     clustering_method = "average", #UPGMA
     cellwidth = 100,
     cellheight = 8,
     treeheight_row = 100,
     fontsize_row = 8,
     fontsize_col = 15,
     scale = 'none',
     annotation_names_row = TRUE,
     labels_row = heatmap_inputs_GLS$row_label
)
dev.off()
```

## 5. CANOPUS class enrichment analysis

Biological questions ask which class of chemical compounds are enriched
in a set of compounds of interest (e.g., DAMs from above). This is
analogue to the Gene Ontology enrichment analysis performed in
transcriptomics. In MMO, NPC and Classyfire terms annotated by Canopus
of SIRIUS are used to perform chemical class enrichment analysis of
given list of features. The enrichment score of each term is calculated
to plot the number of each term and the significance.

``` r
# For a single set of features, a detailed enrichment plot can be generated
# There are two plotting styles available
CanopusListEnrichmentPlot(mmo, DAMs_up$ctrl_vs_sl1.up, pthr = 0.1, outdir = 'plot/sl1_up_enrichment.pdf', height = 6, width = 6)
CanopusListEnrichmentPlot_2(mmo, DAMs_up$ctrl_vs_sl1.up, pthr = 0.1, outdir = 'plot/sl1_up_enrichment_2.pdf', topn = 10, height = 6, width = 6)

# For a list of sets features, a summary enrichment plot can be generated
# The summary enrichment plot can be generated for either a single level of CANOPUS classification (8.2.1) or for all levels (8.2.2)
# For a single level of CANOPUS classification
term_levels <- c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_most_specific', 'ClassyFire_level5', 'ClassyFire_subclass', 'ClassyFire_class', 'ClassyFire_superclass')

CanopusLevelEnrichmentPlot(mmo, DAMs_up, term_level = 'NPC_class', pthr = 0.1, prefix = 'plot/DAMs_up_NPC_class')


# For all levels of CANOPUS classification
# All levels, or only NPC or ClassyFire
CanopusAllLevelEnrichmentPlot(mmo, DAMs_up, term_level = 'NPC', pthr = 0.1, prefix = 'plot/DAMs_up_all_terms', width = 8, height = 12)
CanopusAllLevelEnrichmentPlot(mmo, DAMs_up, term_level = 'ClassyFire', pthr = 0.1, prefix = 'plot/DAMs_up_all_terms', width = 8, height = 12)
```

## 6. Correlation analysis

To find phenotype-linked metabolties, we can use correlation analysis.
In this tutorial, we will use correlation between amount of each feature
and herbivore performance in bioassay to find defense compounds.

``` r
# We will first screen for all features to find whether there is any correlation between amount of each feature and herbivore performance in bioassay
sl_cor <- ScreenFeaturePhenotypeCorrelation(mmo, phenotype = 'sl', groups = c('sl1'), model = 'spearman', normalization = 'Z')
head(sl_cor)
```

## 7. Sharing MMO object

MMO object can be shared with other users by saving it to a file and
loading it from a file.

``` r
# Save MMO object
SaveMMO(mmo, 'mmo.RData')

# Load MMO object
mmo <- LoadMMO('mmo.RData')
```
