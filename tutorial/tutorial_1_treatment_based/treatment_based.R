# pak::pak("phytoecia/eCOMET")
library(ecomet)
library(dplyr)
library(stringr)
#Load demo data
setwd('/Users/dlforrister/Library/CloudStorage/OneDrive-SmithsonianInstitution/One_Drive_BackUps_Local_Mac_Files/CODE_GIT_HUB_2017_Aug_31/eCOMET/tutorial/tutorial_1_treatment_based/')

demo_feature <- "raw_data/feature_table_demo.csv"
demo_metadata <- "raw_data/metadata_demo.csv"
demo_sirius_formula <- "raw_data/canopus_formula_summary.tsv"
demo_sirius_structure <- "raw_data/structure_identifications.tsv"
demo_dreams <- "raw_data/dreams_sim_demo.csv"
gls_db <- "raw_data/custom_DB_glucosinolates.csv"

# Initialize eCOMET object
mmo <- GetMZmineFeature(mzmine_dir=demo_feature, metadata_dir = demo_metadata, group_col = 'group',sample_col = "sample")
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

###########################################################################

# Plot PCA and run permanova
colors <- c("ctrl" = "grey", "sl1" = "#fdcdac", "le1" = "#b3e2cd")
PCAplot(mmo, color = colors, outdir = 'plot/PCA')
PCAplot(mmo, color = colors, outdir = 'plot/PCA_log', label = FALSE, normalization = 'Log') # Data normalization option
PCAplot(mmo, color = colors, outdir = 'plot/PCA_FLV', label = FALSE, filter_feature = TRUE, feature_list = FLVs) # Feature filtering option

###########################################################################

# Plot PLS-DA
PLSDAplot(mmo, color = colors, outdir = 'plot/PLSDA.pdf')
PLSDAplot(mmo, color = colors, outdir = 'plot/PLSDA_meancentered.pdf', normalization = 'Meancentered')

###########################################################################

# Run pairwise comparison
mmo <- PairwiseComp(mmo, group1 = 'ctrl', group2 = 'sl1')
mmo <- PairwiseComp(mmo, group1 = 'ctrl', group2 = 'le1')
# Extract DAMs
DAMs <- GetDAMs(mmo, fc_cutoff = 0.5849625, pval_cutoff = 0.1) # log2(1.5) = 0.5849625
DAMs_up <- DAMs$DAMs_up
DAMs_down <- DAMs$DAMs_down
head(DAMs_up)

###########################################################################

# Venn diagram and UpSet of DAMs
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

###########################################################################

# Volcano plot
VolcanoPlot(mmo, comp = 'ctrl_vs_sl1', outdir = 'plot/Volcano_ctrl_vs_sl1.pdf')
VolcanoPlot(mmo, comp = 'ctrl_vs_le1', outdir = 'plot/Volcano_ctrl_vs_le1.pdf', topk = 0) # Remove label by setting topk = 0

###########################################################################

# Heatmap
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

###########################################################################
# CANOPUS term enrichment analysis
# For a single set of features, a detailed enrichment plot can be generated
# There are two plotting styles available
CanopusListEnrichmentPlot(mmo, DAMs_up$ctrl_vs_sl1.up, pthr = 0.1, outdir = 'plot/sl1_up_enrichment.pdf', height = 6, width = 6)
CanopusListEnrichmentPlot_2(mmo, DAMs_up$ctrl_vs_sl1.up, pthr = 0.1, outdir = 'plot/sl1_up_enrichment_2.pdf', topn = 10, height = 6, width = 6)

# For a list of sets features, a summary enrichment plot can be generated
# The summary enrichment plot can be generated for either a single level of CANOPUS classification (8.2.1) or for all levels (8.2.2)
# For a single level of CANOPUS classification
term_levels <- c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_most_specific', 'ClassyFire_level5', 'ClassyFire_subclass', 'ClassyFire_class', 'ClassyFire_superclass')
CanopusLevelEnrichmentPlot(mmo, DAMs_up, term_level = 'NPC_class', pthr = 0.1, prefix = 'plot/DAMs_up_NPC_class.pdf')


# For all levels of CANOPUS classification
# All levels, or only NPC or ClassyFire
CanopusAllLevelEnrichmentPlot(mmo, DAMs_up, term_level = 'NPC', pthr = 0.1, prefix = 'plot/DAMs_up_all_terms', width = 8, height = 12)
CanopusAllLevelEnrichmentPlot(mmo, DAMs_up, term_level = 'ClassyFire', pthr = 0.1, prefix = 'plot/DAMs_up_all_terms', width = 8, height = 12)



###########################################################################

# Alpha diversity by Hill numbers
# q : Hill number order, incresing q gives more weights to abundant features
# mode : weighted for the chemical structure (use distance argument), unweighted for no chemical weight
alphadiv <- GetAlphaDiversity(mmo, q = 3, mode = 'weighted', distance = 'dreams', filter_feature = FALSE, feature_list = NULL)
alphadiv <- GetAlphaDiversity(mmo, q = 3, mode = 'unweighted', filter_feature = FALSE, feature_list = NULL) # Unweighted alpha diversity
alphadiv <- GetAlphaDiversity(mmo, q = 3, mode = 'weighted', distance = 'dreams', filter_feature = TRUE, feature_list = GLSs) # Feature filtering option
# 9.1.1 Plot the alpha diversity
ggplot(alphadiv, aes(x = group, y = hill_number)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = 0.5) +
  theme_classic() +
  labs(title = 'Alpha Diversity by Hill Numbers', x = 'Group', y = 'Hill Number') +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('plots/dreams_q3_log.pdf', width = 5, height = 5)
# Test for significant differences between groups with ANOVA
anova <- anova_tukey_dunnett(alphadiv, 'hill_number ~ group')
write_anova(anova, 'plots/anova_hill_number.csv')

###########################################################################

# Calculate beta diversity for different normalizations and methods
# For unweighted beta diversity, Bray-Curtis or Jaccard distance can be used
bray <- GetBetaDiversity(mmo, method = 'bray', normalization = 'None', filter_feature = FALSE, feature_list = NULL)
jaccard <- GetBetaDiversity(mmo, method = 'jaccard', normalization = 'None', filter_feature = FALSE, feature_list = NULL)

# For weighted beta diversity, Generalized UniFrac and CSCS can be used
guni <- GetBetaDiversity(mmo, method = 'Gen.Uni', normalization = 'Log', distance = 'dreams', filter_feature = FALSE, feature_list = NULL)
guni.0 <- guni[,,'d_0'] # Generalized UniFrac with alpha 0
guni.05 <- guni[,,'d_0.5'] # Generalized UniFrac with alpha 0.5
guni.1 <- guni[,,'d_1'] # Generalized UniFrac with alpha 1
CSCS <- GetBetaDiversity(mmo, method = 'CSCS', normalization = 'Log', distance = 'dreams', filter_feature = FALSE, feature_list = NULL)

# The beta diversity can be visualized using NMDS or PCoA plots
NMDSplot(mmo, betadiv = bray, prefix = 'plots/NMDS_bray', color = colors)
PCoAplot(mmo, betadiv = guni.05, prefix = 'plots/PCoA_guni.05', color = colors)

# Or, distance against a group can be extracted
group_distances <- CalculateGroupBetaDistance(mmo, beta_div = guni.05, reference_group = 'ctrl', groups = c('le1', 'sl1'))
ggplot(group_distances, aes(x = group, y = distance)) +
      geom_boxplot(outlier.shape = NA) +
      geom_beeswarm(size = 0.5) +
      theme_classic() +
      labs(x = "Group", y = "Beta Diversity")
ggsave('plots/betadiv/group_dist.pdf', height = 6, width = 6)

###########################################################################
# Visualization of individual features by Barplot
GLSs_id <- FeatureToID(mmo, GLSs)
AnovaBarplot(mmo, ID_list = GLS_id, outdir = 'plots/Barplot_GLSs', normalization = 'None')

###########################################################################
# Export features of interest to CSV
ExportFeaturesToCSV(mmo, feature_list = DAMs_up$ctrl_vs_sl1.up, normalization = 'None', output_dir = 'output/sl1_up_features.csv')
