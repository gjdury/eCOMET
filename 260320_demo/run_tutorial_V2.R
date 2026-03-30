# ============================================================================
# Treatment-based Study Tutorial V2 — Verification Script
# Runs all code from the V2 tutorial and saves outputs to 260320_demo/
# ============================================================================

library(ecomet)
library(dplyr)
library(stringr)
library(ggplot2)

cat("[INFO] Starting tutorial V2 verification...\n")

# --- Step 1: Locate tutorial data ---
data_dir <- system.file(
  "extdata/tutorials/treatment_based",
  package = "ecomet"
)
stopifnot(nzchar(data_dir))

demo_feature          <- file.path(data_dir, "feature_table_demo.csv")
demo_metadata         <- file.path(data_dir, "metadata_demo.csv")
demo_sirius_formula   <- file.path(data_dir, "canopus_formula_summary.tsv")
demo_sirius_structure <- file.path(data_dir, "structure_identifications.tsv")
demo_dreams           <- file.path(data_dir, "dreams_sim_demo.csv")
gls_db                <- file.path(data_dir, "custom_DB_glucosinolates.csv")

cat("[INFO] Data paths set. data_dir =", data_dir, "\n")

# --- Step 2: Create mmo object ---
mmo <- GetMZmineFeature(
  mzmine_dir   = demo_feature,
  metadata_dir = demo_metadata,
  group_col    = 'group',
  sample_col   = 'sample'
)
cat("[INFO] mmo object created. Features:", nrow(mmo$feature_data), "\n")
cat("[INFO] Samples:", ncol(mmo$feature_data) - 2, "\n")  # minus id, feature columns

# --- Step 3: Preprocessing ---
mmo <- ReplaceZero(mmo, method = 'one')
cat("[INFO] ReplaceZero done.\n")

mmo <- MassNormalization(mmo)
cat("[INFO] MassNormalization done.\n")

mmo <- MeancenterNormalization(mmo)
mmo <- LogNormalization(mmo)
mmo <- ZNormalization(mmo)
cat("[INFO] All normalizations done.\n")

mmo <- FeaturePresence(mmo, threshold = 1)
cat("[INFO] FeaturePresence done.\n")

# --- Step 4: Annotations ---
mmo <- AddSiriusAnnot(
  mmo,
  canopus_structuredir = demo_sirius_structure,
  canopus_formuladir   = demo_sirius_formula
)
cat("[INFO] AddSiriusAnnot done. Annotations:", nrow(mmo$sirius_annot), "\n")

# Filter CANOPUS
mmo <- filter_canopus_annotations(
  mmo,
  pathway_level = "NPC#pathway",
  threshold     = 0.8,
  suffix        = "NPC_pathway_0.8",
  overwrite     = TRUE
)
cat("[INFO] filter_canopus_annotations done.\n")

# Filter COSMIC
filtered_annot <- mmo$sirius_annot_filtered_NPC_pathway_0.8
filtered_annot$ConfidenceScoreApproximate[
  which(filtered_annot$ConfidenceScoreApproximate == "-Infinity")
] <- 0
mmo$sirius_annot_filtered_NPC_pathway_0.8 <- filtered_annot

Cosmic_Scores <- as.numeric(filtered_annot$ConfidenceScoreApproximate)
cosmic_thr <- quantile(x = Cosmic_Scores, probs = 0.9, na.rm = TRUE)
cat("[INFO] COSMIC top 10% threshold:", cosmic_thr, "\n")

mmo <- filter_cosmic_structure(
  mmo,
  input       = "sirius_annot_filtered_NPC_pathway_0.8",
  cosmic_mode = "approx",
  threshold   = as.numeric(cosmic_thr),
  suffix      = "CANOPUS_0.8_COSMIC_Top_10"
)
cat("[INFO] filter_cosmic_structure done.\n")

# Feature vectors
FLVs <- mmo$sirius_annot %>%
  filter(str_detect(`ClassyFire#most specific class`, "Flavonoid")) %>%
  pull(id)
cat("[INFO] Flavonoid features:", length(FLVs), "\n")

# Custom annotation
mmo <- AddCustomAnnot(mmo, DB = gls_db, mztol = 5, rttol = 0.2)
GLSs <- mmo$custom_annot %>%
  filter(lengths(custom_annot) > 0) %>%
  pull(id)
cat("[INFO] Glucosinolate features:", length(GLSs), "\n")

# Chemical distance
mmo <- AddChemDist(mmo, dreams_dir = demo_dreams)
cat("[INFO] AddChemDist done. DreaMS matrix:", nrow(mmo$dreams.dissim), "x", ncol(mmo$dreams.dissim), "\n")

# --- Step 5: Dimensionality reduction ---
setwd("/home/minsoo/software/eCOMET")

colors <- c("ctrl" = "grey", "sl1" = "#fdcdac", "le1" = "#b3e2cd")

# PCA plots
tryCatch({
  PCAplot(mmo, color = colors, outdir = '260320_demo/plot/PCA/PCA_basic', save_output = TRUE)
  cat("[INFO] PCA basic done.\n")
}, error = function(e) cat("[ERROR] PCA basic:", conditionMessage(e), "\n"))

tryCatch({
  PCAplot(mmo, color = colors, outdir = '260320_demo/plot/PCA/PCA_log',
          label = FALSE, normalization = 'Log', save_output = TRUE)
  cat("[INFO] PCA log done.\n")
}, error = function(e) cat("[ERROR] PCA log:", conditionMessage(e), "\n"))

tryCatch({
  PCAplot(mmo, color = colors, outdir = '260320_demo/plot/PCA/PCA_FLV',
          label = FALSE, filter_id = TRUE, id_list = FLVs, save_output = TRUE)
  cat("[INFO] PCA FLV done.\n")
}, error = function(e) cat("[ERROR] PCA FLV:", conditionMessage(e), "\n"))

# Custom PCA
tryCatch({
  pca_res <- PCAplot(mmo, color = colors,
                     outdir = '260320_demo/plot/PCA/PCA_custom',
                     save_output = FALSE)
  custom_pca <- ggplot(pca_res$df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 3) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    labs(title = "Custom PCA Plot")
  ggsave("260320_demo/plot/PCA/PCA_custom.pdf", custom_pca, width = 8, height = 6)
  cat("[INFO] Custom PCA done.\n")
}, error = function(e) cat("[ERROR] Custom PCA:", conditionMessage(e), "\n"))

# PLS-DA
tryCatch({
  PLSDAplot(mmo, color = colors, outdir = '260320_demo/plot/PLSDA/PLSDA_basic.pdf',
            save_output = TRUE)
  cat("[INFO] PLS-DA basic done.\n")
}, error = function(e) cat("[ERROR] PLS-DA basic:", conditionMessage(e), "\n"))

tryCatch({
  PLSDAplot(mmo, color = colors, outdir = '260320_demo/plot/PLSDA/PLSDA_meancentered.pdf',
            normalization = 'Meancentered', save_output = TRUE)
  cat("[INFO] PLS-DA meancentered done.\n")
}, error = function(e) cat("[ERROR] PLS-DA meancentered:", conditionMessage(e), "\n"))

# --- Step 6: DAMs ---
mmo <- PairwiseComp(mmo, group1 = 'ctrl', group2 = 'sl1')
mmo <- PairwiseComp(mmo, group1 = 'ctrl', group2 = 'le1')
cat("[INFO] PairwiseComp done.\n")

DAMs <- GetDAMs(mmo, fc_cutoff = 0.5849625, pval_cutoff = 0.1)
DAMs_up   <- DAMs$DAMs_up
DAMs_down <- DAMs$DAMs_down
cat("[INFO] DAMs extracted.\n")
for (nm in names(DAMs_up)) cat("  DAMs_up$", nm, ":", length(DAMs_up[[nm]]), "\n")
for (nm in names(DAMs_down)) cat("  DAMs_down$", nm, ":", length(DAMs_down[[nm]]), "\n")

# --- Step 7: DAM Visualization ---
tryCatch({
  VolcanoPlot(mmo, comp = 'ctrl_vs_sl1',
              outdir = '260320_demo/plot/Volcano/Volcano_ctrl_vs_sl1.pdf',
              save_output = TRUE)
  cat("[INFO] Volcano ctrl_vs_sl1 done.\n")
}, error = function(e) cat("[ERROR] Volcano ctrl_vs_sl1:", conditionMessage(e), "\n"))

tryCatch({
  VolcanoPlot(mmo, comp = 'ctrl_vs_le1',
              outdir = '260320_demo/plot/Volcano/Volcano_ctrl_vs_le1.pdf',
              topk = 0, save_output = TRUE)
  cat("[INFO] Volcano ctrl_vs_le1 done.\n")
}, error = function(e) cat("[ERROR] Volcano ctrl_vs_le1:", conditionMessage(e), "\n"))

# Venn and Upset
tryCatch({
  VennInput <- list(
    sl1.up = DAMs_up$ctrl_vs_sl1.up,
    le1.up = DAMs_up$ctrl_vs_le1.up
  )
  library(ggvenn)
  venn_plot <- ggvenn(VennInput, stroke_size = 0.5, set_name_size = 4,
                      show_percentage = FALSE) +
    theme(legend.position = "none")
  ggsave("260320_demo/plot/Venn/Venn_Upreg.pdf", venn_plot, height = 5, width = 5)
  cat("[INFO] Venn done.\n")
}, error = function(e) cat("[ERROR] Venn:", conditionMessage(e), "\n"))

tryCatch({
  library(UpSetR)
  pdf("260320_demo/plot/Upset/Upset_Upreg.pdf", 7, 5)
  upset(fromList(VennInput), nsets = 10, nintersects = 20,
        order.by = 'freq', mainbar.y.label = 'Features in Set',
        line.size = 1, point.size = 4, shade.color = 'white',
        text.scale = 1, show.numbers = FALSE)
  dev.off()
  cat("[INFO] UpSet done.\n")
}, error = function(e) { try(dev.off(), silent = TRUE); cat("[ERROR] UpSet:", conditionMessage(e), "\n") })

# --- Step 8: Heatmaps ---
library(pheatmap)

tryCatch({
  heatmap_inputs <- GenerateHeatmapInputs(
    mmo, summarize = 'fold_change', control_group = 'ctrl',
    normalization = 'None', distance = 'dreams'
  )
  # Handle Inf/NaN in fold-change matrix
  fc_mat <- as.matrix(heatmap_inputs$FC_matrix)
  finite_max <- max(abs(fc_mat[is.finite(fc_mat)]), na.rm = TRUE)
  fc_mat[is.infinite(fc_mat) & fc_mat > 0] <-  finite_max
  fc_mat[is.infinite(fc_mat) & fc_mat < 0] <- -finite_max
  fc_mat[is.nan(fc_mat)] <- 0
  pheatmap(
    mat = fc_mat,
    clustering_distance_rows = heatmap_inputs$dist_matrix,
    clustering_method = "average",
    cellwidth = 100, cellheight = 0.3,
    treeheight_row = 100, fontsize_row = 3, fontsize_col = 15,
    scale = 'none',
    filename = "260320_demo/plot/Heatmap/Heatmap_FC_dreams.pdf",
    width = 10, height = 10
  )
  cat("[INFO] Heatmap FC dreams done.\n")
}, error = function(e) cat("[ERROR] Heatmap FC dreams:", conditionMessage(e), "\n"))

tryCatch({
  heatmap_inputs_z <- GenerateHeatmapInputs(
    mmo, summarize = 'mean', normalization = 'Z', distance = 'dreams'
  )
  pheatmap(
    mat = heatmap_inputs_z$FC_matrix,
    clustering_method = "average",
    cellwidth = 100, cellheight = 0.3,
    treeheight_row = 100, fontsize_row = 3, fontsize_col = 15,
    scale = 'none',
    filename = "260320_demo/plot/Heatmap/Heatmap_Mean_Z_clustering.pdf",
    width = 10, height = 10
  )
  cat("[INFO] Heatmap Mean Z done.\n")
}, error = function(e) cat("[ERROR] Heatmap Mean Z:", conditionMessage(e), "\n"))

tryCatch({
  if (length(GLSs) > 0) {
    heatmap_inputs_GLS <- GenerateHeatmapInputs(
      mmo, summarize = 'mean', normalization = 'Z', distance = 'dreams',
      filter_id = TRUE, id_list = GLSs
    )
    pheatmap(
      mat = heatmap_inputs_GLS$FC_matrix,
      clustering_method = "average",
      cellwidth = 100, cellheight = 8,
      treeheight_row = 100, fontsize_row = 8, fontsize_col = 15,
      scale = 'none', labels_row = heatmap_inputs_GLS$row_label,
      filename = "260320_demo/plot/Heatmap/Heatmap_Mean_Z_GLS.pdf",
      width = 10, height = 10
    )
    cat("[INFO] Heatmap GLS done.\n")
  } else {
    cat("[WARN] No glucosinolate features found, skipping GLS heatmap.\n")
  }
}, error = function(e) cat("[ERROR] Heatmap GLS:", conditionMessage(e), "\n"))

# --- Step 9: Enrichment ---
tryCatch({
  CanopusListEnrichmentPlot(
    mmo, DAMs_up$ctrl_vs_sl1.up, pthr = 0.1,
    outdir = '260320_demo/plot/Enrichment/sl1_up_enrichment.pdf',
    height = 6, width = 6
  )
  cat("[INFO] CanopusListEnrichmentPlot done.\n")
}, error = function(e) cat("[ERROR] CanopusListEnrichmentPlot:", conditionMessage(e), "\n"))

tryCatch({
  CanopusListEnrichmentPlot_2(
    mmo, DAMs_up$ctrl_vs_sl1.up, pthr = 0.1,
    outdir = '260320_demo/plot/Enrichment/sl1_up_enrichment_top10.pdf',
    topn = 10, height = 6, width = 6
  )
  cat("[INFO] CanopusListEnrichmentPlot_2 done.\n")
}, error = function(e) cat("[ERROR] CanopusListEnrichmentPlot_2:", conditionMessage(e), "\n"))

tryCatch({
  CanopusLevelEnrichmentPlot(
    mmo, DAMs_up, term_level = 'NPC_class', pthr = 0.1,
    outdir = '260320_demo/plot/Enrichment/DAMs_up_NPC_class'
  )
  cat("[INFO] CanopusLevelEnrichmentPlot done.\n")
}, error = function(e) cat("[ERROR] CanopusLevelEnrichmentPlot:", conditionMessage(e), "\n"))

tryCatch({
  CanopusAllLevelEnrichmentPlot(
    mmo, DAMs_up, term_level = 'NPC', pthr = 0.1,
    outdir = '260320_demo/plot/Enrichment/DAMs_up_all_NPC',
    width = 8, height = 12
  )
  cat("[INFO] CanopusAllLevelEnrichmentPlot NPC done.\n")
}, error = function(e) cat("[ERROR] CanopusAllLevelEnrichmentPlot NPC:", conditionMessage(e), "\n"))

tryCatch({
  CanopusAllLevelEnrichmentPlot(
    mmo, DAMs_up, term_level = 'ClassyFire', pthr = 0.1,
    outdir = '260320_demo/plot/Enrichment/DAMs_up_all_ClassyFire',
    width = 8, height = 12
  )
  cat("[INFO] CanopusAllLevelEnrichmentPlot ClassyFire done.\n")
}, error = function(e) cat("[ERROR] CanopusAllLevelEnrichmentPlot ClassyFire:", conditionMessage(e), "\n"))

# --- Step 10: Correlation ---
tryCatch({
  sl_cor <- ScreenFeaturePhenotypeCorrelation(
    mmo, phenotype = 'sl', groups = c('sl1'),
    model = 'spearman', normalization = 'Z'
  )
  cat("[INFO] Correlation screening done. Results:", nrow(sl_cor), "features\n")
  write.csv(sl_cor, "260320_demo/output/sl_correlation_results.csv", row.names = FALSE)
}, error = function(e) cat("[ERROR] Correlation:", conditionMessage(e), "\n"))

# --- Step 11: Save mmo ---
tryCatch({
  SaveMMO(mmo, '260320_demo/output/mmo.RData')
  cat("[INFO] mmo saved to 260320_demo/output/mmo.RData\n")
}, error = function(e) cat("[ERROR] SaveMMO:", conditionMessage(e), "\n"))

cat("\n[DONE] Tutorial V2 verification complete.\n")
cat("[INFO] sessionInfo():\n")
print(sessionInfo())
