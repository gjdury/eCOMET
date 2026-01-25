# Package index

## Core dataset processing

- [`GetMZmineFeature()`](https://phytoecia.github.io/eCOMET/reference/GetMZmineFeature.md)
  : Import MZmine feature table and metadata to create a mmo object
- [`AddSiriusAnnot()`](https://phytoecia.github.io/eCOMET/reference/AddSiriusAnnot.md)
  : Adding annotation from SIRIUS to the mmo object
- [`AddChemDist()`](https://phytoecia.github.io/eCOMET/reference/AddChemDist.md)
  : Add chemical distance matrices to the mmo object
- [`AddCustomAnnot()`](https://phytoecia.github.io/eCOMET/reference/AddCustomAnnot.md)
  : Add custom annotations to an mmo object
- [`LogNormalization()`](https://phytoecia.github.io/eCOMET/reference/LogNormalization.md)
  : Log-normalize the peak area in the mmo object
- [`MassNormalization()`](https://phytoecia.github.io/eCOMET/reference/MassNormalization.md)
  : Use sample mass in the metadata file to normalize the peak area
- [`MeancenterNormalization()`](https://phytoecia.github.io/eCOMET/reference/MeancenterNormalization.md)
  : Mean-center the peak area in the mmo object
- [`ZNormalization()`](https://phytoecia.github.io/eCOMET/reference/ZNormalization.md)
  : Z-normalize the peak area in the mmo object
- [`ReplaceZero()`](https://phytoecia.github.io/eCOMET/reference/ReplaceZero.md)
  : Replace zero and NA values in the mmo object
- [`ReorderGroups()`](https://phytoecia.github.io/eCOMET/reference/ReorderGroups.md)
  : Reorder samples in the mmo object based on group order
- [`SwitchGroup()`](https://phytoecia.github.io/eCOMET/reference/SwitchGroup.md)
  : Switch the group column in the mmo object
- [`FeaturePresence()`](https://phytoecia.github.io/eCOMET/reference/FeaturePresence.md)
  : Convert feature abundances to presence / absence

## Filterng MMO and associated MGF

- [`filter_mmo()`](https://phytoecia.github.io/eCOMET/reference/filter_mmo.md)
  : Filter an mmo object by samples, groups, and/or features
- [`filter_mgf_to_mmo()`](https://phytoecia.github.io/eCOMET/reference/filter_mgf_to_mmo.md)
  : Filter an MGF file to keep only spectra for features present in
  mmo\$feature_data\$id
- [`annotate_feature_info_ms2_from_mgf()`](https://phytoecia.github.io/eCOMET/reference/annotate_feature_info_ms2_from_mgf.md)
  : Annotate mmo\$feature_info with MS2 presence and MS2 block counts
  from an MGF
- [`filter_canopus_annotations()`](https://phytoecia.github.io/eCOMET/reference/filter_canopus_annotations.md)
  : Filter CANOPUS / SIRIUS annotations in an ecomet mmo object by
  probability threshold
- [`filter_cosmic_structure()`](https://phytoecia.github.io/eCOMET/reference/filter_cosmic_structure.md)
  : Filter SIRIUS structure (CSI:FingerID) annotations by COSMIC
  confidence score

## Basic statistics and visualization

- [`PairwiseComp()`](https://phytoecia.github.io/eCOMET/reference/PairwiseComp.md)
  : Perform pairwise comparison between two groups in the mmo object
- [`VolcanoPlot()`](https://phytoecia.github.io/eCOMET/reference/VolcanoPlot.md)
  : Volcano plot for visualizing differential metabolite analysis
  results
- [`PCAplot()`](https://phytoecia.github.io/eCOMET/reference/PCAplot.md)
  : Plots PCA and performs PERMANOVA
- [`PLSDAplot()`](https://phytoecia.github.io/eCOMET/reference/PLSDAplot.md)
  : PLS-DA plot with feature loadings
- [`AnovaBarPlot()`](https://phytoecia.github.io/eCOMET/reference/AnovaBarPlot.md)
  : Generate barplots for each feature and perform ANOVA
- [`GenerateHeatmapInputs()`](https://phytoecia.github.io/eCOMET/reference/GenerateHeatmapInputs.md)
  : Generate input files to be used for pheatmap from the mmo object
- [`GetDAMs()`](https://phytoecia.github.io/eCOMET/reference/GetDAMs.md)
  : Generates lists of DAMs (Differentially Accumulated Metabolites) for
  each comparison in the mmo object
- [`GetGroupMeans()`](https://phytoecia.github.io/eCOMET/reference/GetGroupMeans.md)
  : Calculate group means from the mmo object
- [`GetLog2FoldChange()`](https://phytoecia.github.io/eCOMET/reference/GetLog2FoldChange.md)
  : Calculate log2 fold change for a given control group
- [`HCplot()`](https://phytoecia.github.io/eCOMET/reference/HCplot.md) :
  HCplot

## Chemical class analysis

- [`CanopusLevelEnrichmentAnal()`](https://phytoecia.github.io/eCOMET/reference/CanopusLevelEnrichmentAnal.md)
  : Enrichment analysis for Canopus-predicted terms
- [`CanopusListEnrichmentPlot()`](https://phytoecia.github.io/eCOMET/reference/CanopusListEnrichmentPlot.md)
  : Generate a plot for enrichment analysis of Canopus-predicted terms
- [`CanopusListEnrichmentPlot_2()`](https://phytoecia.github.io/eCOMET/reference/CanopusListEnrichmentPlot_2.md)
  : Generate a plot for enrichment analysis of Canopus-predicted terms
  across multiple levels
- [`CanopusLevelEnrichmentPlot()`](https://phytoecia.github.io/eCOMET/reference/CanopusLevelEnrichmentPlot.md)
  : Generate a plot for enrichment analysis of Canopus-predicted terms
  at a specific level using a list of vectors of features
- [`CanopusAllLevelEnrichmentPlot()`](https://phytoecia.github.io/eCOMET/reference/CanopusAllLevelEnrichmentPlot.md)
  : Generate a plot for enrichment analysis of Canopus-predicted terms
  across all levels
- [`MSEA()`](https://phytoecia.github.io/eCOMET/reference/MSEA.md) :
  Metabolite Set Enrichment Analysis (MSEA)

## Chemical diversity Analysis

- [`GetRichness()`](https://phytoecia.github.io/eCOMET/reference/GetRichness.md)
  : GetRichness
- [`GetAlphaDiversity()`](https://phytoecia.github.io/eCOMET/reference/GetAlphaDiversity.md)
  : GetAlphaDiversity
- [`GetFunctionalHillNumber()`](https://phytoecia.github.io/eCOMET/reference/GetFunctionalHillNumber.md)
  : GetFunctionalHillNumber
- [`GetHillNumbers()`](https://phytoecia.github.io/eCOMET/reference/GetHillNumbers.md)
  : GetHillNumbers
- [`GetBetaDiversity()`](https://phytoecia.github.io/eCOMET/reference/GetBetaDiversity.md)
  : GetBetaDiversity
- [`NMDSplot()`](https://phytoecia.github.io/eCOMET/reference/NMDSplot.md)
  : NMDSplot
- [`PCoAplot()`](https://phytoecia.github.io/eCOMET/reference/PCoAplot.md)
  : PCoAplot
- [`CalculateGroupBetaDistance()`](https://phytoecia.github.io/eCOMET/reference/CalculateGroupBetaDistance.md)
  : CalculateGroupBetaDistance
- [`GetSpecializationIndex()`](https://phytoecia.github.io/eCOMET/reference/GetSpecializationIndex.md)
  : GetSpecializationIndex
- [`PlotNPCStackedBar()`](https://phytoecia.github.io/eCOMET/reference/PlotNPCStackedBar.md)
  : PlotNPCStackedBar
- [`BootCumulRichnessAUC()`](https://phytoecia.github.io/eCOMET/reference/BootCumulRichnessAUC.md)
  : BootCumulRichnessAUC
- [`BootstrapCumulativeRichness()`](https://phytoecia.github.io/eCOMET/reference/BootstrapCumulativeRichness.md)
  : BootstrapCumulativeRichness
- [`CalcNormalizedAUC()`](https://phytoecia.github.io/eCOMET/reference/CalcNormalizedAUC.md)
  : CalcNormalizedAUC
- [`CalculateCumulativeRichness()`](https://phytoecia.github.io/eCOMET/reference/CalculateCumulativeRichness.md)
  : CalculateCumulativeRichness
- [`CalculateNullCumulativeRichness()`](https://phytoecia.github.io/eCOMET/reference/CalculateNullCumulativeRichness.md)
  : CalculateNullCumulativeRichness

## Metadata regression analysis

- [`FeaturePhenotypeCorrelation()`](https://phytoecia.github.io/eCOMET/reference/FeaturePhenotypeCorrelation.md)
  : FeaturePhenotypeCorrelation
- [`ScreenFeaturePhenotypeCorrelation()`](https://phytoecia.github.io/eCOMET/reference/ScreenFeaturePhenotypeCorrelation.md)
  : Screen feature-phenotype correlation
- [`GetPerformanceFeatureCorrelation()`](https://phytoecia.github.io/eCOMET/reference/GetPerformanceFeatureCorrelation.md)
  : GetPerformanceFeatureCorrelation
- [`GetPerformanceFeatureLMM()`](https://phytoecia.github.io/eCOMET/reference/GetPerformanceFeatureLMM.md)
  : GetPerformanceFeatureLMM
- [`GetPerformanceFeatureRegression()`](https://phytoecia.github.io/eCOMET/reference/GetPerformanceFeatureRegression.md)
  : GetPerformanceFeatureRegression
- [`PlotFoldchangeResistanceQuad()`](https://phytoecia.github.io/eCOMET/reference/PlotFoldchangeResistanceQuad.md)
  : PlotFoldchangeResistanceQuad
- [`PlotFoldchangeResistanceRegression()`](https://phytoecia.github.io/eCOMET/reference/PlotFoldchangeResistanceRegression.md)
  : PlotFoldchangeResistanceRegression
- [`PlotFoldchangeResistanceRegression_t()`](https://phytoecia.github.io/eCOMET/reference/PlotFoldchangeResistanceRegression_t.md)
  : PlotFoldchangeResistanceRegression_t

## Output

- [`ExportFeaturesToCSV()`](https://phytoecia.github.io/eCOMET/reference/ExportFeaturesToCSV.md)
  : ExportFeaturesToCSV
- [`SaveMMO()`](https://phytoecia.github.io/eCOMET/reference/SaveMMO.md)
  : Save entire mmo object to a file (RDS)
- [`LoadMMO()`](https://phytoecia.github.io/eCOMET/reference/LoadMMO.md)
  : Load an mmo object previously saved with SaveMMO

## misc

- [`FeatureToID()`](https://phytoecia.github.io/eCOMET/reference/FeatureToID.md)
  : Convert feature names to IDs in the mmo object

- [`IDToFeature()`](https://phytoecia.github.io/eCOMET/reference/IDToFeature.md)
  : Convert feature IDs to names in the mmo object

- [`anova_tukey_dunnett()`](https://phytoecia.github.io/eCOMET/reference/anova_tukey_dunnett.md)
  : Perform ANOVA and Tukey's HSD test on the mmo object

- [`write_anova()`](https://phytoecia.github.io/eCOMET/reference/write_anova.md)
  : Write results of anova_tukey_dunnett to a CSV file

- [`permanova_stat()`](https://phytoecia.github.io/eCOMET/reference/permanova_stat.md)
  : Perform PERMANOVA and pairwise comparisons

- [`GetDistanceMat()`](https://phytoecia.github.io/eCOMET/reference/GetDistanceMat.md)
  : Get the distance matrix from the mmo object based on the specified
  distance metric

- [`GetNormFeature()`](https://phytoecia.github.io/eCOMET/reference/GetNormFeature.md)
  : Retrieve feature data from the mmo object, with normalization
  options

- [`filter_mmo()`](https://phytoecia.github.io/eCOMET/reference/filter_mmo.md)
  : Filter an mmo object by samples, groups, and/or features

- [`print(`*`<mmo>`*`)`](https://phytoecia.github.io/eCOMET/reference/print.mmo.md)
  :

  Print method for mmo objects Provides a clean, human-readable overview
  of an `mmo` list object instead of dumping the entire list when the
  object is printed in the console.
