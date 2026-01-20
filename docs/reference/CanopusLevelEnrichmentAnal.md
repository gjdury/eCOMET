# Enrichment analysis for Canopus-predicted terms

This function performs enrichment analysis for Canopus-predicted terms
on a given list of features.

## Usage

``` r
CanopusLevelEnrichmentAnal(
  mmo,
  list_test,
  pthr = 0.1,
  sig = TRUE,
  term_level = "NPC_pathway",
  representation = "greater",
  pval = "pval"
)
```

## Arguments

- mmo:

  The mmo object with sirius annotation and normalized data

- list_test:

  A vector containing names of features to analyze

- pthr:

  The threshold for adjusted p-value to be considered significant
  (default: 0.1)

- sig:

  A logical value indicating whether to return only significant terms
  (default: TRUE)

- term_level:

  The level of term to use for enrichment analysis Options are
  'NPC_pathway', 'NPC_superclass', 'NPC_class', 'ClassyFire_superclass',
  'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', or
  'ClassyFire_most_specific' (default: 'NPC_pathway')

- representation:

  The representation type for enrichment analysis. Options are 'greater'
  for overrepresentation (default: 'greater')

- pval:

  pvalue options-pval or fdr (default: 'pval')

## Value

A data frame containing the enrichment results, including term level,
term name, subset count, total count, fold enrichment, p-value, and
adjusted p-value (FDR)

## Examples

``` r
if (FALSE) {
# Perform enrichment analysis for a list of features using NPC_pathway level
sig_terms <- CanopusLevelEnrichmentAnal(
 mmo, list_test = c("feature1", "feature2"), pthr = 0.1,
 sig = TRUE, term_level = 'NPC_pathway', representation = 'greater'
)
# Perform enrichment analysis for a list of features using ClassyFire_class level and return all terms
all_terms <- CanopusLevelEnrichmentAnal(
 mmo, list_test = c("feature1", "feature2"), pthr = 0.1,
 sig = FALSE, term_level = 'ClassyFire_class', representation = 'greater'
)
}
```
