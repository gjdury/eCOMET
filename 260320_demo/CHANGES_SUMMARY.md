# Treatment-based Study Tutorial V2 — Changes Summary

**Date:** 2026-03-30
**Original file:** `vignettes/Treatment-based_study_tutorial.Rmd`
**New file:** `vignettes/Treatment-based_study_tutorial_V2.Rmd`

---

## 1. Developer comments removed

All internal developer TODO comments and discussion notes were removed:

- `#MS - I discussed this with Sedio Group...` (line 109)
- `#MS-same don't replace add...` (line 118)
- `#Can you add a little more annotation...` (line 125-128)
- `MS-I am still getting duplicates...` (line 141)
- `### MS- There are still duplicates...` (line 153)
- `#MS - I think that we should explain the logic...` (line 262)
- `#Need to explain how to handle duplicates!` (line 149)

These were in-progress discussion notes not suitable for a public tutorial.

---

## 2. Structure reorganized

| V1 Section | V2 Section | Change |
|---|---|---|
| Unnamed "Step 1" paragraph | **Step 1: Locate tutorial data** | Proper heading with explanation |
| "1. Create an mmo object" | **Step 2: Create an mmo object** | Expanded with sub-sections for each mmo component |
| Scattered normalization chunks | **Step 3: Preprocessing and normalization** | Grouped with table of normalization methods |
| SIRIUS/custom/dreams annotation scattered | **Step 4: Add annotations** (4.1–4.6) | Unified section with 6 clear sub-steps |
| "2. Plot dimensionality reduction" | **Step 5: Dimensionality reduction** | Split into PCA (5.1) and PLS-DA (5.2) |
| "3. Identify DAMs" | **Step 6: Identify DAMs** | Cleaner explanation with convention note |
| Volcano/Venn/Upset interleaved | **Step 7: Visualization of DAMs** | Dedicated section |
| "4. Heatmap" | **Step 8: Heatmap** | Three sub-sections (FC, Z-score, targeted GLS) |
| "5. CANOPUS enrichment" | **Step 9: CANOPUS class enrichment** | Single-list and multi-comparison sub-sections |
| "6. Correlation" | **Step 10: Correlation analysis** | Same content, better framing |
| "7. Sharing MMO" | **Step 11: Save and share** | Same content |
| (none) | **Summary table** | New: overview table of all steps/functions |

---

## 3. Explanatory content added

### Background section
- Added table describing the three experimental groups (ctrl, sl1, le1) with replicates
- Explained what a treatment-based study is and its goals

### mmo object structure (Step 2)
- Explained each component (`feature_data`, `feature_info`, `metadata`, `pairwise`) individually with context

### Preprocessing (Step 3)
- Table summarizing all normalization methods, their result location, and purpose
- Explained *why* zero replacement is needed before log/fold-change calculations
- Described `MassNormalization()` formula

### Annotations (Step 4)
- Added note on SIRIUS duplicate handling
- Included guidance on choosing CANOPUS threshold with link to SIRIUS documentation
- Explained COSMIC confidence score interpretation with reference to the Nature paper
- Described the `filter_id` / `id_list` API pattern for feature filtering

### Heatmap (Step 8)
- Added options table for `GenerateHeatmapInputs()`
- **Bug fix:** Added Inf/NaN handling for fold-change matrices (log2(x/0) = Inf when control has zero abundance)
- Documented the `pheatmap(filename=...)` pattern to avoid PDF contamination from `GenerateHeatmapInputs(distance=...)` side effects

### Enrichment (Step 9)
- Added table of all available classification levels

---

## 4. Code fixes and updates

| Issue | V1 | V2 |
|---|---|---|
| Feature filter API | `filter_feature = TRUE, feature_list = FLVs` (deprecated) | `filter_id = TRUE, id_list = FLVs` (current API) |
| Variable name inconsistency | `FLVs_features` used but defined as `FLVs` later | Consistent `FLVs` throughout |
| Heatmap Inf crash | No handling | Added `Inf/NaN` capping to `finite_max` before `pheatmap()` |
| Heatmap PDF pattern | `pdf()/dev.off()` wrapper | `pheatmap(filename = ...)` to avoid side-effect contamination |
| Feature vector extraction | `pull(feature)` | `pull(id)` (consistent with mmo ID-based API) |
| Output paths | `'plot/PCA'` (relative, no demo prefix) | `'260320_demo/plot/PCA/PCA_basic'` (organized) |
| FeaturePresence comment | `# Add Zscore` (wrong comment) | Correct description |

---

## 5. Style improvements

- Consistent use of `eval=FALSE` in all code chunks
- Markdown tables instead of inline lists for structured information
- Blockquote callouts (`> **Note:**`, `> **Important:**`) for warnings and tips
- External links to SIRIUS documentation formatted as proper hyperlinks
- All output paths organized under `260320_demo/plot/{type}/` structure

---

## 6. Output files generated

All outputs verified successfully in `260320_demo/`:

```
260320_demo/
  plot/
    PCA/          — 3 PCA plots (basic, log, FLV) + custom + PERMANOVA CSVs
    PLSDA/        — 2 PLS-DA plots (basic, meancentered)
    Volcano/      — 2 volcano plots (sl1, le1) + data CSVs
    Venn/          — 1 Venn diagram (upregulated DAMs)
    Upset/         — 1 UpSet plot (upregulated DAMs)
    Heatmap/       — 3 heatmaps (FC+dreams, Z-score, GLS targeted)
    Enrichment/    — 5 enrichment plots (single-list x2, NPC class, all NPC, all ClassyFire)
  output/
    mmo.RData              — saved mmo object
    sl_correlation_results.csv — phenotype correlation results
```
