# Welcome to eCOMET: A Tool for Mass Spectrometry–Based Ecological Metabolomics

<!-- badges: start -->
[![CRAN Status](https://www.r-pkg.org/badges/version/ecomet)](https://CRAN.R-project.org/package=ecomet)
[![License: AGPL-3](https://img.shields.io/badge/License-AGPL--3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![GitHub issues](https://img.shields.io/github/issues/phytoecia/eCOMET)](https://github.com/phytoecia/eCOMET/issues)
[![GitHub stars](https://img.shields.io/github/stars/phytoecia/eCOMET)](https://github.com/phytoecia/eCOMET/stargazers)
<!-- badges: end -->

**eCOMET** is an R package for processing and analyzing mass spectrometry–based metabolomics data in ecological and evolutionary contexts. It establishes a standardized pipeline that integrates MS1 feature–abundance tables with MS2 spectral similarity to generate common data products, including principal component analyses (PCA), chemical dendrograms where tips represent compounds, principal coordinates analyses (PCoA), and differential accumulation analyses. Unlike existing metabolomics toolkits that focus primarily on biomedical or cheminformatics applications, ecomet is designed explicitly for ecometabolomics. It emphasizes workflows that link metabolomic variation to ecological data, making it easier to move from raw mass spectrometry files to reproducible, comparative analyses of plant chemical diversity.

---

## Key Features

- Import and normalize MS1 feature–abundance tables
- Integrate MS2 spectral similarity for chemical diversity analysis
- Perform differential accumulation analysis to find DAMs (Differentially accumulated metabolites)
- Visualize results with PCA, PLS-DA, PCoA, NMDS, volcano plots, heatmaps, and more
- Quantify chemical diversity using alpha and beta diversity metrics, utilizing spectral similarity based on modified cosine score, MS2DeepScore, and DreaMS
- Chemical class enrichment analysis using CANOPUS annotation
- Screening molecular features related to sample metadata (e.g., temperature, height, herbivory)

<img width="1921" height="1080" alt="Screenshot 2025-07-24 at 6 59 47 PM" src="https://github.com/user-attachments/assets/517b33f9-7b66-4067-bda8-a66ae0a1d99c" />
---

## Installation
The eCOMET uses pairwiseAdonis for PERMANOVA; as it is not in CRAN, users may have to install it separately:
https://github.com/pmartinezarbizu/pairwiseAdonis

Install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("phytoecia/eCOMET")
library(ecomet)
```
## Preparing data for analysis
We recommend reading [Articles-Intro](./articles/Intro.html) before starting, to understand how the raw data from mass spectrometers can be prepared for analyses.  
We currently support data processed by MZMine (version 4 or later).

## Tutorials
eCOMET can be used for chemical ecological studies from diverse spectra of interests. We prepared two case studies:
1. [Treatment-based_study](./articles/Treatment-based_study_tutorial.html)
2. [Observational_field_study](./articles/Interspecific_study_tutorial.html)

---

For more details, see the documentation and vignettes.  
If you have questions or suggestions, visit the [GitHub repository](https://github.com/phytoecia/eCOMET), where we are actively developing. 
