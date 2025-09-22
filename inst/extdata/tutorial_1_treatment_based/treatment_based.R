# pak::pak("phytoecia/eCOMET")
library(ecomet)

demo_feature <- system.file("extdata/tutorial_1_treatment_based/raw_data/250922_mzwiz_full_feature_table.csv", package = "ecomet")
demo_metadata <- system.file("extdata/tutorial_1_treatment_based/raw_data/250922_demo_metadata.csv", package = "ecomet")
demo_sirius_formula <- system.file("extdata/tutorial_1_treatment_based/raw_data/canopus_formula_summary.tsv", package = "ecomet")
demo_sirius_structure <- system.file("extdata/tutorial_1_treatment_based/raw_data/structure_identificaations.tsv", package = "ecomet")
demo_dreams <- system.file("extdata/tutorial_1_treatment_based/raw_data/250922_mzwiz_edges_dreams.csv", package = "ecomet")

mmo <- GetMZmineFeature(mzmine_dir=demo_feature, metadata_dir = demo_metadata, group_col = 'group')
