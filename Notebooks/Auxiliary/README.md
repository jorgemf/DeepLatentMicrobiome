The auxiliary notebooks are divided into:

* `model_aggregated_*`: microbiome prediction at different taxonomic levels (class, family, genus, order, phylum)

* `model_reference_*domainFeatures*`: relevant feature analysis (1, 2, 3 or all domain features)

* `model_transferLearning*`: transfer learning analysis. In a 100 samples subset from [Walters et al., 2018](https://doi.org/10.1073/pnas.1800918115) with 2, 3 or all features; or in [Maarastawi et al.](https://doi.org/10.3389/fmicb.2018.01295) dataset with differente features.

* `compute_performance_per_sample_test_set`: computation for supplementary figure.

* `reference_model_novel_predictions`: base notebook for the not-installation user-friendly interface https://tinyurl.com/DeepLatentMicrobiome-Maize