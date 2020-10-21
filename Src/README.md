Code source files:

* `*.py`: python libraries with the functions called in the [Notebooks/](../Notebooks/). They code all the functionalities of the Deep Latent Microbiome system.

* `*.r`: R scripts that mainly analyses the results of the notebooks, and generates the figures for the manuscript and supplementary material. There is an independent script for each section or figure. In addition, `aggregation_taxomomic_levels.r` and `parsing_Maarastawi2018.r` pre-process input data files, for OTU aggregation at different taxonomic levels and transfer learning dataset, respectively.