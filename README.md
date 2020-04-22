# Deep Latent Microbiome

Microbial composition changes are associated with disease and performance/productivity. Discovering relevant microbiome pattern using data-driven analysis requires thousands of metagenomic samples.  Unfortunately, obtaining large numbers of metagenomic samples is expensive, and their degree of complexity limits classical association analyses. Here we propose a novel approach to microbiome representation based on Deep Learning using heterogenous autoencoders to condense the microbiome OTU vector into a latent space.  We show that encoded microbiomes can be expanded into an accurate reconstruction of the original species distribution.  Moreover, through the latent space, microbial species composition can be predicted with relatively high accuracy from the environmental features alone. This makes it possible to predict the structure of a microbiome under hypothetical scenarios, such as future climate states. We also show that transfer learning can be used to predict microbiomes in scenarios where sequencing resources are limited.

This *Deep Latent Microbiome* repository contains software (including Jupyter notebooks), data and results files that reproduce the study described in:

**Beatriz García-Jiménez, Jorge Muñoz, Sara Cabello, Joaquín Medina and Mark D. Wilkinson ; Predicting microbiomes through a deep latent space (Under review) **
 
 
***

The notebooks could be divided into:

* autoencoder_results: summary of results of experiments for selection of hyperparameters

* reference_model_predictions_and_analysis: build selected reference model (experiment 351), run analysis about microbiome reconstruction and prediction in test set, and prediction in novel ecosystems

* model_reference_\*domainFeatures\*: relevant feature analysis

* model_aggregated_\*: microbiome prediction at different taxonomic levels

* model_transferLearning: transfer learning analysis