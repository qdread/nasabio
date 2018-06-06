# nasabio

Readme documentation last updated by QDR on 06 June 2018

## Description of repo

This repository is dedicated to code for data processing and analysis for the NASA Biodiversity working group. Currently, almost all the scripts here are written in R, with a few shell scripts to run R code remotely on MSU's cluster. However most of the shell scripts are written on an as-needed basis directly on the cluster and so are not copied to the repository. Almost all the scripts here are specific to BBS and FIA datasets and would need to be edited if we decide to use other datasets.

There is a more detailed description of the code pipeline for calculating biodiversity and geodiversity metrics. It's found on the repo at `methods/r_pipeline.md`.

## Other locations referred to by the repo

- Most of the data and shell scripts are found on MSU's cluster at `/mnt/research/nasabio/`. That space has a `data` folder and a `code` folder. The `data` folder is organized by dataset and it should be easy to find things there. The `code` folder mostly has shell scripts and R scripts that were copied from this repository to the cluster so that code could be run remotely. 
- Most of the scripts that create figures or pdfs write them to the `NASABiodiversityWG/Figures/` folder on Google Drive, or one of its subfolders.
- Some of the scripts refer to Quentin's dropbox because it was convenient to copy some files there to do some exploratory analysis in R. All those files should also be on MSU's cluster.

## How the repo is organized

The repo now contains a lot of subdirectories roughly sorted by task:

- **alphadiversity**: Some preliminary scripts needed to get distance matrices for functional and phylogenetic diversity for birds. These are no longer used in the most recent version of the workflow.
- **betadiversity**: Some scripts to calculate beta diversity. The most important scripts here are `beta_part.r`, `beta_part_finalindex.r`, and `pairwise_beta_focal.r`. They all have the actual functions used to calculate beta-diversity. There are several scripts because there is still uncertainty about which methods will be used in the final analysis.
- **figs**: Scripts for plotting figures and maps.
- **imputation**: Scripts related to the trait imputation project.
- **methods**: RMarkdown and Markdown documents, and a few R scripts, with demonstrations of the different methods we are deciding between, as well as documents describing our methods in beautiful prose.
- **occupancy**: Scripts related to the BBS occupancy models.
- **prep_diversity_files**: This folder contains scripts that are used in the diversity workflow. Given the size of the data, there is a lot of preprocessing that needs to be done before actually calculating diversity.
- **readmes**: This contains some markdown files with technical documentation about the code.
- **run_compile_diversity**: This is also part of the diversity workflow, which includes scripts that calculate diversity metrics on the MSU cluster then compile the outputs into CSV files and export them.
- **spatial_data_extraction**: This folder contains scripts related to extracting data out of the rasters with geodiversity variables and calculating different summary metrics on them.
- **specieslists**: CSV files with species lists needed for analysis.
- **stats**: Scripts that actually attempt to do some statistical analysis on all the crunched numbers. This is also where the code is stored to do stratified random sampling on the plots for model fitting.
- **trait_phylo_data_processing**: Scripts that do some preprocessing on data sources like TRY and phylogenies that we got from external sources to calculate functional and phylogenetic diversity.
- **OLD**: This folder contains old and no longer used code; it's probably not necessary to refer to it.
 