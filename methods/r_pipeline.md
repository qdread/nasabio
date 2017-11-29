# Pipeline for NASA bioxgeo

Last updated by QDR on 29 Nov 2017

This is documentation explaining which R scripts should be used to extract biodiversity and geodiversity information from different data sources in the NASA bioxgeo project, and calculate summary metrics on them. In many cases, the R scripts mentioned will only calculate the metrics for a very small slice of the data. To get the metrics for all the data, you need to run the script many times in parallel on the MSU cluster, one for each slice. Each time, it will save an .R object or .CSV file. The files for each slice then need to be combined into a single output .CSV in yet another script.

The names of the scripts are prefixed with the subdirectory on the Git repo where they are saved. 

## 1. Biodiversity

### 1.0 Define functions

There are two methods currently being used to calculate beta-diversity. The newer method uses the measure based on the multisite S&oslash;rensen dissimilarity index. In addition, the newer method uses Baselga's formulas to partition taxonomic and phylogenetic beta-diversity into a portion due to nestedness and a portion due to species turnover. Unfortunately, this method does not support functional diversity, only TD and PD. Because of this, the old method is still part of the workflow. In the old method, we calculate pairwise beta-diversity among all sites within the radius for TD, PD, and FD, then average them to get a value for the circle. 

- `betadiversity/pairwise_beta_focal.r`: This script has the function definitions for the old method beta-diversity. The function currently used is `singlepair_beta()` which calculates pairwise dissimilarity between any two sites (TD, PD, and FD). This script also contains a function called `diversity_3ways()` which calculates either alpha or gamma diversity (TD, PD, and FD) for a single site (alpha) or a site-by-species matrix (gamma).
- `betadiversity/beta_part_finalindex.r`: This script has the function definition for new method beta-diversity. It uses all the indices that we agreed on.
- `betadiversity/beta_part.r`: This script has a function that will calculate beta-diversity partitioning using a lot of different methods, which was used to compare the results from the different methods so that we could agree upon one to use.

### 1.1 BBS

All BBS diversity metrics are *incidence-only*. We do not use any information about abundance. Originally, a different value of diversity was found for each year then averaged across the time period. Now, incidences are pooled across ten years of data for each site to get a single community composition value for the site across time. This is a very conservative method that means it is not necessary to run the models correcting for imperfect detection. The old and new scripts are all in the Git repo. There are also some even older scripts I wrote originally that do separate diversity calculations for individual stops rather than at the route level, but I think we decided not to pursue that. Here I have only listed the names of the most up-to-date scripts that find diversity at the route level pooled across years.

#### 1.1.1 Prepare matrices

- `prep_diversity_files/bbsbeta_byroute_prep.r`: This script loads in a site-by-species matrix for BBS (this matrix was generated using code written for AquaXTerra project). It throws out nocturnal birds which are not analyzed, then calculates distances among plots to determine which plots will fall in the same neighborhoods. At the end of the script the R workspace is saved to `/mnt/research/nasabio/data/bbs`.
- `prep_diversity_files/bbs_oneyear_getidx.r`: This script is run in parallel on the cluster. For each site by radius combination, it generates a site-by-species matrix containing the focal route and all its neighboring BBS routes in the radius. (The older version of this script that is done for each year separately is `prep_diversity_files/bbs_getidx.r`.) The output files are written to `/mnt/research/nasabio/data/bbs/mats`.

#### 1.1.2 Find alpha-diversity

- `run_compile_diversity/bbs1year/bbs_allplotsalpha.r`: Calculates alpha diversity (TD, PD, and FD) for each BBS route. Run in parallel on cluster. The temporary output files are written to `/mnt/research/nasabio/data/bbs/diversity`.

#### 1.1.3 Find gamma-diversity

- `run_compile_diversity/bbs1year/bbs_allradiigamma.r`: Calculates gamma diversity (TD, PD, and FD) for each BBS route for each radius. Run in parallel on cluster. The temporary output files are written to `/mnt/research/nasabio/data/bbs/diversity`.

#### 1.1.4 Find beta-diversity

- `run_compile_diversity/bbs1year/bbs_allplotsbetapart.r`: New method (multisite index for each route and each radius, does not include FD). The temporary output files are written to `/mnt/research/nasabio/data/bbs/diversity`.
- `run_compile_diversity/bbs1year/bbs_allpairsbeta.r`: Old method (pairwise indices for each pair of routes, includes FD). The temporary output files are written to `/mnt/research/nasabio/data/bbs/diversity`.

#### 1.1.5 Combine slices and aggregate by radius

- `run_compile_diversity/bbs1year/compilebbs1year.r`: Combines the many output files from the alpha, beta, and gamma runs and exports single CSV files to `/mnt/research/nasabio/data/bbs`.

### 1.2 FIA

We calculate FIA diversity metrics for both incidence and abundance. We can do this because our estimate of abundance is more robust than for BBS. In all cases, basal area of a species represents its abundance (counting stems is a poor measure of abundance).

#### 1.2.1 Prepare matrices

- `prep_diversity_files/fiabeta_prep.r`: This script loads the raw FIA tree data and generates the site-by-species matrix. It also loads phylogenetic and functional information. The resulting workspace is saved to `/mnt/research/nasabio/data/fia`. Because of security, the FIA true coordinates are not saved in this workspace. Only the neighbor identities are needed to calculate diversity, not the true plot locations.
- `prep_diversity_files/fia_getidx.r`: This script is run in parallel on the cluster. For each site by radius combination, it generates a site-by-species matrix containing the focal plot and all its neighboring BBS plots in the radius. The output files are written to `/mnt/research/nasabio/data/fia/mats`.
- `prep_diversity_files/fiaprep_compilemats.r`: This script combines the slices in the previous script into single files.

#### 1.2.2 Find alpha-diversity

- `run_compile_diversity/fia_allplotsalpha.r`: Calculates alpha diversity (TD, PD, and FD) for each FIA plot. Run in parallel on cluster. The temporary output files are written to `/mnt/research/nasabio/data/fia/diversity`.

#### 1.2.3 Find gamma-diversity

- `run_compile_diversity/run_compile_diversity/fia_allradiigamma.r`: Calculates gamma diversity (TD, PD, and FD) for each FIA plot for each radius. Run in parallel on cluster. The temporary output files are written to `/mnt/research/nasabio/data/fia/diversity`.

#### 1.2.4 Find beta-diversity

- `run_compile_diversity/fia_allpairsbeta.r`: New method (multisite index for each plot and each radius, does not include FD). The temporary output files are written to `/mnt/research/nasabio/data/fia/diversity`.
- `run_compile_diversity/fia_betapart_smallslice.r`: Old method (pairwise indices for each pair of plots, includes FD). The temporary output files are written to `/mnt/research/nasabio/data/fia/diversity`.

#### 1.2.5 Combine slices and aggregate by radius

Each script below exports final .csv files to `/mnt/research/nasabio/data/fia`.

- `run_compile_diversity/fia_alpharadius.r`: Combines alpha-diversity slices and aggregates by radius.
- `run_compile_diversity/fia_compilegamma.r`: Combines gamma-diversity slices.
- `run_compile_diversity/fia_compilebetapart.r`: Combines new-method beta-diversity slices.
- `run_compile_diversity/fiabeta_radius.r`: Combines old-method beta-diversity slices (this script may need to be updated for the unfuzzed plots).

## 2. Geodiversity

For now this is just the mean, standard deviation, minimum, and maximum for continuous variables, and the richness (number of unique values) and diversity (Shannon entropy) for categorical variables.

### 2.0 Define functions

- `spatial_data_extraction/extractbox.r`: Defines function `extractBox()` to extract a square raster from the large latitude-longitude rasters that include either the entire contiguous USA or the entire world. Next defines function `statsByRadius()` for continuous variables and `diversityByRadius()` for categorical variables treated as non-ordinal. The latter functions load the entire cropped square raster into RAM, take subsets of the various circle sizes you input, and calculate summary statistics on them. **Note**: I edited these on 29 November. Previously, these functions used pre-calculated distance tables which gave erroneous results. Any geodiversity statistics from before 29 November 2017 are wrong.

### 2.1 Extract pixels from different radii and calculate summary metrics

There is a separate R script for each raster that was extracted. For the most part I used the command line to create virtual rasters (.vrt) from all the .tif files, which makes it easier to read them in parallel. They are all in the `spatial_data_extraction` folder. Temporary output files are written to `climstats`, `elevstats`, and `geostats` subdirectories in either the `data/bbs` or `data/fia` directories on the cluster.

### 2.2 Combine slices into single output files

- `spatial_data_extraction/compileelev_bioclim.r`: Combines all extracted geodiversity metrics for both BBS and FIA and exports them to single .csv files, one for BBS and one for FIA. The master geodiversity .csv file is written to the `data/bbs` or `data/fia` directory on the cluster.
