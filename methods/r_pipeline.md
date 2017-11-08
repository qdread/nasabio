# Pipeline for NASA bioxgeo

Last updated by QDR on 08 Nov 2017 (still incomplete but work in progress . . . )

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

- `prep_diversity_files/bbsbeta_byroute_prep.r`: This script loads in a site-by-species matrix for BBS (this matrix was generated using code written for AquaXTerra project). It throws out nocturnal birds which are not analyzed, then calculates distances among plots to determine which plots will fall in the same neighborhoods. 
- `prep_diversity_files/bbs_oneyear_getidx.r`: This script is run in parallel on the cluster. For each site by radius combination, it generates a site-by-species matrix containing the focal route and all its neighboring BBS routes in the radius. (The older version of this script that is done for each year separately is `prep_diversity_files/bbs_getidx.r`.)

#### 1.1.2 Find alpha-diversity

- `run_compile_diversity/bbs1year/bbs_allplotsalpha.r`

#### 1.1.3 Find gamma-diversity

- `run_compile_diversity/bbs1year/bbs_allradiigamma.r`

#### 1.1.4 Find beta-diversity

- `run_compile_diversity/bbs1year/bbs_allplotsbetapart.r`: New method
- `run_compile_diversity/bbs1year/bbs_allpairsbeta.r`: Old method

#### 1.1.5 Combine slices and aggregate by radius

- `run_compile_diversity/bbs1year/compilebbs1year.r`

### 1.2 FIA

We calculate FIA diversity metrics for both incidence and abundance. Our estimate of abundance is more robust than for BBS. In all cases, basal area of a species represents its abundance (counting stems is a poor measure of abundance).

#### 1.2.1 Prepare matrices

- `prep_diversity_files/fiabeta_prep.r`
- `prep_diversity_files/fia_getidx.r`

#### 1.2.2 Find alpha-diversity

- `run_compile_diversity/fia_allplotsalpha.r`

#### 1.2.3 Find gamma-diversity

- `run_compile_diversity/run_compile_diversity/fia_allradiigamma.r`

#### 1.2.4 Find beta-diversity

- `run_compile_diversity/fia_allpairsbeta.r`: New method
- `run_compile_diversity/fia_allplotsbetapart.r`: Old method

#### 1.2.5 Combine slices and aggregate by radius

- `fia_alpharadius.r`
- `fiabeta_radius.r`

## 2. Geodiversity

For now this is just the mean, standard deviation, minimum, and maximum for continuous variables, and the richness (number of unique values) and diversity (Shannon entropy) for categorical variables.

### 2.0 Define functions

- `spatial_data_extraction/precalcdist.r`: Defines functions to precalculate distances and create square matrices that are TRUE where a point falls within a circle centered at the center of the matrix, and FALSE otherwise, then convert those matrices into tables for many different radii. This is used to quickly extract all the pixels within a circle of a given radius without having to run those calculations each time.

### 2.1 Extract pixels from different radii and calculate summary metrics

### 2.2 Combine slices into single output files

- `spatial_data_extraction/compileelev_bioclim.r`: Combines all extracted geodiversity metrics for both BBS and FIA and exports them to single .CSV files, one for BBS and one for FIA.
