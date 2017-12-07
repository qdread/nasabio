# Methods: calculating taxonomic, functional, and phylogenetic diversity for BBS and FIA

Author: QDR  
Project: NASA Biodiversity (and AquaXTerra)
Date created: 05 March 2017
Last modified: 04 October 2017

**Edits on 04 October 2017**: added percent ambiguous individuals, edited occupancy modeling description, added consensus tree based on 1000 trees, edit description of FIA traits, say we pooled BBS occurrence across years, also made changes to the diversity methods updating them with the newer methods. This reflects our consensus on what the best methods are.

## 1. Breeding Bird Survey: preparation of dataset

### 1.1 Survey Data

We obtained the most recent USGS Breeding Bird Survey (BBS) data from https://www.pwrc.usgs.gov/bbs/. Breeding Bird Survey data consists of point counts taken every 800m along prescribed routes. Most of the ~5000 routes are surveyed once yearly at the time of peak breeding bird abundance. See *references* for more detailed description of the survey protocols. We located each route at the centroid of the points defining the transect.

**Quality control**: The BBS protocol is known to have poor success detecting nocturnal species; we excluded all nocturnal species from any diversity calculations (nocturnality was defined based on information taken from the traits databases described below). Because the number of observations of each species at each BBS point count is not a reliable index of abundance, all the functional, taxonomic, and phylogenetic diversity indices for breeding birds are presence-only metrics. Furthermore, not all BBS observations represent individuals that can be unambiguously identified as a single bird species. We dealt with this issue using the following rules, separately for each unidentified individual: 
- If the individual was coded as a subspecies, we assigned it to the parent species.
- If the individual was coded as a hybrid of two species, we assigned it randomly to one of the two parent species.
- If the individual was coded as possibly belonging to two or more species, we assigned it randomly to one of those species.
- If the individual was only identified to the genus or family level, we assigned it randomly to a species in that genus or family.
The ambiguous individuals only comprised a small percentage of the total (about 1.8%), so this should have a relatively small effect on the values of functional and phylogenetic diversity we calculated. Any other discrepancies in species names were resolved so that the most recent taxonomy is used.

#### 1.1.1 Occupancy Modeling

As of 6 July, QDR is working on implementing the occupancy model following the model described in Jarzyna & Jetz 2016, and using the JAGS model specification provided by Jarzyna. It is possible that any inference we may make about variation in bird diversity in space and time is biased because bird species have a significant probability of escaping detection at every point along the 50-stop transect. What is more, the probability that a species will be missed is not random with respect to habitat type (probably higher in routes with denser vegetation), phylogeny (probably higher in taxa with a lot of small, inconspicuous species), and functional traits (probably many of the traits we use to calculate FD are confounded with being small and inconspicuous to surveyors). The model assumes that there is a latent variable which represents the true occupancy of a species at a route. It is possible that the species will be undetected in a subset of the sites where it is truly present. For example, if a bird species is only observed at 1 of the 50 stops at many of the routes where it is observed at all, it is likely that it went undetected at many routes where it was not observed. In contrast, if a bird species is observed at 20 or 30 of the stops at most of the routes where it is observed at all, it is likely that it probably isn't present if not observed. If we choose to implement this model, it would mean fitting the latent variable to the true occupancies and some covariates, and then using the modeled values to calculate diversity rather than the raw observations. 

As of 4 October, we have addressed this partially by pooling incidences from 1997-2016 in the BBS dataset to single values for each route. We are still going to explore this further if it is warranted, but the current consensus is to acknowledge it in the paper but not actually run the correction. The method of Jarzyna and Jetz is too problematic.

### 1.2 Trait Data

We obtained species mean trait values from the following two sources: the Amniote Life History Database (Myhrvold et al. 2015), which has morphological and life-history trait values for each of the world's bird species, and EltonTraits (Wilman et al. 2014), which has foraging traits that define the resource niche of each of the world's bird species in continuous space. For some of these species-level trait values, no actual observations are available, so trait values were imputed based on taxonomy; see the above references for details on how the database authors interpolated the missing values. We used the following traits: proportion of the diet consisting of invertebrates, birds/mammals, reptiles/amphibians, fish, scavenged meat, fruit, nectar, and seeds; proportion of foraging time spent in the water below and above the surf line, on the ground, in the understory, at mid-height, in the canopy, and in the air; pelagic status; time to maturity of females and males; clutch size; number of clutches per year; average adult body mass; average and maximum longevity; birth weight; egg mass; time of incubation; age at fledging.

### 1.3 Phylogenetic Data

We obtained a recent phylogeny of all bird species (Jetz et al. 2012). The phylogeny consists of a posterior distribution of 10000 trees; we randomly chose 1000 trees from this posterior distribution and generated a consensus tree with branch lengths using the R package *phytools*. We used the branch lengths from this consensus tree to calculate the distance-based phylogenetic diversity indices described below.

## 2. Forest Inventory and Analysis: preparation of dataset

### 2.1 Survey Data

We obtained USFS Forest Inventory and Analysis (FIA) data for the Pacific Northwest region. Forest plots surveyed according to the FIA protocol consist of four subplots, each circles with 7.3m radius, located 36.6 m from one another in a three-pointed star pattern. See *references* for more detailed description of the survey protocols. Each tree in each subplot is identified to species, and its diameter at breast height is recorded. Using the diameters to calculate basal area of each individual tree, we summed the basal areas within each species to estimate the relative abundance of each species in each subplot. Any discrepancies in species names were resolved so that the most recent taxonomy is used.

### 2.2 Trait Data

We obtained all available trait data for all species in our survey dataset from the TRY database (http://try-db.org), and supplemented these traits with data collected by Jens Stevens (unpublished). We restricted our analysis to continuousl traits that have a documented link to species performance and niche and that have at least one available observation for the majority of species in our dataset. These criteria applied to the following traits: specific leaf area, leaf C content per area, leaf N content per area, leaf N content per dry mass, leaf P content per dry mass, leaf C:N ratio, leaf N:P ratio, leaf lifespan, photosynthetic rate per leaf area, photosynthetic rate per leaf dry mass, leaf thickness, litter decomposition rate, plant lifespan, plant shade tolerance category, rooting depth, seed dry mass, proportion vessel area per unit stem cross-sectional area, specific stem density, and stomatal conductance per leaf area. We used phylogenetic imputation implemented in the R package *Rphylopars* to impute missing values. For the final analysis, we used the following six traits that had a value for almost all species to minimize the number of imputed traits: bark thickness, specific leaf area, specific stem density, seed dry mass, rooting depth, and plant lifespan.

Things to add to trait data methods when done: 
- we can expand the traits list with the added traits that John's undergrads are compiling.
- We have provisionally developed a hybrid phylogenetic-spatial imputation method for filling holes in the TRY dataset, based on Jay Jain's work.

### 2.3 Phylogenetic Data

We obtained a phylogeny of all trees in the FIA survey area. Details on how the phylogeny was assembled are available at (cite Kevin Potter's two pubs). For this phylogeny, only a single tree was available rather than a distribution of trees, so all phylogenetic diversity metrics that we calculated for FIA trees are based on this single tree assuming that it represents a consensus tree. 

## 3. Diversity calculations for both datasets

For both FIA and BBS, we calculated diversity metrics based on species presence at each plot or route, respectively. We also calculated abundance (basal area)-weighted diversity metrics for FIA. Because abundance estimates are not robust for BBS, we did not calculate abundance-weighted metrics. For FIA, we used the most recent survey as a single time point for each plot. For BBS, we pooled all incidence data from 1997-2016 into a single presence-absence value for each species at the route level and calculated diversity indices based on this pooled presence to address the issue of imperfect detection. If a species was present at any of the fifty point counts at any time during the ten years, it was treated as present.

We calculated alpha, beta, and gamma diversity at a number of different radii around each FIA plot and BBS route by taking the median diversity of all plots in the radius, including the focal plot (alpha), the mean pairwise diversity (with appropriate transformations if needed) of all pairs of plots in the radius, including the focal plot (beta), and the aggregated diversity of all plots in the radius as if they were a single community (gamma). For FIA plots, the radii were 1, 5, 7.5, 10, 20, 50, 75, 100, 150, 200, and 300 km. For BBS routes, because the average distance between route centroids was greater than the distance between FIA plots, the radii were 50, 75, 100, 150, 200, and 300 km.

### 3.1 Alpha and gamma diversity

We calculated taxonomic, functional, and phylogenetic alpha-diversity (diversity of a local community) for tree and bird communities. For the bird communities, we calculated diversity indices for communities aggregated at the route level (aggregating the 50 point counts along one route). For the tree communities, we calculated diversity indices for communities aggregated at the plot level (aggregating the four subplots making up one plot). For each plot/route and radius, we calculated alpha diversity within that radius by taking the median diversity value for all plots or routes (including the focal plot) located inside the circle defined by the radius around the focal plot/route. For gamma diversity, we aggregated all the plots/routes within the focal circle to a single community, and calculated taxonomic, functional, and phylogenetic diversity of that community.

#### 3.1.1 Taxonomic

We calculated species richness by totaling the number of species in each community. We did not calculate any abundance-weighted metrics of taxonomic diversity for BBS. However, for FIA, we also calculated basal-area-weighted Shannon diversity as follows: $`H' = -\sum_{i=1}^{R} p_i \ln p_i`$, or the natural logarithm of the true diversity with q = 1, where R is species richness and p<sub>i</sub> is the basal area of species i.

#### 3.1.2 Functional

**New way**

We calculated two distance-based metrics of functional diversity using the `comdist()` and `comdistnt()` functions in the R package *picante*. These metrics are based on the mean pairwise functional distance among all pairs of species in the community, and the mean nearest-neighbor distance for each species in the community, respectively. We chose to use these distance-based metrics for consistency, because our phylogenetic diversity metrics are also distance-based using mean and nearest-neighbor distances, and our beta-diversity metrics are also all distance-based.   

**Old way**

We calculated distance-based metrics of functional diversity using the dbFD() function in the R package FD (cite). These metrics, including FRic, FEve, FDiv, and FDis, are all calculated by generating pairwise distances between species in n-dimensional trait space, where n is the number of functional traits. The method (Gower distance) is able to interpolate missing values, provided that at least several traits are available for each species. We excluded all species that had less than 3 traits measured from the functional diversity calculations. Functional richness (FRic) describes the total volume occupied by a community in trait space, functional evenness (FEve) describes how evenly the community fills the volume of trait space that it occupies, functional divergence (FDiv) is the average distance of each species from the centroid of trait space, and functional dispersion (FDis) is similar but the location of the centroid is weighted by the relative abundance of species (cite Kuebbing et al., in prep, or references therein including Villeger et al.). We excluded communities with less than 3 unique species, as well as communities for which we had trait data for less than 75% of the individuals. We present abundance-weighted functional diversity metrics for FIA communities, but presence-only metrics for BBS communities. 

#### 3.1.3 Phylogenetic

We calculated two distance-based metrics of phylogenetic diversity using functions in the R package *picante*: mean pairwise phylogenetic distance (MPD) and mean nearest-taxon phylogenetic distance (MNTD). We conducted 99 randomizations of the phylogenetic distance matrix and calculated the z-score of the observed phylogenetic distances relative to the distribution of phylogenetic distances of the randomized matrices. For the FIA dataset, we calculated both abundance-weighted and presence-only versions of MPD and MNTD, but for the BBS dataset we only calculated the presence-only version. 

### 3.2 Beta diversity 

We calculated taxonomic, functional, and phylogenetic beta-diversity (turnover of diversity among local communities) for tree and bird communities. Beta-diversity is defined as the variation in community composition across multiple local communities. To determine beta-diversity at a point, it is necessary to define the kernel or radius within which variation in community composition is taken into account. For the FIA dataset, we aggregated species abundances of each plot and calculated beta-diversity for each plot at a number of different radii around the focal plot; as the radius increases, the number of pairwise comparisons among plots also increases as more plots fall within the kernel. 

Things to do here, if desired:
- Use rarefaction-based estimate of effective number of communities as the denominator in the multi-site mean pairwise dissimilarity metric, rather than the raw number of sites sampled

#### 3.2.1 Taxonomic

**Old way** 

We calculated taxonomic beta-diversity using the d() function from the R package vegetarian (diversity partitioning method). We also calculated beta-diversity with the pairwise dissimilarity method using the vegdist() function from the R package vegan and taking the mean of the pairwise Jaccard dissimilarity of all local communities within a particular radius of the focal plot. As before, we calculated both abundance-weighted and presence-only metrics for FIA communities, and only presence-only metrics for BBS communities. 

**New way**

We calculated taxonomic beta-diversity using the `betapart()` family of functions from the R package *betapart*. We opted to use the pairwise distance-based method rather than Whittaker's diversity-partitioning method because it allows us to use more consistent methods across taxonomic, functional, and phylogenetic diversity calculations. We calculated multi-site indices based on the Sorensen dissimilarities among communities. Each of these represent the mean pairwise dissimilarity among plots or routes within a given radius. We were additionally interested in the relative contributions of nestedness and turnover to beta-diversity at different spatial scales. There are two commonly used methods to partition the multi-site index (cite Baselga, cite Podani). We opted to use Baselga's method.  

#### 3.2.2 Functional

We calculated functional beta-diversity using the `comdist()` and `comdistnt()` function from the R package picante. We used the trait databases described above to generate pairwise distance matrices for all species in each survey. We imputed missing trait values as described above. We calculated both pairwise and nearest-taxon functional beta-diversity, both abundance-weighted and presence-only for FIA, and just presence-only for BBS. Refer to Swenson et al. 2011, PLoS One, for definitions of the metrics we calculated.

#### 3.2.3 Phylogenetic

**Old way**

We calculated phylogenetic beta-diversity with the same set of metrics as we used for functional beta-diversity (including pairwise and nearest-neighbor, and abundance-weighted metrics for FIA plots only). We used the phylogenies described above to generate pairwise distance matrices for all species in each survey.

**New way**

We calculated phylogenetic beta-diversity using the `betapart()` family of functions, in a similar manner as for taxonomic diversity. We also partitioned phylogenetic beta-diversity into nestedness and turnover components, following Baselga's method, as we did for taxonomic diversity.

