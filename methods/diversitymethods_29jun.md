# Methods: calculating taxonomic, functional, and phylogenetic diversity for BBS and FIA

Author: QDR  
Project: NASA Biodiversity (and AquaXTerra)
Date created: 05 March 2017
Last modified: 29 June 2017

## 1. Breeding Bird Survey: preparation of dataset

### 1.1 Survey Data

We obtained the most recent USGS Breeding Bird Survey (BBS) data from https://www.pwrc.usgs.gov/bbs/. Breeding Bird Survey data consists of point counts taken every 800m along prescribed routes. Most of the ~5000 routes are surveyed once yearly at the time of peak breeding bird abundance. See *references* for more detailed description of the survey protocols. We located each route at the centroid of the points defining the transect.

**Quality control**: The BBS protocol is known to have poor success detecting nocturnal species; we excluded all nocturnal species from any diversity calculations (nocturnality was defined based on information taken from the traits databases described below). Because the number of observations of each species at each BBS point count is not a reliable index of abundance, all the functional, taxonomic, and phylogenetic diversity indices for breeding birds are presence-only metrics. Furthermore, not all BBS observations represent individuals that can be unambiguously identified as a single bird species. We dealt with this issue using the following rules, separately for each unidentified individual: 
- If the individual was coded as a subspecies, we assigned it to the parent species.
- If the individual was coded as a hybrid of two species, we assigned it randomly to one of the two parent species.
- If the individual was coded as possibly belonging to two or more species, we assigned it randomly to one of those species.
- If the individual was only identified to the genus or family level, we assigned it randomly to a species in that genus or family.
The ambiguous individuals only comprised a small percentage of the total (*look this up*), so this should have a relatively small effect on the values of functional and phylogenetic diversity we calculated.

### 1.2 Trait Data

We obtained species mean trait values from the following two sources: the Amniote Life History Database (Myhrvold et al. 2015), which has morphological and life-history trait values for each of the world's bird species, and EltonTraits (Wilman et al. 2014), which has foraging traits that define the resource niche of each of the world's bird species in continuous space. For some of these species-level trait values, no actual observations are available, so trait values were interpolated based on taxonomy; see the above references for details on how the database authors interpolated the missing values. **add exactly what traits we ended up using**

### 1.3 Phylogenetic Data

We obtained a recent phylogeny of all bird species (Jetz et al. 2012). The phylogeny consists of a posterior distribution of 10000 trees; we randomly chose 10 trees from this posterior distribution and calculated all phylogenetic diversity indices separately for each of the 10 trees. We took the mean of each diversity index across the 10 trees for each point.

## 2. Forest Inventory and Analysis: preparation of dataset

### 2.1 Survey Data

We obtained USFS Forest Inventory and Analysis (FIA) data for the Pacific Northwest region. Forest plots surveyed according to the FIA protocol consist of four subplots, each circles with 7.3m radius, located 36.6 m from one another in a three-pointed star pattern. See *references* for more detailed description of the survey protocols. Each tree in each subplot is identified to species, and its diameter at breast height is recorded. Using the diameters to calculate basal area of each individual tree, we summed the basal areas within each species to estimate the relative abundance of each species in each subplot. 

### 2.2 Trait Data

We obtained all available trait data for all species in our survey dataset from the TRY database (http://try-db.org). We restricted our analysis to continuous traits that have a documented link to species performance and niche and that have at least one available observation for the majority of species in our dataset. These criteria applied to the following traits: **list of traits here**.  

### 2.3 Phylogenetic Data

We obtained a recent phylogeny of seed plants (Qian and Jin 2016). For this phylogeny, only a single tree was available, so all phylogenetic diversity metrics that we calculated for FIA trees are based on this single tree. Every genus was represented in the phylogeny. If a species present in FIA plots was not present in the phylogeny, we used the R package phytools (cite) to add it to the root of the appropriate genus, assuming that it is equally distantly related to all other species in that genus. 

## 3. Diversity calculations for both datasets

### 3.1 Alpha and gamma diversity

We calculated taxonomic, functional, and phylogenetic alpha-diversity (diversity of a local community) for tree and bird communities. For the bird communities, we calculated diversity indices for communities aggregated at the route level (aggregating the 50 point counts along one route). For the tree communities, we calculated diversity indices for communities aggregated at the plot level (aggregating the four subplots making up one plot). For each plot/route and radius, we calculated alpha diversity within that radius by taking the median diversity value for all plots or routes (including the focal plot) located inside the circle defined by the radius around the focal plot/route. For gamma diversity, we aggregated all the plots/routes within the focal circle to a single community, and calculated taxonomic, functional, and phylogenetic diversity of that community. 

#### 3.1.1 Taxonomic

We calculated species richness by totaling the number of species in each community. We did not calculate any abundance-weighted metrics of taxonomic diversity for BBS. However, for FIA, we also calculated Shannon diversity as follows: **insert equation here**.

#### 3.1.2 Functional

We calculated distance-based metrics of functional diversity using the dbFD() function in the R package FD (cite). These metrics, including FRic, FEve, FDiv, and FDis, are all calculated by generating pairwise distances between species in n-dimensional trait space, where n is the number of functional traits. The method is able to interpolate missing values, provided that at least several traits are available for each species. Functional richness (FRic) describes the total volume occupied by a community in trait space, functional evenness (FEve) describes how evenly the community fills the volume of trait space that it occupies, functional divergence (FDiv) is the average distance of each species from the centroid of trait space, and functional dispersion (FDis) is similar but the location of the centroid is weighted by the relative abundance of species (cite Kuebbing et al., in prep, or references therein including Villeger et al.). We excluded communities with less than 3 unique species, as well as communities for which we had trait data for less than 75% of the individuals. We present abundance-weighted functional diversity metrics for FIA communities, but presence-only metrics for BBS communities. 

#### 3.1.3 Phylogenetic

We calculated two distance-based metrics of phylogenetic diversity using functions in the R package picante: mean pairwise phylogenetic distance (MPD) and mean nearest-taxon phylogenetic distance (MNTD). We conducted 999 randomizations of the phylogenetic distance matrix and calculated the z-score of the observed phylogenetic distances relative to the distribution of phylogenetic distances of the randomized matrices. For the FIA dataset, we calculated both abundance-weighted and presence-only versions of MPD and MNTD, but for the BBS dataset we only calculated the presence-only version. 

### 3.2 Beta diversity 

We calculated taxonomic, functional, and phylogenetic beta-diversity (turnover of diversity among local communities) for tree and bird communities. Beta-diversity is defined as the variation in community composition across multiple local communities. To determine beta-diversity at a point, it is necessary to define the kernel or radius within which variation in community composition is taken into account. For the FIA dataset, we calculated beta-diversity among the four subplots within a plot. We also aggregated species abundances of each plot and calculated beta-diversity for each plot at a number of different radii around the focal plot; as the radius increases, the number of pairwise comparisons among plots also increases as more plots fall within the kernel. We calculated beta-diversity at the following radii for the FIA plots: 5, 10, 20, 50, 100 km. For the BBS plots, **I haven't yet determined how to implement beta-diversity for BBS, will add to methods when I do.**  

#### 3.2.1 Taxonomic

We calculated taxonomic beta-diversity using the d() function from the R package vegetarian (diversity partitioning method). We also calculated beta-diversity with the pairwise dissimilarity method using the vegdist() function from the R package vegan and taking the mean of the pairwise Jaccard dissimilarity of all local communities within a particular radius of the focal plot. As before, we calculated both abundance-weighted and presence-only metrics for FIA communities, and only presence-only metrics for BBS communities.

#### 3.2.2 Functional

We calculated functional beta-diversity using the comdist() function from the R package picante. We used the trait databases described above to generate pairwise distance matrices for all species in each survey. We used Gower distance to deal with missing trait values, also as described above. We calculated both pairwise and nearest-taxon functional beta-diversity, both abundance-weighted and presence-only for FIA, and just presence-only for BBS. Refer to Swenson et al. 2011, PLoS One, for definitions of the metrics we calculated.

#### 3.2.3 Phylogenetic

We calculated phylogenetic beta-diversity with the same set of metrics as we used for functional beta-diversity. We used the phylogenies described above to generate pairwise distance matrices for all species in each survey. 

