# Methods: calculating taxonomic diversity for FIA

Appendix for the conceptual paper. This does not have any information on the BBS dataset, on functional and phylogenetic diversity, and on the Baselga/Podani diversity partitioning methods, because none of the above are part of the conceptual paper. See the more complete methods document for those details.

Author: QDR  
Project: NASA Biodiversity  
Date created: 10 July 2017

## 1. Preparation of Forest Inventory and Analysis dataset

We obtained USFS Forest Inventory and Analysis (FIA) data for the Pacific Northwest region. Forest plots surveyed according to the FIA protocol consist of four subplots, each circles with 7.3m radius, located 36.6 m from one another in a three-pointed star pattern. See `https://www.fia.fs.fed.us/library/field-guides-methods-proc/docs/2016/core_ver7-1_10_2016-opt.pdf` for more detailed description of the survey protocols. Each tree in each subplot is identified to species, and its diameter at breast height is recorded. Using the diameters to calculate basal area of each individual tree, we summed the basal areas within each species to estimate the relative abundance of each species in each subplot. Any discrepancies in species names were resolved so that the most recent taxonomy is used.

## 2. Diversity calculations

We calculated diversity metrics based on species presence at each plot. We also calculated abundance (basal area)-weighted diversity metrics. We used the most recent survey as a single time point for each plot. 

We calculated alpha, beta, and gamma diversity at a number of different radii around each FIA plot by taking the median diversity of all plots in the radius, including the focal plot (alpha), the median pairwise diversity of all pairs of plots in the radius, including the focal plot (beta), and the aggregated diversity of all plots in the radius as if they were a single community (gamma). The radii for which we calculated diversities were 1, 5, 7.5, 10, 20, 50, 75, 100, 150, 200, and 300 km; a subset of these results are presented in the current manuscript. 

### 2.1 Alpha and gamma diversity

We calculated taxonomic alpha-diversity (diversity of a local community) for FIA tree communities. We calculated diversity indices for communities aggregated at the plot level (aggregating the four subplots making up one plot). For each plot and radius, we calculated alpha diversity within that radius by taking the median diversity value for all plots or routes (including the focal plot) located inside the circle defined by the radius around the focal plot. For gamma diversity, the diversity of a region that consists of multiple local communities, we aggregated all the plots within the focal circle to a single community, and calculated taxonomic, functional, and phylogenetic diversity of that community. We calculated non-abundance-weighted alpha and gamma diversity (species richness) by totaling the number of species in each community. We also calculated basal-area-weighted Shannon alpha and gamma diversity as follows: `$H' = -\sum_{i=1}^{R} p_i \ln p_i$`, or the natural logarithm of the true diversity with q = 1, where R is species richness and p<sub>i</sub> is the basal area of species i.

### 2.2 Beta diversity 

We calculated taxonomic beta-diversity (turnover of diversity among local communities) for FIA tree communities. Beta-diversity is defined as the variation in community composition across multiple local communities. To determine beta-diversity at a point, it is necessary to define the kernel or radius within which variation in community composition is taken into account. For the FIA dataset, we calculated beta-diversity among the four subplots within a plot. We also aggregated species abundances of each plot and calculated beta-diversity for each plot at a number of different radii around the focal plot; as the radius increases, the number of pairwise comparisons among plots also increases as more plots fall within the kernel. 

We calculated taxonomic beta-diversity using the beta.div.comp() function from the R package adespatial. We opted to use the pairwise distance-based method rather than Whittaker's diversity-partitioning method because it allows us to use more consistent methods across taxonomic, functional, and phylogenetic diversity calculations. We calculated multi-site indices based on both the Sorensen and Jaccard pairwise distances. Each of these represent the mean pairwise dissimilarity among plots or routes within a given radius. We were additionally interested in the relative contributions of nestedness and turnover to beta-diversity at different spatial scales. There are two commonly used methods to partition the multi-site index (Baselga & Orme 2010, GEB; Podani et al. 2013, Ecological Complexity). We used both methods to partition beta-diversity into these two components, again as implemented with the beta.div.comp() function.

*Potential future modification to methods*: We might eventually use rarefaction-based estimate of effective number of communities as the denominator in the multi-site mean pairwise dissimilarity metric, rather than the raw number of sites sampled. Also, we might eventually subsample the plots in larger radii to control for the fact that the number of pairwise comparisons increases with increasing radius.


