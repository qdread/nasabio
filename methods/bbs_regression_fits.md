# BBS diversity versus geodiversity regressions

QDR, 19 Dec. 2017

Here are some new plots using the pooled incidence from 2007-2016 to get route-level diversity values for BBS, regressed against geodiversity as measured by standard deviation of elevation within radii, ranging from 50 km-200 km. It turns out that elevation standard deviation actually is a very strong predictor of beta- and gamma-diversity for BBS. The relationship with gamma-diversity is more scale-dependent than for beta-diversity. These plots are made with the same method as the ones that are now in the conceptual paper for FIA (the main difference is that these are based on presence-absence, not abundance). First, I threw out any BBS routes that were within 200 km of Mexico or Canada because there would be some missing areas that should factor in to the diversity around those points. Then, I subsampled the remaining ~3000 routes to make sure that 200-km radii around the points did not overlap, to avoid spatial autocorrelation. As it turns out, you can fit 35 to 40 non-overlapping circles of that kind in the continental USA, so each of the 999 subsamples I took was used to fit a regression with approximately n = 35 or so. Just like with the FIA data, the alpha and gamma are simple linear regressions while the beta is a beta-regression with logit link to account for the fact that the response variable is bounded between 0 and 1. I was surprised at how well elevation variability predicts BBS diversity especially at the larger spatial scales. See the figures below.

![Taxonomic alpha-diversity](/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_diversity_maps/bbs_alpha_regressions.png)  
**Figure 1** BBS taxonomic alpha-diversity (richness) versus geodiversity (elevation SD) across the USA, by radius.

![Taxonomic beta-diversity](/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_diversity_maps/bbs_beta_regressions.png)  
**Figure 2** BBS taxonomic beta-diversity versus geodiversity (elevation SD) across the USA, by radius.

![Taxonomic gamma-diversity](/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_diversity_maps/bbs_gamma_regressions.png)  
**Figure 3** BBS taxonomic gamma-diversity (richness) versus geodiversity (elevation SD) across the USA, by radius.
 
![R-squared of fits](/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_diversity_maps/bbs_r2s.png)  
**Figure 4** R-squared values (or pseudo-r-squared values for beta regression) of the fits.