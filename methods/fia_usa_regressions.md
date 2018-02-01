# FIA diversity maps and biodiversity-geodiversity regressions: whole USA

QDR, 16 January 2018

All the stats and maps we generated so far for FIA data have been for the Pacific Northwest (PNW) only because (1) it was the first region where we could get access to the true plot locations and (2) it is a relatively manageable area to work with in terms of both low species richness and not having to deal with all 135,000 forested FIA plots across the continental United States. I mapped out alpha, beta, and gamma biodiversity (taxonomic only for starters) across the 22,000 forested FIA plots in the PNW, and the first guess at a relevant geodiversity variable which was standard deviation of elevation within all the different radii and also ran regressions using elevation geodiversity to predict the 3 flavors of biodiversity. You can see those results in the conceptual paper: in the PNW, there is a strong positive relationship with gamma diversity especially at larger scales (50-100 km), a moderately strong positive relationship with beta diversity at intermediate to large scales, and little to no relationship with alpha diversity at any scale.

After we got access to the true coordinates for the rest of the USA, I calculated all the biodiversity and elevation geodiversity metrics for the remaining plots. As you can see from the maps, the story is not as clean as when you are looking at just the PNW. For one thing, the USA spans a number of different biomes that have very different baseline tree species diversity. The gamma and beta diversity trends within a region are probably fairly consistent: for example the Appalachians have higher beta+gamma diversity than the surrounding lowlands, and the western mountain ranges also do. However if you look at the entire USA in a single regression those two within-biome patterns are thrown together and there is no overall trend. My next goal, once the rest of the covariates are calculated, is to add additional predictors that will split up the East and West coast so that we can tease out those patterns. But for now, here are some of the maps and regression plots for you to look at.

## Maps

I am only showing the 50 km maps because it illustrates the pattern best. I also made maps and analyses for different metrics and for incidence-based versus abundance-based metrics, but the results were pretty similar. You can see that the beta and gamma relationship is not really the same for the East Coast forests. There are some very flat and homogeneous areas with very high beta and gamma diversity, not to mention that the entire scale of it is different on the East Coast (East Coast forests with low tree diversity are more diverse than the high-diversity forests on the West Coast).

![Taxonomic alpha-diversity map](/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_diversity_maps/fia_usa_alpha_shannon_50_km.png)  
**Figure 1** Map of FIA taxonomic alpha-diversity (effective species number, q=1, or "true diversity" equivalent of Shannon diversity).

![Taxonomic beta-diversity map](/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_diversity_maps/fia_usa_beta_sorensen_50_km.png)  
**Figure 2** Map of FIA taxonomic beta-diversity (S&oslash;rensen, by abundance).

![Taxonomic gamma-diversity map](/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_diversity_maps/fia_usa_gamma_shannon_50_km.png)  
**Figure 3** Map of FIA taxonomic gamma-diversity (effective species number, q=1, or "true diversity" equivalent of Shannon diversity).

![Elevation geodiversity map](/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_diversity_maps/fia_usa_stdev_elev_50_km.png)  
**Figure 4** Map of FIA elevation geodiversity (SD of elevation).

## Regression plots

The regressions were done by doing many subsamples of plots with non-overlapping 100-km radii. This usually results in about 150 plots that cover the USA. You can see that there is little to no relationship for any of the diversity levels, because of the confused signal when putting all the USA with its varied biomes into a single analysis. 

![Taxonomic alpha-diversity regression](/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_exploratory_plots/USA/fia_alpha_regressions.png)  
**Figure 5** FIA taxonomic alpha-diversity (effective species number, q=1, or "true diversity" equivalent of Shannon diversity) versus geodiversity (elevation SD) across the USA, by radius.

![Taxonomic beta-diversity regression](/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_exploratory_plots/USA/fia_beta_regressions.png)  
**Figure 6** FIA taxonomic beta-diversity (S&oslash;rensen, by abundance) versus geodiversity (elevation SD) across the USA, by radius.

![Taxonomic gamma-diversity regression](/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_exploratory_plots/USA/fia_gamma_regressions.png)  
**Figure 7** FIA taxonomic gamma-diversity (effective species number, q=1, or "true diversity" equivalent of Shannon diversity) versus geodiversity (elevation SD) across the USA, by radius.
 
![R-squared of fits](/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_exploratory_plots/USA/fia_lm_rsquared.png)  
**Figure 8** R-squared values (or pseudo-r-squared values for beta regression) of the fits.