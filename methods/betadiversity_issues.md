# Beta-diversity issues for NASABioxgeo project

QDR, 28 June 2017

In our call on 23 June, we discussed some methodological/technical snafus that have come up with beta-diversity. Some of it gets into too much technical detail to be easily explained over the phone, so I decided to quickly write up the issues to see if anyone had feedback on them. Also, since the call, I heard back from Philip Barton who wrote a paper on scaling beta-diversity (in our google drive folder at `NASABiodiversityWG/Literature/Beta diversity`). I am incorporating his insightful input into this document.

Overall, Barton seemed to think that there aren't any major issues with using the data we have, though he probably isn't too familiar with the particular datasets. He said:

> Looking through your email I see no major issues for you to proceed with your analysis, largely because you are dealing with bird community data thats typically of quite high quality in terms of completeness. The issue, if I understand correctly, is that you are planning to use data from a breeding bird survey, and this may not have been done in a systematic and standardised way? Is this right? I would think that choosing the best subset of samples that are most similar in terms of survey effort yet spatially representative of your study area is the first challenge and probably the most important step.

> Without knowing how far you are into this beta diversity business, you might also like to look at the work of Anne Chao in Taiwan, Andres Baselga in Spain, Hanna Tuomisto in Finland, and Karel Mokany and Simon Ferrier here in Australia. They do some great work looking at community turnover and use a variety of techniques. My first thought when reading about your project was that Generalised Dissimilarity Modelling would be appropriate. How you compare the patterns and drivers across scales really comes down to how you aggregate samples, treat your data relative to grain and extent, and what spatial and environmental predictors you have available.

## Issue 1. Sample size increases as radius increases.

From our preliminary data, it looks like the absolute value of beta-diversity steadily increases with increasing radius. 

![beta div figure](file:///C:\\Users\\Q\\google_drive\\NASABiodiversityWG\\Figures\\betadiv\\data_paper_maps\\fig2_bottomrow_fits.png)

You can see in the figure, which shows FIA data, that the 100 km radius has no low values of pairwise beta-diversity. This may be a true pattern, but it also may have to do with the increased sample size as the radius increases. Many more plots are being compared in each data point in the 100 km figure. In later versions, we will use the multi-site Sorensen index which is suggested by Barton to deal with this. However, these plots are using a multi-site Jaccard index which really should not differ much. Another possible way to deal with the increased number of plots would be to take many subsamples of plots in the larger radii equal to the number of plots in the smaller radius, and calculate beta-diversity on those reduced sets so that sample size would be the same regardless of radius. If anyone has any ideas about that, that'd be great.

## Issue 2. Rarefaction.

The Barton paper states the following on rarefaction.

> We recognize that when scaling curves are constructed from empirical samples, as will be necessary in practice, then the number of sampling units will often incompletely represent the ‘true’ number of communities, and will require standardization by rarefaction or extrapolation (Colwell et al., 2012). This must be considered prior to the calculation of a normalized differentiation measure, such as one minus the multiple-site Sørensen index (Chao et al., 2012), and will improve comparability of beta-diversity values across different studies.

I asked Barton whether this meant to use rarefaction to estimate the effective number of communities, and use that as the denominator when calculating mean pairwise beta-diversity. He seemed to agree with this and also said:

> Yes, I would agree with your thinking here. However, because you will using bird data, this is perhaps not a big issue because bird surveys tend to generate very complete data sets relative to most other taxa. I would first run some accumulation curves and generate some richness estimates and just see how much discrepancy there is between estimated and observed richness in your samples. Also consider removing singleton species, and use only presence/absence data. For birds, there are other issues as well, such as whether water birds or raptors should be excluded from terrestrial areas or samples at fine spatial grains. I would emphasise that this kind of ecological robustness of your samples is just as important as the statistical robustness and issues of sampling completeness.

## Issue 3. Consistency of our methods across TD, PD, and FD.

Our goal in this project is to look at the taxonomic, functional, and phylogenetic aspects of beta-diversity. However, a lot of the technical methods papers on beta-diversity that we are relying on only talk about taxonomic diversity. We are unsure how to make sure our methods are consistent across the three "flavors" of diversity. We are using Barton's recommended methods for taxonomic diversity but still using methods from the vegan package to do functional and phylogenetic diversity.

Barton's opinion on the topic:

> I have never performed any kind of rarefaction prior to analysis of FD, and I am less familiar with phylogenetic beta diversity methods. These data all derive from your species survey data, so again, a conservative approach might be to first ensure you have removed singletons etc, and you are confident that your raw data is treated appropriately and of high quality before delving into rarefaction of PD or FD. This paper might provide further advice: https://link.springer.com/chapter/10.1007/978-3-319-22461-9_10

## Issue 4. Accounting for imperfect detection in BBS.

In the call, we also discussed the problem of imperfect detection in the BBS data, which is related to beta-diversity because it is possible that it may bias our results. There could be some bird species present at each route that weren't detected because they're lurking stealthily in the bushes. A way that it has been addressed is by doing occupancy modeling to determine detection probabilities for each species at each site. Instead of working with the raw observed presence-absence values from the survey, we would work instead with the modeled occupancy values. This is a hyper-conservative approach that would be very labor-intensive to implement. See the paper by Jarzyna and Jetz for how it was done in a recent paper, on our google drive in the folder `NASABiodiversityWG/Literature/BBS` (the appendices are also in that folder). Quentin also has the code that Jarzyna used to fit the models, if anyone is interested.

The consensus we had was that it's probably not necessary, though we could always do it if reviewers force us to. Looking at the methods in the Jarzyna and Jetz paper, it seems like they had to make quite a few assumptions to do the occupancy modeling. If you look at their appendices, you can see figures where they compared diversity and turnover estimates with and without the occupancy modeling. There are no statistics provided for whether the estimates are significantly different, but it does not seem like skipping the occupancy models would have changed the overall inference much. If anyone would like to volunteer to work on the occupancy models, that would be much appreciated.