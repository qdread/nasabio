# NASABioXGeo metadata

Column descriptions
  
Created by: QDR  
Last modified: 07 Feb 2018

## Biodiversity

The same metrics were calculated for both BBS and FIA, with the exception that we calculated both incidence-based and abundance-based metrics for FIA, but only incidence-based metrics for BBS. FIA plots' composition is the composition at the most recent survey, while BBS routes' incidence for each bird species is based on pooled incidence from the past 10 surveys (2007-2016).

Most of the column names are in multiple CSV files so I list the column names for each CSV, followed by descriptions of the columns all together.

### List of CSV files

#### BBS

- `bbs_alphadiv.csv`

	Raw alpha-diversity values for each route (for our analyses, we aren't really using these, instead using the radius-averaged values).  
	**Column names**: "rteNo","lon","lat","lon_aea","lat_aea","richness","shannon","evenness","MPD_pa_z","MNTD_pa_z","MPD_z","MNTD_z","MPDfunc_pa_z","MNTDfunc_pa_z","MPDfunc_z","MNTDfunc_z"

- `bbs_alpha_1year.csv`

	Average alpha-diversity values for each route-radius combination (includes all routes within the radius).  
	**Column names**: "rteNo","lon","lat","lon_aea","lat_aea","radius","richness","MPD_pa_z","MNTD_pa_z","MPDfunc_pa_z","MNTDfunc_pa_z"

- `bbs_betatdpdfd_1year.csv`

	Average beta-diversity values for each route-radius combination (includes pairwise of all routes within the radius).  
	**Column names**: "rteNo","lon","lat","lon_aea","lat_aea","radius","beta_td_pairwise_pa","beta_td_sorensen_pa","beta_td_pairwise","beta_td_sorensen","beta_td_shannon","beta_pd_pairwise_pa","beta_pd_pairwise_pa_z","beta_pd_nt_pa","beta_pd_nt_pa_z","beta_pd_pairwise","beta_pd_pairwise_z","beta_pd_nt","beta_pd_nt_z","beta_fd_pairwise_pa","beta_fd_pairwise_pa_z","beta_fd_nt_pa","beta_fd_nt_pa_z","beta_fd_pairwise","beta_fd_pairwise_z","beta_fd_nt","beta_fd_nt_z"


- `bbs_gamma_1year.csv`

	Gamma-diversity values for each route-radius combination (includes all routes within the radius).  
	**Column names**: "rteNo","lon","lat","lon_aea","lat_aea","radius","richness","shannon","evenness","MPD_pa_z","MNTD_pa_z","MPD_z","MNTD_z","MPDfunc_pa_z","MNTDfunc_pa_z","MPDfunc_z","MNTDfunc_z"

#### FIA

- `fiausa_alphadiv.csv`

	Raw alpha-diversity values for each plot (for our analyses, we aren't really using these, instead using the radius-averaged values).  
	**Column names**: "PLT_CN","richness","shannon","evenness","MPD_pa_z","MNTD_pa_z","MPD_z","MNTD_z","MPDfunc_pa_z","MNTDfunc_pa_z","MPDfunc_z","MNTDfunc_z"

- `fiausa_alpha.csv`

	Average alpha-diversity values for each plot-radius combination (includes all plots within the radius).  
	**Column names**: "PLT_CN","radius","richness","shannon","evenness","MPD_z","MNTD_z","MPDfunc_z","MNTDfunc_z","MPD_pa_z","MNTD_pa_z","MPDfunc_pa_z","MNTDfunc_pa_z"

- `fiausa_betatd.csv`

	Average beta-diversity values for each plot-radius combination (includes pairwise of all plots within the radius). *This is taxonomic only, because functional and phylogenetic are still being calculated!*  
	**Column names**: "PLT_CN","radius","beta_td_pairwise_pa","beta_td_sorensen_pa","beta_td_pairwise","beta_td_sorensen"

- `fiausa_gamma.csv`

	Gamma-diversity values for each plot-radius combination (includes all plots within the radius).  
	**Column names**: "PLT_CN","radius","richness","shannon","evenness","MPD_pa_z","MNTD_pa_z","MPD_z","MNTD_z","MPDfunc_pa_z","MNTDfunc_pa_z","MPDfunc_z","MNTDfunc_z"

### Column descriptions

- `radius`: radius, in km, around the focal point where the data used to find the metrics or summary statistics were taken from

#### BBS plot info

- `rteNo`: Numerical ID of BBS route
- `lon`: longitude of BBS route centroid in decimal degrees
- `lat`: latitude of BBS route centroid in decimal degrees
- `lon_aea`: longitude of BBS route centroid in Albers equal-area projection for continental USA
- `lat_aea`: latitude of BBS route centroid in Albers equal-area projection for continental USA

#### FIA plot info

- `PLT_CN`: Numerical ID of FIA plot
- *No location information is in the FIA files because locations are confidential. It is stored elsewhere.*

#### Diversity 

Any metric with the abbrevation `pa` means presence-absence (a.k.a. incidence-based). Anything abundance-based will be `NA` for BBS. All z-scores are evaluated against a null distribution of diversity metrics found for that community.

##### Alpha and gamma

- `richness`: Total number of species present (taxonomic, incidence-based)
- `shannon`: Shannon diversity (taxonomic, abundance-based)
- `evenness`: Shannon evenness (taxonomic, abundance-based) (min 0, max 1)
- `MPD_pa_z`: Mean pairwise phylogenetic distance z-score (phylogenetic, incidence-based)
- `MNTD_pa_z`: Mean nearest-neighbor phylogenetic distance z-score (phylogenetic, incidence-based)
- `MPD_z`: Mean pairwise phylogenetic distance z-score (phylogenetic, abundance-based)
- `MNTD_z`: Mean nearest-neighbor phylogenetic distance z-score (phylogenetic, abundance-based)
- `MPDfunc_pa_z`: Mean pairwise functional distance z-score (phylogenetic, incidence-based)
- `MNTDfunc_pa_z`: Mean nearest-neighbor functional distance z-score (phylogenetic, incidence-based)
- `MPDfunc_z`: Mean pairwise functional distance z-score (phylogenetic, abundance-based)
- `MNTDfunc_z`: Mean nearest-neighbor functional distance z-score (phylogenetic, abundance-based)

##### Beta

- `beta_td_pairwise_pa`: Taxonomic, incidence-based, Jaccard index (min 0, max 1)
- `beta_td_sorensen_pa`: Taxonomic, incidence-based, Sorensen index (min 0, max 1)
- `beta_td_pairwise`: Taxonomic, abundance-based, Jaccard index (min 0, max 1)
- `beta_td_sorensen`: Taxonomic, abundance-based, Sorensen index (min 0, max 1)
- `beta_td_shannon`: Taxonomic, abundance-based, using method from *vegetarian* package
- `beta_pd_pairwise_pa`: Raw value of mean pairwise phylogenetic beta-diversity, incidence-based
- `beta_pd_pairwise_pa_z`: Z-score of mean pairwise phylogenetic beta-diversity, incidence-based
- `beta_pd_nt_pa`: Raw value of mean nearest-neighbor phylogenetic beta-diversity, incidence-based
- `beta_pd_nt_pa_z`: Z-score of mean nearest-neighbor phylogenetic beta-diversity, incidence-based
- `beta_pd_pairwise`: Raw value of mean pairwise phylogenetic beta-diversity, abundance-based
- `beta_pd_pairwise_z`: Z-score of mean pairwise phylogenetic beta-diversity, abundance-based 
- `beta_pd_nt`: Raw value of mean nearest-neighbor phylogenetic beta-diversity, abundance-based
- `beta_pd_nt_z`: Z-score of mean nearest-neighbor phylogenetic beta-diversity, abundance-based
- `beta_fd_pairwise_pa`: Raw value of mean pairwise functional beta-diversity, incidence-based
- `beta_fd_pairwise_pa_z`: Z-score of mean pairwise functional beta-diversity, incidence-based
- `beta_fd_nt_pa`: Raw value of mean nearest-neighbor functional beta-diversity, incidence-based
- `beta_fd_nt_pa_z`: Z-score of mean nearest-neighbor functional beta-diversity, incidence-based
- `beta_fd_pairwise`: Raw value of mean pairwise functional beta-diversity, abundance-based
- `beta_fd_pairwise_z`: Z-score of mean pairwise functional beta-diversity, abundance-based 
- `beta_fd_nt`: Raw value of mean nearest-neighbor functional beta-diversity, abundance-based
- `beta_fd_nt_z`: Z-score of mean nearest-neighbor functional beta-diversity, abundance-based

## Geodiversity

The same summary statistics were calculated for FIA plot centers and BBS route centroids, for each radius. Different geo variables all come in different native resolutions. For the ones that had a finer native resolution than 5 km pixel size, the variable is given in its native resolution and also aggregated into 5 km pixels. For each continuous variable (at both resolutions) and each radius, the following summary statistics are provided: 

- mean
- standard deviation
- minimum
- maximum
- mean terrain ruggedness index (TRI)
- mean topographic roughness

The TRI and roughness are calculated for all variables, not just elevation.

For the categorical variables (geological age and soil type), the following summary statistics are provided:

- mode (most common type within the radius)
- number of unique types (analogous to richness)
- Shannon diversity of proportions of types

### Types of variables

The variables fall in these categories:

- **Climate**: The 19 bioclim variables calculated directly from MODIS temperature data, as well as 8 "biocloud" variables we developed to be roughly analogous to some of the bioclim variables and calculated from cloudiness data.
- **Topography**: Elevation variables and things calculated from them.
- **Vegetation**: Variables associated with the Dynamic Habitat Index.
- **Geology/Soils**: Geological age and soil type
- **Human impacts**: Two different variables describing level of human influence (human footprint index and brightness of artificial lights)

### List of CSV files

### Column descriptions
