# Explore why precipitation has opposite effects on tree PD versus FD and TD by looking at species compositions
# QDR/NASAbioXgeo/21 aug 2018

# precipitation has a negative coefficient on PD but a positive one on TD and FD. Why?

# Load FIA and find the plots with high PD and low precip, as well as low PD and high precip

# Load FIA
# --------

load('/mnt/research/nasabio/temp/fia_spatial_mm_dat_50k.RData')

library(dplyr)
fia <- left_join(fiabio, fiageo)

# Scatterplot of phylogenetic diversity vs precipitation
# ------------------------------------------------------

with(fia, plot(bio12_5k_50_mean, alpha_phy_pa))
# Not really linear but you can see a big cluster in the upper left.

# Find plots that are in the top 10% of phylogenetic diversity and bottom 10% of precipitation
# --------------------------------------------------------------------------------------------

highdiv_lowprecip <- fia %>%
  filter(alpha_phy_pa >= quantile(alpha_phy_pa, 0.9, na.rm = TRUE) & bio12_5k_50_mean <= quantile(bio12_5k_50_mean, 0.1, na.rm = TRUE))
# Identifies 84 high phy diversity and low precipitation plots

lowdiv_highprecip <- fia %>%
  filter(alpha_phy_pa <= quantile(alpha_phy_pa, 0.1, na.rm = TRUE) & bio12_5k_50_mean >= quantile(bio12_5k_50_mean, 0.9, na.rm = TRUE))
# Identifies 13 low phy diversity and high precipitation plots

# Load raw FIA data so that we can get the species lists for each group of plots
# ------------------------------------------------------------------------------

fiaraw <- read.csv('/mnt/research/nasabio/data/fia/treedata10nov/finley_trees_continental_US_most_recent_evaluations_nov8_2017.csv', stringsAsFactors = FALSE)
sptable <- read.csv('/mnt/research/nasabio/data/fia/treedata10nov/lookup_table_allfia.csv', stringsAsFactors = FALSE) %>%
  select(FIA.Code, Common.name, binomial)

fiaraw_highdiv <- fiaraw %>%
  filter(PLT_CN %in% highdiv_lowprecip$PLT_CN)

fiaraw_lowdiv <- fiaraw %>%
  filter(PLT_CN %in% lowdiv_highprecip$PLT_CN)

highdiv_spp <- fiaraw_highdiv %>%
  select(SPCD) %>%
  rename(FIA.Code = SPCD) %>%
  left_join(sptable)

lowdiv_spp <- fiaraw_lowdiv %>%
  select(SPCD) %>%
  rename(FIA.Code = SPCD) %>%
  left_join(sptable)

sort(table(highdiv_spp$binomial))
sort(table(lowdiv_spp$binomial))