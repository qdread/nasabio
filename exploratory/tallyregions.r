# Get rid of regions with less than 5 points in them.

# Count up how many sites are in each region

load('C:/Users/Q/Dropbox/projects/nasabiodiv/bbs_spatial_mm_dat.RData')

bbs_n_bcr <- bbsgeo %>% group_by(BCR) %>% summarize(n = n())
bbs_n_tnc <- bbsgeo %>% group_by(TNC) %>% summarize(n = n()) # Three have only one 
bbs_n_huc <- bbsgeo %>% group_by(HUC4) %>% summarize(n = n()) # 17 have less than 5. We could get rid of them?

load('C:/Users/Q/Dropbox/projects/nasabiodiv/fia_spatial_mm_dat.RData')
fia_n_bcr <- fiageo %>% group_by(BCR) %>% summarize(n = n())
fia_n_tnc <- fiageo %>% group_by(TNC) %>% summarize(n = n()) 
fia_n_huc <- fiageo %>% group_by(HUC4) %>% summarize(n = n()) 

countitself <- function(x) sapply(x, function(z) sum(x %in% z))

bbstallies <- bbsgeo %>%
  select(rteNo, HUC4, BCR, TNC) %>%
  mutate(n_BCR = countitself(BCR),
         n_TNC = countitself(TNC),
         n_HUC4 = countitself(HUC4),
         remove = n_BCR < 5 | n_TNC < 5 | n_HUC4 < 5)

fiatallies <- fiageo %>%
  select(PLT_CN, HUC4, BCR, TNC) %>%
  mutate(n_BCR = countitself(BCR),
         n_TNC = countitself(TNC),
         n_HUC4 = countitself(HUC4),
         remove = n_BCR < 5 | n_TNC < 5 | n_HUC4 < 5)

# Redo tallies
bbstallies %>%
  filter(!remove) %>%
  mutate(nnew_BCR = countitself(BCR), nnew_TNC = countitself(TNC), nnew_)

# Winnow down datasets, removing the offending entries

