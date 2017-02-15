# Calculate diversity from FIA PNW.

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'

fiapnw <- read.csv(file.path(fp, 'finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)

library(dplyr)
library(vegan)

# Get rid of dead trees. STATUSCD = 1 means alive, 2 means dead
# Calculate diversity by plot and subplot

# Calculate basal area at subplot level

fiasums_subplot <- fiapnw %>%
  filter(STATUSCD == 1) %>%
  mutate(ba = pi * (DIA/200)^2) %>% # basal area in m2.
  group_by(STATECD, COUNTYCD, PLOT, SUBP, SPCD) %>%
  summarize(basalarea = sum(ba),
            n = n())
  
# Calculate diversity at subplot level

subplot_diversity <- fiasums_subplot %>% ungroup %>%
  group_by(STATECD, COUNTYCD, PLOT, SUBP) %>%
  summarize(richness = length(unique(SPCD)),
            shannon_basalarea = diversity(basalarea, index = 'shannon'),
            evenness_basalarea = shannon_basalarea/log(richness),
            shannon_n = diversity(n, index = 'shannon'),
            evenness_n = shannon_n/log(richness))

# Calculate basal area at plot level

fiasums_plot <- fiapnw %>%
  filter(STATUSCD == 1) %>%
  mutate(ba = pi * (DIA/200)^2) %>% # basal area in m2.
  group_by(STATECD, COUNTYCD, PLOT, SPCD) %>%
  summarize(basalarea = sum(ba),
            n = n())

# Calculate diversity at plot level

plot_diversity <- fiasums_plot %>% ungroup %>%
  group_by(STATECD, COUNTYCD, PLOT) %>%
  summarize(richness = length(unique(SPCD)),
            shannon_basalarea = diversity(basalarea, index = 'shannon'),
            evenness_basalarea = shannon_basalarea/log(richness),
            shannon_n = diversity(n, index = 'shannon'),
            evenness_n = shannon_n/log(richness)) 

save(plot_diversity, subplot_diversity, file = file.path(fp, 'fia_diversitymetrics.RData'))
   