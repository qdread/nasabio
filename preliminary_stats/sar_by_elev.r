# Species area curve slope: fit SAR to FIA plots at different radii
# Then use the z coefficient from SAR as response variable in geodiversity regression

# Load gamma diversity and elevational diversity
fia_gamma <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/fia_gammadiv.csv', stringsAsFactors = FALSE)
ed <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/fia_elev_stats_noalaska.csv', stringsAsFactors = FALSE)


# Function to fit SAR (S = cA^z), return intercept, slope, and R^2
fit_sar <- function(A, S) {
  sar <- lm(I(log(S)) ~ I(log(A)))
  rsq <- summary(sar)$r.sq
  return(data.frame(c = sar$coef[1], z = sar$coef[2], r.squared = rsq))
}

library(dplyr)

fia_sars <- fia_gamma %>%
  mutate(area = pi * radius^2) %>%
  filter(richness > 0) %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT, lat, lon) %>%
  do(fit_sar(A = .$area, S = .$richness))

# Look at how the SAR slope at each point is related to elevation diversity at different radii

# The SAR slopes are the same every time, but I used different x values for the different radii just to look.

library(cowplot)
ed %>%
  filter(radius >= 5) %>%
  left_join(fia_sars) %>%
  ggplot(aes(x = sd_elev, y = z)) +
    geom_point() +
    stat_smooth(color = 'red', method = 'lm') +
    facet_wrap(~ radius, labeller = labeller(radius = function(x) paste(x, 'km'))) +
    panel_border(colour = 'black') +
    theme(strip.background = element_blank()) +
    scale_x_continuous(breaks = c(0,500,1000)) +  
    labs(x = 'Standard deviation of elevation', y = 'z (slope of SAR)')
  
# It looks like the SAR slope is lower where elevational variability is higher.
# This means species accumulate more slowly in more "geodiverse" areas, interesting.

ed %>%
  filter(radius >= 5) %>%
  left_join(fia_sars) %>%
  group_by(radius) %>%
  summarize(rsq = summary(lm(z ~ sd_elev))$r.sq)


# Added 22 Sep: relationship between SAR intercept and ED, and between SAR intercept and SAR slope.

ed %>%
  filter(radius >= 5) %>%
  left_join(fia_sars) %>%
  ggplot(aes(x = sd_elev, y = c)) +
  geom_point() +
  stat_smooth(color = 'red', method = 'lm') +
  facet_wrap(~ radius, labeller = labeller(radius = function(x) paste(x, 'km'))) +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  scale_x_continuous(breaks = c(0,500,1000)) +  
  labs(x = 'Standard deviation of elevation', y = 'c (log intercept of SAR)')

ggplot(fia_sars, aes(x = c, y = z)) +
  geom_point() +
  panel_border(colour = 'black') +
  labs(x = 'c (intercept)', y = 'z (slope)')

fia_gamma %>%
  filter(radius == 5) %>%
  left_join(fia_sars) %>%
  ggplot(aes(x = richness, y = c)) +
  geom_point()

fia_gamma %>%
  filter(radius == 5) %>%
  left_join(fia_sars) %>%
  ggplot(aes(x = richness, y = z)) +
  geom_point()
