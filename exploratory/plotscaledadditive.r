# Plots for scaled additive, unscaled additive, and dissim beta coefficients for BBS

library(tidyverse)

load('~/Dropbox/projects/nasabiodiv/modelfits/testscalingadditivebetafit.RData')

scaled_fixef <- fit$coef %>% 
  filter(effect %in% 'fixed') %>%
  separate(parameter, into = c('response', 'parameter'), sep= '_', extra = 'merge') %>%
  select(response, parameter, stat, value) %>%
  group_by(response, parameter) %>%
  spread(stat, value)

model_fixef_all <- rbind(data.frame(model = 'scaled additive', scaled_fixef), data.frame(model = 'original', model_fixef_orig), data.frame(model = 'unscaled additive', model_fixef)) %>%
  mutate(schnignificant = (Q2.5>0 & Q97.5>0) | (Q2.5<0 & Q97.5<0),
         response = case_when(
           grepl('func',response) ~ 'functional',
           grepl('phy',response) ~ 'phylogenetic',
           TRUE ~ 'taxonomic'
         ))

param_labels <- c('climate: temp mean', 'climate: precip mean', 'geodiv: GPP', 'geodiv: elevation', 'geodiv: geological age', 'geodiv: soil type')

ggplot(model_fixef_all %>% filter(!parameter %in% 'Intercept', response %in% 'taxonomic', !model %in% 'unscaled additive'), aes(x = parameter, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange(aes(color = model), position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'blue') +
  scale_x_discrete(labels = param_labels) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = c('slateblue', 'indianred'), name = 'beta-diversity type', labels = c('scaled additive', 'dissimilarity')) +
  theme(legend.position = 'bottom')
