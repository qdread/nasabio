# Load plot data that was summarized on hpcc
fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/fia_unfuzzed'
load(file.path(fp, 'fiafitplotdat.R'))
fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_exploratory_plots'


library(cowplot)

# Plot the r2 values.
r2_lm_quant$radius_plot <- r2_lm_quant$radius + rep(c(-1,0,1), each = 5)

ggplot(r2_lm_quant, aes(x = radius_plot, group = interaction(radius, diversity_type), color = diversity_type)) +
 # geom_segment(aes(xend = radius_plot, y = r2_q025, yend = r2_q975)) +
  geom_segment(aes(xend = radius_plot, y = r2_q25, yend = r2_q75), size = 1) +
  geom_point(aes(y = r2), size = 2) +
  scale_color_discrete(name = 'Diversity', labels = c('alpha','beta','gamma')) +
  labs(x = 'Radius (km)', y = expression(R^2)) +
  theme(legend.position = c(0.1, 0.9))
  
# Get rid of values out of range.
xranges <- biogeo %>% group_by(radius) %>% summarize_at(vars(elevation_sd), funs(sdmin = min, sdmax = max), na.rm = TRUE)
pred_val_lm_quant <- pred_val_lm_quant %>%
  ungroup %>%
  left_join(xranges) %>%
  filter(x <= sdmax)

# Median (pseudo) r2s for plotting
r2_lm_quant <- r2_lm_quant %>%
  ungroup %>%
  group_by(diversity_type, radius) %>%
  mutate(r2expr = as.character(eval(substitute(expression(R^2 == x), list(x = round(r2, 2))))))

hexfill <- scale_fill_gradient(low = 'gray90', high = 'black')
hextheme <- theme(strip.background = element_blank())
border <- panel_border(colour = 'black')
radlabel <- labeller(radius = function(x) paste(as.integer(x), 'km'))

alphaplot <- ggplot(biogeo %>% mutate(alpha_diversity=exp(alpha_diversity))) + 
  facet_grid(~ radius, scales = 'free', labeller = radlabel) +
  geom_hex(aes(x = elevation_sd, y = alpha_diversity)) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'alpha_diversity'), 
              aes(x = x, ymin = pred_y_q025, ymax = pred_y_q975), fill = 'red', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'alpha_diversity'), 
            aes(x = x, y = pred_y), color = 'red', size = 1.5) +
  geom_text(data = subset(r2_lm_quant, diversity_type == 'alpha_diversity'),
            aes(x = -Inf, y = Inf, label = r2expr),
            parse = TRUE, hjust = -0.5, vjust = 2) +
  hexfill + hextheme + border + labs(x = 'Elevation standard deviation', y = 'Alpha-diversity')

betaplot <- ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free', labeller = radlabel) +
  geom_hex(aes(x = elevation_sd, y = beta_diversity)) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'beta_diversity'), 
              aes(x = x, ymin = pred_y_q025, ymax = pred_y_q975), fill = 'red', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'beta_diversity'), 
            aes(x = x, y = pred_y), color = 'red', size = 1.5) +
  geom_text(data = subset(r2_lm_quant, diversity_type == 'beta_diversity'),
            aes(x = Inf, y = -Inf, label = r2expr),
            parse = TRUE, hjust = 1.1, vjust = -2) +
  hexfill + hextheme + border + labs(x = 'Elevation standard deviation', y = 'Beta-diversity')

gammaplot <- ggplot(biogeo %>% mutate(gamma_diversity = exp(gamma_diversity))) + 
  facet_grid(~ radius, scales = 'free', labeller = radlabel) +
  geom_hex(aes(x = elevation_sd, y = gamma_diversity)) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'gamma_diversity'), 
              aes(x = x, ymin = pred_y_q025, ymax = pred_y_q975), fill = 'red', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'gamma_diversity'), 
            aes(x = x, y = pred_y), color = 'red', size = 1.5) +
  geom_text(data = subset(r2_lm_quant, diversity_type == 'gamma_diversity'),
            aes(x = -Inf, y = Inf, label = r2expr),
            parse = TRUE, hjust = -0.5, vjust = 2) +
  hexfill + hextheme + border + labs(x = 'Elevation standard deviation', y = 'Gamma-diversity')

ggsave(file.path(fpfig, 'fia_alpha_regressions.png'), alphaplot, height = 4, width = 12, dpi = 300)
ggsave(file.path(fpfig, 'fia_beta_regressions.png'), betaplot, height = 4, width = 12, dpi = 300)
ggsave(file.path(fpfig, 'fia_gamma_regressions.png'), gammaplot, height = 4, width = 12, dpi = 300)
