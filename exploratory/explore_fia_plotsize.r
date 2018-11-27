# Exploratory script to see if anything is different with the different FIA plots
# Use raw data

fia <- read.csv('/mnt/research/nasabio/data/fia/treedata10nov/finley_trees_continental_US_most_recent_evaluations_nov8_2017.csv', stringsAsFactors = FALSE)

library(dplyr)
library(psych)

fia_rich <- fia %>%
	group_by(TPA_UNADJ, PLT_CN) %>%
	summarize(richness = length(unique(SPCD)),
			  n = n())

qs <- function(x) data.frame(q025 = quantile(x, 0.025), q25 = quantile(x, 0.25), q50 = quantile(x, 0.5), q75 = quantile(x, 0.75), q975 = quantile(x, 0.975))			  
			  
fia_rich %>%
	filter(!is.na(TPA_UNADJ)) %>%
	group_by(TPA_UNADJ) %>%
	do(describe(.$richness))

fia_rich %>%
	filter(!is.na(TPA_UNADJ)) %>%
	group_by(TPA_UNADJ) %>%
	do(qs(.$richness))
	
fia %>%
	filter(!is.na(TPA_UNADJ)) %>%
	group_by(TPA_UNADJ) %>%
	summarize(n_plot = length(unique(PLT_CN)),
			  n_state = length(unique(STATECD)))	


library(ggplot2)

p <- ggplot(fia_rich, aes(x = richness, group = factor(TPA_UNADJ), fill = factor(TPA_UNADJ))) +
	geom_density(adjust = 2) + 
	theme_bw() +
	scale_x_log10() +
	scale_y_continuous(limits = c(0, 2), expand = c(0, 0))
	
ggsave('/mnt/research/nasabio/figs/richness_by_plotsize.png', p, height = 5, width = 7, dpi = 300)

# See whether they are contained in each other
tpa <- sort(unique(fia$TPA_UNADJ))

plot_big <- unique(fia$PLT_CN[fia$TPA_UNADJ == tpa[1]])
plot_med <- unique(fia$PLT_CN[fia$TPA_UNADJ == tpa[2]])
plot_small <- unique(fia$PLT_CN[fia$TPA_UNADJ == tpa[3]])

table(plot_big %in% plot_med)
table(plot_small %in% plot_med)

# Make a table by state of whether each plot has a macroplot
plot_names <- c('macro','normal','micro')
fia <- fia %>%
	mutate(plot_type = plot_names[match(TPA_UNADJ, tpa)])
	
micro_macro <- fia %>%
	group_by(STATECD, PLT_CN) %>%
	summarize(micro = 'micro' %in% plot_type,
			  macro = 'macro' %in% plot_type)
			  
micro_macro_bystate <- micro_macro %>%
	ungroup %>%
	group_by(STATECD) %>%
	summarize(n_with_micro = sum(micro),
			  n_with_macro = sum(macro),
			  n_with_neither = sum(!micro & !macro))
			  
# Ou sont les macroplots d'antan? 
# (map with fuzzed locations)
map_dat <- fia %>%
	group_by(STATECD, PLT_CN, FUZZ_LON, FUZZ_LAT) %>%
	summarize(micro = 'micro' %in% plot_type,
			  macro = 'macro' %in% plot_type)

the_map <- ggplot(map_dat, aes(x = FUZZ_LON, y = FUZZ_LAT, color = macro)) +
	borders('state') +
	geom_point() +
	theme_bw() +
	coord_map()
	