# Exploratory data analysis of BBS and MODIS data sets
# Start with by-route data, then switch to by-stop data (larger)

# Modified 15 Feb 2017: Generate a subset of plots with the key relationships on them, with some nicer aesthetics, for the powerpoint.
# Note: At this point, the data are still hosted on the aquaxterra space but can later be moved to the nasabio space.

fp <- '/mnt/research/aquaxterra/DATA/raw_data/BBS'
bm <- read.csv(file.path(fp, 'bbs_allcovariates_byroute.csv'), stringsAsFactors = FALSE)
bd <- read.csv(file.path(fp, 'bbs_div_byroute.csv'), stringsAsFactors = FALSE)

# Join BBS diversity and MODIS/PRISM data by route number and year.

# Convert MODIS df to wide format, similar to the BBS df.

library(dplyr)
library(reshape2)

bmwide <- dcast(bm, rteNo + year ~ variable)

# This joined data frame contains all needed to run exploratory correlational analysis (incl. raw, spatial, and time-series)
bmd <- full_join(bmwide, bd)
bmd2 <- filter(bmd, year >= 2001, year <= 2011)

# For raw analysis, average out by year
bmdyear <- bmd2[,1:51] %>% group_by(rteNo) %>% summarize_all(.funs = mean, na.rm = TRUE) 

# Interannual CV of some of the variables.
bmdcvs <-
bmd2[,1:51] %>% group_by(rteNo) %>% summarize(bio01cv = sd(bio01+273, na.rm=T)/mean(bio01+273, na.rm=T),
											  bio12cv = sd(bio12, na.rm=T)/mean(bio12, na.rm=T),
											  laicv = sd(lai_max_aeaproj, na.rm=T)/mean(lai_max_aeaproj, na.rm=T),
											  nppcv = sd(log10(MOD17A3_Science_NPP), na.rm=T)/mean(log10(MOD17A3_Science_NPP), na.rm=T))

bmdyear <- left_join(bmdyear, bmdcvs)											  
										
# Log transform NPP

bmdyear <- mutate(bmdyear, npplog = log10(MOD17A3_Science_NPP))
										
# Group variable column indices by type (split different diversity types into groups as well).
bioclim_idx <- 3:21
modisnlcd_idx <- 22:29
tdiv_idx <- 30:33
fdiv_idx <- 34:37
pdiv_idx <- c(39,40,44,46,50)

pairwise_bygroup <- function(dat, x, y) {
	pairwise_cors <- matrix(nrow = length(x), ncol = length(y))
	dimnames(pairwise_cors) <- list(names(dat)[x], names(dat)[y])

	pb <- txtProgressBar(0, length(x)*length(y), style=3)

	for (i in 1:length(x)) {
		for (j in 1:length(y)) {
			pairwise_cors[i, j] <- cor(dat[, x[i]], dat[, y[j]], use = 'complete.obs')
			setTxtProgressBar(pb, j + length(y)*(i-1))
		}
	}

	close(pb)
	return(pairwise_cors)
} 

# Go fishing!
bc_td <- pairwise_bygroup(bmdyear, bioclim_idx, tdiv_idx)
bc_fd <- pairwise_bygroup(bmdyear, bioclim_idx, fdiv_idx)
bc_pd <- pairwise_bygroup(bmdyear, bioclim_idx, pdiv_idx)
mn_td <- pairwise_bygroup(bmdyear, modisnlcd_idx, tdiv_idx)
mn_fd <- pairwise_bygroup(bmdyear, modisnlcd_idx, fdiv_idx)
mn_pd <- pairwise_bygroup(bmdyear, modisnlcd_idx, pdiv_idx)

# Modified pairs plots
pairplot_bygroup <- function(dat, x, y, xlabels, ylabels, filename, directory) {
	require(ggplot2)
	require(gridExtra)
	
	plot_list <- list()
	
	for (i in 1:length(y)) {
		for (j in 1:length(x)) {
			th <- theme_bw() + theme(panel.grid = element_blank())
			if (i != length(y)) th <- th + theme(axis.title.x = element_blank())
			if (j != 1) th <- th + theme(axis.title.y = element_blank())
			plot_list[[length(plot_list) + 1]] <- ggplot(dat, aes_string(x = names(dat)[x[j]], y = names(dat)[y[i]])) + geom_point(alpha = 0.25) + stat_smooth() + labs(x=xlabels[j], y=ylabels[i]) + th
		}
	}
	
	png(file.path(directory, filename), height = length(y)*2.5, width = length(x)*2.5, units = 'in', res = 200)
	grid.arrange(grobs = plot_list, nrow = length(y))
	dev.off()

} 



# pairplot_bygroup(bmdyear, bioclim_idx, tdiv_idx, 'bioclim_taxdiv_pairs.png')
# pairplot_bygroup(bmdyear, bioclim_idx, fdiv_idx, 'bioclim_fundiv_pairs.png')
# pairplot_bygroup(bmdyear, bioclim_idx, pdiv_idx, 'bioclim_phydiv_pairs.png')
# pairplot_bygroup(bmdyear, modisnlcd_idx, tdiv_idx, 'modisnlcd_taxdiv_pairs.png')
# pairplot_bygroup(bmdyear, modisnlcd_idx, fdiv_idx, 'modisnlcd_fundiv_pairs.png')
# pairplot_bygroup(bmdyear, modisnlcd_idx, pdiv_idx, 'modisnlcd_phydiv_pairs.png')

# 15 Feb: Pruned down indices.

# Predictors (x)
temp_idx <- c(3, 13, 6, 52)
temp_lbl <- c('Temperature mean', 'Winter temperature', 'Temperature seasonality', 'Temperature yearly CV')
precip_idx <- c(14, 16, 17, 53)
precip_lbl <- c('Rainfall mean', 'Driest month rainfall', 'Rainfall seasonality', 'Rainfall yearly CV')
prod_idx <- c(24, 54, 56, 55)
prod_lbl <- c('LAI mean', 'LAI yearly CV', 'NPP mean (log)', 'NPP yearly CV')

# Responses (y)
tdiv_idx <- c(30, 31)
tdiv_lbl <- c('Species richness', 'Shannon diversity')
fdiv_idx <- c(34, 37)
fdiv_lbl <- c('Functional richness', 'Functional dispersion')
pdiv_idx <- c(39, 44)
pdiv_lbl <- c('Total PD', 'Pairwise PD z-score')

fp_fig <- '/mnt/research/nasabio/figs/pairplots'

pairplot_bygroup(bmdyear, temp_idx, tdiv_idx, temp_lbl, tdiv_lbl, 'temperature_taxdiv.png', fp_fig)
pairplot_bygroup(bmdyear, precip_idx, tdiv_idx, precip_lbl, tdiv_lbl, 'precip_taxdiv.png', fp_fig)
pairplot_bygroup(bmdyear, prod_idx, tdiv_idx, prod_lbl, tdiv_lbl, 'productivity_taxdiv.png', fp_fig)
pairplot_bygroup(bmdyear, temp_idx, fdiv_idx, temp_lbl, fdiv_lbl, 'temperature_funcdiv.png', fp_fig)
pairplot_bygroup(bmdyear, precip_idx, fdiv_idx, precip_lbl, fdiv_lbl, 'precip_funcdiv.png', fp_fig)
pairplot_bygroup(bmdyear, prod_idx, fdiv_idx, prod_lbl, fdiv_lbl, 'productivity_funcdiv.png', fp_fig)
pairplot_bygroup(bmdyear, temp_idx, pdiv_idx, temp_lbl, pdiv_lbl, 'temperature_phydiv.png', fp_fig)
pairplot_bygroup(bmdyear, precip_idx, pdiv_idx, precip_lbl, pdiv_lbl, 'precip_phydiv.png', fp_fig)
pairplot_bygroup(bmdyear, prod_idx, pdiv_idx, prod_lbl, pdiv_lbl, 'productivity_phydiv.png', fp_fig)

# Plot of relationship among diversity variables

library(GGally)
png(file.path(fp_fig, 'bbsdiv_ggpairs.png'), height=10, width=10, res=400, units='in')
ggpairs(bmdyear[, c('richness','shannon','FRic','FDis','PD','mpd.obs.z')],
        diag=list(continuous=wrap('barDiag', bins=15))) + theme_bw()
dev.off()
