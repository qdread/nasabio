# Calculate median beta-diversity of fia by radius.
# Modified 25 June to include PD and FD.
# Modified 09 Sep to load a different object
# Modified 04 Dec for the unfuzzed coordinates
# Modified 05 Dec to use the mean rather than the median
# Modified 06 Dec to use arcsin sqrt
# Modified 08 Jan 2018: for taxonomic only, and load only the slice that's needed at the time (parallel job).
# New version created 01 Mar 2018: keep only natural plots, don't use plantation plots
# Another version created 04 Apr 2018: do all metrics separately now that TD, PD, and FD are all done. 
# Yet another version created 28 Nov 2018: new version with new FIA data (macroplots removed) and slurm IDs
# Load FIA coordinates

metric <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Determine if arcsin sqrt transformation needs to be applied (if variable is on proportion scale)
prop_vars <- c('beta_td_pairwise', 'beta_td_sorensen', 'beta_td_pairwise_pa', 'beta_td_sorensen_pa', 'beta_fd_pairwise', 'beta_fd_pairwise_pa', 'beta_fd_nt', 'beta_fd_nt_pa')
all_vars <- c('beta_td_pairwise_pa', 'beta_td_sorensen_pa',
							'beta_td_pairwise', 'beta_td_sorensen',
							'beta_td_shannon', 
							'beta_pd_pairwise_pa', 'beta_pd_pairwise_pa_z',
							'beta_pd_nt_pa', 'beta_pd_nt_pa_z',
							'beta_pd_pairwise', 'beta_pd_pairwise_z',
							'beta_pd_nt', 'beta_pd_nt_z',
							'beta_fd_pairwise_pa', 'beta_fd_pairwise_pa_z',
							'beta_fd_nt_pa', 'beta_fd_nt_pa_z',
							'beta_fd_pairwise', 'beta_fd_pairwise_z',
							'beta_fd_nt', 'beta_fd_nt_z')
metric_name <- all_vars[metric]
is_prop <- metric_name %in% prop_vars

load('/mnt/home/qdr/data/fiaworkspace_spatial_wholeusa_2018.r')

###############################################
# subset to keep only the natural plots
plantation <- read.csv('/mnt/research/nasabio/data/fia/plotcond/plantation.csv')
fiacoords <- left_join(fiacoords, plantation)

###############################################

# Load the correct metric.
print(metric)
load(paste0('/mnt/research/nasabio/data/fia/diversity/usa/metric_betatd100_', metric, '.r'))


# For each year and route number, get the median beta diversity within each radius.
# The pairwise table was constructed to only go up to 100 km, so that is all we will have.
radii <- c(5, 10, 20, 50, 75, 100) # in km

library(sp)
library(purrr)

# Function with 2 arguments. x is element i of metric_list, while y is row i of fiacoords
neighbordivfromlist <- function(x, y, is_prop) {
	neighbors <- fiacoords[x[,1],]
	neighbordists <- spDistsN1(pts = cbind(neighbors$lon, neighbors$lat), pt = c(y$lon, y$lat), longlat = TRUE)
	commdat <- list()
	for (i in 1:length(radii)) {
		beta_incircle <- x[,2][neighbordists <= radii[i] & !neighbors$plantation]
		beta_incircle <- beta_incircle[is.finite(beta_incircle)]
		if (is_prop) {
			beta_mean <- sin(mean(asin(sqrt(beta_incircle))))^2
		} else {
			beta_mean <- mean(beta_incircle)
		}
		commdat[[i]] <- c(PLT_CN = y$PLT_CN, radius = radii[i], beta = beta_mean)
						  
	}
	as.data.frame(do.call('rbind', commdat))
}

fia_beta <- map2(metric_list, split(fiacoords, 1:nrow(fiacoords)), neighbordivfromlist, is_prop = is_prop)
fia_beta <- bind_rows(fia_beta)
names(fia_beta)[3] <- metric_name

write.csv(fia_beta, file = paste0('/mnt/research/nasabio/data/fia/diversity/usa2018/beta100_metric_', metric, '.csv'), row.names = FALSE)

