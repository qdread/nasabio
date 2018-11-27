# in each slice, load calculated fia matrix.
# Calculate gamma diversity within different radii.

# Split into 1000 groups.
# Edited 21 Dec for whole USA
# Edited 22 Dec not to use precalculated matrices
# New version created on 28 Feb 2018: filter the plots by whether they are plantation or not, in order to calculate the gammas.
# Edited 26 Nov 2018: file paths update, update for slurm
# Edited 27 Nov 2018: redo plantation subsetting to fix the mismatched ordering. Better coding practice to not assume the two are both sorted the same way, anyway.

load('/mnt/ffs17/groups/nasabio/fiaworkspace_nospatial_wholeusa_2018.r')
load('/mnt/home/qdr/data/fiaworkspace_spatial_wholeusa_2018.r')
source('/mnt/research/nasabio/code/pairwise_beta_focal.r')
source('/mnt/research/nasabio/code/nofuncspp.r')

library(sp)
library(vegan)
source('/mnt/research/nasabio/code/fixpicante.r')

nnull <- 99
trydist <- as.matrix(trydist)

radii <- c(5, 10, 20, 50, 75, 100, 150, 200, 300)
n_slices <- 1000
slice <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

#################################
# Subset fiaplotmat and fiacoords to keep only "natural" plots.

plantation <- read.csv('/mnt/research/nasabio/data/fia/plotcond/plantation.csv')
natural_plots <- plantation$PLT_CN[!plantation$plantation]
natural_plot_idx <- fiacoords$PLT_CN %in% natural_plots
fiacoords <- fiacoords[natural_plot_idx, ]
fiaplotmat <- fiaplotmat[natural_plot_idx, ]


#################################

# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,nrow(fiaplotmat),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

# Declare structures to hold data
# Edit 25 Oct: make each slice smaller by getting rid of some of the NAs.
pb <- txtProgressBar(rowidxmin, rowidxmax, style = 3)
gamma_div <- array(NA, dim = c(length(rowidxmin:rowidxmax), length(radii), 11))

cnames <- c('richness', 'shannon', 'evenness',
            'MPD_pa_z', 'MNTD_pa_z',
            'MPD_z', 'MNTD_z',
            'MPDfunc_pa_z', 'MNTDfunc_pa_z',
            'MPDfunc_z', 'MNTDfunc_z')

null_result <- rep(NA, length(cnames))

for (p in rowidxmin:rowidxmax) {
  setTxtProgressBar(pb, p)
  # Distances between target plot and all other plots.
	dist_p <- spDistsN1(pts=with(fiacoords, cbind(lon, lat)), pt = c(fiacoords$lon[p], fiacoords$lat[p]), longlat = TRUE)
	
  for (r in 1:length(radii)) {
	
	neighbs <- fiaplotmat[dist_p <= radii[r], ]
	
    gamma_div[p - rowidxmin + 1, r, ] <- tryCatch(diversity_3ways(m = neighbs, flavor = 'gamma', 
                                                                  dotd = TRUE, dopd = TRUE, dofd = TRUE, abundance = TRUE,
                                                                  pddist = fiadist, fddist = trydist,
                                                                  nnull = nnull,
                                                                  phylo_spp = fullphylo$tip.label, func_problem_spp = nofuncspp),
                                                  error = function(e) null_result)	
  }
  
}

close(pb)

dimnames(gamma_div)[[3]] <- cnames
dimnames(gamma_div)[[2]] <- paste('r',radii,sep='_')

save(gamma_div, file = paste0('/mnt/research/nasabio/data/fia/diversity/usa2018/gamma_', slice, '.r'))
