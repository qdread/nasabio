# FIA taxonomic and phylogenetic beta diversity calculation for the cluster.

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/'
fiapnw <- read.csv(file.path(fp, 'finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)

library(ape)
load(file.path(fp, 'phytophylo_fia.r'))
fiadist <- cophenetic(fiaphytophylo)

fiataxa <- read.csv('specieslists/fia_taxon_lookuptable.csv', stringsAsFactors = FALSE)
pnw_codes <- unique(fiapnw$SPCD)
all(pnw_codes %in% fiataxa$FIA.Code) #yes

# Scientific names of FIA species in PNW plots.
pnw_species <- fiataxa[match(pnw_codes, fiataxa$FIA.Code), c('Genus','Species')]
pnw_scinames <- paste(pnw_species$Genus, pnw_species$Species, sep = '_')
pnw_scinames <- gsub(' ', '', pnw_scinames) #Remove extraneous spaces

# Correction for subspecies that are not in the tree.
pnw_species$sciname <- pnw_scinames
pnw_species$sciname[pnw_species$sciname == 'Chrysolepis_chrysophyllavar.chrysophylla'] <- 'Chrysolepis_chrysophylla'
pnw_species$sciname[pnw_species$sciname == 'Populus_balsamiferassp.Trichocarpa'] <- 'Populus_trichocarpa'

pnw_species <- subset(pnw_species, !grepl('Tree', sciname)) # unidentified species
#pnw_species$sciname <- gsub('_', ' ', pnw_species$sciname)
pnw_species <- pnw_species[order(pnw_species$sciname), ]

library(dplyr)

radii <- c(1000,2000,3000,4000,5000,10000)
# Initialize data structures for observed metrics
fia_shannonbetadiv <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_meanpairwisedissim <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_meanpairwisedissim_pa <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_phypairwise <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_phypairwise_pa <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_phynt <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_phynt_pa <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_phypairwise_z <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_phypairwise_pa_z <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_phynt_z <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_phynt_pa_z <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))

fia_nneighb <- matrix(0, nrow = nrow(fiaalbers), ncol = length(radii))


pb2 <- txtProgressBar(0, length(radii) * nrow(fiaalbers), style = 3)
i <- 0

# Number of simulations for null model
nnull <- 999

for (r in 1:length(radii)) {
  fianhb_r <- getNeighbors(fiaalbers, radius = radii[r])
  for (p in 1:nrow(fiaalbers)) {
    i <- i + 1
    if (class(fianhb_r[[p]]) == 'data.frame') {
      # Subset out the data frame with the nearest neighbors
      plotcns <- fiaalbers[c(p, fianhb_r[[p]]$idx), ]$PLT_CN
      dat_p <- subset(fiasums_plot, PLT_CN %in% plotcns)
      # Convert into a site x species matrix
      sppids <- sort(unique(dat_p$SPCD))
      mat_p <- dat_p %>% group_by(PLT_CN) %>% do(x = sapply(sppids, function(z) sum(.$basalarea[.$SPCD == z])))
      mat_p <- do.call('rbind', mat_p$x)
      
      if(!is.null(mat_p)) {
        if(nrow(mat_p) > 1) {
      # Fix the species names to match the phylogeny, and get rid of the unknown species.
      sppnames <- pnw_species$sciname[match(sppids, pnw_codes)]
      dimnames(mat_p)[[1]] <- 1:nrow(mat_p)
      dimnames(mat_p)[[2]] <- sppnames
      mat_p <- mat_p[, dimnames(mat_p)[[2]] %in% fiaphytophylo$tip.label, drop = FALSE]
      
      
      # Calculate beta-diversity for that matrix.
      
        fia_shannonbetadiv[p, r] <- d(abundances = mat_p, lev = 'beta', wts = FALSE, q = 1)
        fia_meanpairwisedissim[p, r] <- mean(vegdist(x = mat_p, binary = FALSE, method = 'jaccard'))
        fia_meanpairwisedissim_pa[p, r] <- mean(vegdist(x = mat_p, binary = TRUE, method = 'jaccard'))
        
        # Must catch errors in comdist for when the number of columns is zero
        if (ncol(mat_p) > 1) {
        fia_phypairwise[p, r] <- mean(comdist(comm = mat_p, dis = fiadist, abundance.weighted = TRUE))
        fia_phypairwise_pa[p, r] <- mean(comdist(comm = mat_p, dis = fiadist, abundance.weighted = FALSE))
        fia_phynt[p, r] <- mean(comdistnt(comm = mat_p, dis = fiadist, abundance.weighted = TRUE))
        fia_phynt_pa[p, r] <- mean(comdistnt(comm = mat_p, dis = fiadist, abundance.weighted = FALSE))
        }
        else {
          fia_phypairwise[p, r] <- 0
          fia_phypairwise_pa[p, r] <- 0
          fia_phynt[p, r] <- 0
          fia_phynt_pa[p, r] <- 0
        }
        fia_nneighb[p, r] <- nrow(mat_p) - 1
        
        # Null models by scrambling distance matrix
        phypairwise_null <- phypairwise_pa_null <- phynt_null <- phynt_pa_null <- rep(NA, nnull)
        
        for (sim in 1:nnull) {
          nullidx <- sample(1:nrow(fiadist))
          fiadistnull <- fiadist
          dimnames(fiadistnull)[[1]] <- dimnames(fiadistnull)[[1]][nullidx]
          dimnames(fiadistnull)[[2]] <- dimnames(fiadistnull)[[2]][nullidx]
          
          phypairwise_null[sim] <- mean(comdist(comm = mat_p, dis = fiadistnull, abundance.weighted = TRUE))
          phypairwise_pa_null[sim] <- mean(comdist(comm = mat_p, dis = fiadistnull, abundance.weighted = FALSE))
          phynt_null[sim] <- mean(comdistnt(comm = mat_p, dis = fiadistnull, abundance.weighted = TRUE))
          phynt_pa_null[sim] <- mean(comdistnt(comm = mat_p, dis = fiadistnull, abundance.weighted = FALSE))
        }
        
        fia_phypairwise_z[p, r] <- (fia_phypairwise[p,r] - mean(phypairwise_null, na.rm=T))/sd(phypairwise_null, na.rm=T)
        fia_phypairwise_pa_z[p, r] <- (fia_phypairwise_pa[p,r] - mean(phypairwise_pa_null, na.rm=T))/sd(phypairwise_pa_null, na.rm=T)
        fia_phynt_z[p, r] <- (fia_phynt[p,r] - mean(phynt_null, na.rm=T))/sd(phynt_null, na.rm=T)
        fia_phynt_pa_z[p, r] <- (fia_phynt_pa[p,r] - mean(phynt_pa_null, na.rm=T))/sd(phynt_pa_null, na.rm=T)
        }
      }
    }
    setTxtProgressBar(pb2, i)
  }
}

close(pb2)

save(fia_shannonbetadiv, fia_nneighb, file = file.path(fp, 'fia_taxonomicbetadiv.RData'))
