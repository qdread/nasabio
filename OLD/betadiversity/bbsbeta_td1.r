# Run on small chunk.

r <- 1

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

load(paste0('/mnt/research/nasabio/data/bbs/mats/mat_',task,'_',r,'.r',sep='_'))

bbs_meanpairwisedissim_pa <- rep(NA, length(x))
bbs_nneighb <- rep(NA, length(x))

pb2 <- txtProgressBar(0, length(x), style = 3)

library(vegan)
library(vegetarian)

for (p in 1:length(x)) {
   
    if(!is.null(x[[p]])) {
      if(nrow(x[[p]]) > 1) {
		
		mat_p <- x[[p]]
		
        dimnames(mat_p)[[1]] <- 1:nrow(mat_p)
        
        # Calculate beta-diversity for that matrix.
        bbs_meanpairwisedissim_pa[i] <- mean(vegdist(x = mat_p, binary = TRUE, method = 'jaccard'))
        bbs_nneighb[i] <- nrow(mat_p) - 1
        
        
      }
    }
  
  setTxtProgressBar(pb2, p)
}

close(pb2)

# Compile all of these values into a single data frame and save.
bbs_betadiv <- data.frame(nneighb = bbs_nneighb,
						  beta_td_pairwise_presence = bbs_meanpairwisedissim_pa
						  )

write.csv(bbs_betadiv, file = paste0('/mnt/research/nasabio/data/bbs/bbs_betatd_',task,'_',r,'.csv'), row.names = FALSE)			