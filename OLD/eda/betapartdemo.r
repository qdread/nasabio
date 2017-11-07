# source('~/GitHub/NEON/code/baselga.r')
# 
# betadiv_baselga <- lapply(comm_bybout, function(x) {
#   ddply(x, .(year(date)), function(z) {
#     1/as.numeric(Simpson.multi(z[,-(1:2)] > 0)[2])
#   })         
# })

library(betapart)
data(bbsData)

beta.multi(bbs2000, index.family = 'sorensen')
beta.multi.abund(bbs2000, index.family = 'bray')
bbspairdmat <- beta.pair(bbs2000, index.family = 'sorensen')
bbspairdmat_abund <- beta.pair.abund(bbs2000, index.family = 'bray')

# Can do pairwise for single pairs
beta.pair(bbs2000[c(14,17,40),], index.family = 'sorensen')
beta.pair.abund(bbs2000[c(14,17),], index.family = 'bray')

# Resampling in subset of sites? To get confidence interval
beta.sample(bbs2000, index.family = 'sorensen', sites=3, samples = 10)
beta.sample.abund(bbs2000, index.family = 'bray', sites=3, samples = 10)

# This also has functional and phylogenetic diversity!
# functional.beta.multi(x, traits, index.family = 'sorensen', warning.time=FALSE)