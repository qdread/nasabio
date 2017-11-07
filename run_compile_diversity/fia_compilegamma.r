# Compile FIA gamma diversity (unfuzzed)
# 26 Oct 2017

gamma_div_list <- list()

for (i in 1:100) {
  load(paste0('/mnt/research/nasabio/data/fia/diversity/unfuzzed/gamma_', i, '.r'))
  gamma_div_list[[i]] <- gamma_div
}

# Bind into one large array then convert into a data frame.
library(abind)
gamma_div_array <- abind(gamma_div_list, along = 1)

# Stack into a very long data frame then recast
library(reshape2)
gamma_div_melt <- melt(gamma_div_array, varnames = c('plot', 'radius', 'diversity_type'))
gamma_div_melt$radius <- as.numeric(substr(as.character(gamma_div_melt$radius), start = 3, stop = nchar(as.character(gamma_div_melt$radius))))
gamma_div_cast <- dcast(gamma_div_melt, formula = plot + radius ~ diversity_type)

# Add identifying information
plotmetadata <- read.csv('/mnt/research/nasabio/data/fia/fianocoords.csv')

gamma_div_cast <- cbind(plotmetadata, gamma_div_cast[, !names(gamma_div_cast) %in% 'plot'])
write.csv(gamma_div_cast, file = '/mnt/research/nasabio/data/fia/fia_gammadiv.csv')
