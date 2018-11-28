# Compile FIA gamma diversity (unfuzzed)
# 26 Oct 2017
# Edited 31 Dec for whole USA
# Edited 28 Nov 2018: update all file paths for updated FIA diversity metrics with macroplots removed.

gamma_div_list <- list()

for (i in 1:1000) {
  load(paste0('/mnt/research/nasabio/data/fia/diversity/usa2018/gamma_', i, '.r'))
  gamma_div_list[[i]] <- gamma_div
}

# Bind into one large array then convert into a data frame.
library(abind)
gamma_div_array <- abind(gamma_div_list, along = 1)

# Stack into a very long data frame then recast
library(dplyr)
library(reshape2)
gamma_div_melt <- melt(gamma_div_array, varnames = c('plot', 'radius', 'diversity_type'))
gamma_div_melt$radius <- as.numeric(substr(as.character(gamma_div_melt$radius), start = 3, stop = nchar(as.character(gamma_div_melt$radius))))
gamma_div_cast <- dcast(gamma_div_melt, formula = plot + radius ~ diversity_type)

# Add identifying information
plotmetadata <- read.csv('/mnt/research/nasabio/data/fia/fianocoords_wholeusa_2018.csv')
plantation <- read.csv('/mnt/research/nasabio/data/fia/plotcond/plantation.csv')
plotmetadata <- left_join(plotmetadata, plantation) %>% filter(!plantation)

gamma_div_cast <- cbind(PLT_CN = plotmetadata$PLT_CN[rep(1:nrow(plotmetadata), each = dim(gamma_div_array)[2])], gamma_div_cast[, !names(gamma_div_cast) %in% 'plot'])
write.csv(gamma_div_cast, file = '/mnt/research/nasabio/data/fia/biodiversity_CSVs/updated_nov2018/fiausa_natural_gamma.csv', row.names = FALSE)
