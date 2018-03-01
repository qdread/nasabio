files <- paste0('https://apps.fs.usda.gov/fia/datamart/CSV/', state.abb, '_COND.zip')
for (file in files) system2('wget', args = file)
for (file in paste0(state.abb, '_COND.zip')) system2('unzip', args = file)
conds <- lapply(paste0(state.abb, '_COND.csv'), read.csv, stringsAsFactors = FALSE)
condsdf <- do.call(rbind, conds)
write.csv(condsdf, file = 'site_cond.csv', row.names = FALSE)

condsdf <- read.csv('/mnt/research/nasabio/data/fia/plotcond/site_cond.csv', stringsAsFactors = FALSE)

library(dplyr)
plantation <- condsdf %>%
  group_by(PLT_CN) %>%
  summarize(plantation = any(STDORGCD == 1))

source('/mnt/research/nasabio/code/loadfiaall.r')

table(fiacoords$PLT_CN %in% plantation$PLT_CN) # All are accounted for.
fiacoords <- left_join(fiacoords, plantation)
plantation$plantation[is.na(plantation$plantation)] <- FALSE

write.csv(fiacoords[,c('PLT_CN','plantation')], file = '/mnt/research/nasabio/data/fia/plotcond/plantation.csv', row.names = FALSE)