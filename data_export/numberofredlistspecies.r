# IUCN Red List

library(XLConnect)
birdredlist <- readWorksheetFromFile('C:/Users/Q/google_drive/NASABiodiversityWG/Trait_Data/redlist.xlsx', sheet='birds')
treeredlist <- readWorksheetFromFile('C:/Users/Q/google_drive/NASABiodiversityWG/Trait_Data/redlist.xlsx', sheet='trees')

# Determine which BBS routes and FIA plots contain these species.
fp <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Trait_Data'
write.csv(birdredlist, file=file.path(fp,'birdredlist.csv'), row.names = FALSE)
write.csv(treeredlist, file=file.path(fp,'treeredlist.csv'), row.names = FALSE)

# FIA

fiataxa <- read.csv('~/FIA/FIA10nov/lookup_table_allfia.csv', stringsAsFactors = FALSE)

treeredlist$Name <- gsub('var.', '', treeredlist$Name)
treeredlist$Name <- gsub('ssp.', '', treeredlist$Name)
treeredlist$Name <- gsub('\\s+', '_', treeredlist$Name)

fia_iucn <- treeredlist$Name[treeredlist$Name %in% fiataxa$binomial | treeredlist$Name %in% fiataxa$binomial_forphylo]

fia_iucn_table <- fiataxa[fiataxa$binomial %in% fia_iucn,]

# BBS

bbstaxa <- read.csv('~/GitHub/aquaxterra/data/specieslist.csv', stringsAsFactors = FALSE)

bbs_iucn <- birdredlist$Name[birdredlist$Name %in% bbstaxa$Latin_Name_clean]

bbs_iucn_table <- bbstaxa[bbstaxa$Latin_Name_clean %in% bbs_iucn,]

fp <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Trait_Data'
write.csv(bbs_iucn_table, file=file.path(fp,'birdredlist_inbbs.csv'), row.names = FALSE)
write.csv(fia_iucn_table, file=file.path(fp,'treeredlist_infia.csv'), row.names = FALSE)

### On cluster

#FIA
load('/mnt/research/nasabio/data/fia/fiaworkspace_nospatial_wholeusa.r') # Very large file.
fiaiucn <- read.csv('/mnt/research/nasabio/data/fia/treeredlist_infia.csv', stringsAsFactors = FALSE)

# Look at how many redlisted species are in each FIA plot
library(dplyr)

fiaiucnmat <- fiaplotmat[, dimnames(fiaplotmat)[[2]] %in% fiaiucn$binomial]
fiaiucnsum <- apply(fiaiucnmat>0, 1, sum)

fiacoords <- read.csv('~/data/allfia.csv')

write.csv(data.frame(PLT_CN = fiacoords$CN, n_IUCN = fiaiucnsum), file = '/mnt/research/nasabio/data/fia/IUCN_by_plot.csv', row.names = FALSE)

# Number of plots by species
byplottable <- apply(fiaiucnmat>0, 2, sum)
write.csv( data.frame(species=names(byplottable), n_plots=byplottable), file = '/mnt/research/nasabio/data/fia/redlist_FIAplotsperspecies.csv', row.names= FALSE )

#BBS
load('/mnt/research/nasabio/data/bbs/bbsworkspace_byroute.r')
bbsiucn <- read.csv('/mnt/research/nasabio/data/bbs/birdredlist_inbbs.csv', stringsAsFactors = FALSE)

bbsiucnmat <- fixedbbsmat_byroute[, dimnames(fixedbbsmat_byroute)[[2]] %in% bbsiucn$Latin_Name_clean]
bbsiucnsum <- cbind(bbscov, n_IUCN = apply(bbsiucnmat>0, 1, sum))
bbsiucnsum <- subset(bbsiucnsum, year >= 1997)

write.csv(bbsiucnsum, file = '/mnt/research/nasabio/data/bbs/IUCN_by_plot.csv', row.names = FALSE)

byyeartable <- bbsiucnsum %>% group_by(year) %>% do(data.frame(n_spp = 0:7, n_routes = sapply(0:7, function(x) sum(.$n_IUCN == x))))

write.csv(byyeartable, file = '/mnt/research/nasabio/data/bbs/redlist_BBSroutes.csv', row.names = FALSE)

# For each species, how many routes does it appear in each year?
byspeciestable <- bbscov %>%
  cbind(bbsiucnmat) %>%
  filter(year >= 1997) %>%
  group_by(year) %>%
  summarize_at(-(1:6), funs(sum(.>0)))

write.csv(byspeciestable, file = '/mnt/research/nasabio/data/bbs/redlist_BBSroutesperspecies.csv', row.names = FALSE)