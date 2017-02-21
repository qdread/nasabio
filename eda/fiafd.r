# Functional diversity of FIA

# First get list of species IDs for TRY.
tryspp <- read.csv(file.path(fp, 'tryspp.csv'), stringsAsFactors = FALSE)

trymatch <- pnw_species$sciname %in% tryspp$AccSpeciesName
pnw_species$sciname[!trymatch] # All match except for the "unknowns"

tryids <- tryspp$AccSpeciesID[match(pnw_species$sciname, tryspp$AccSpeciesName)]
write.table(t(na.omit(tryids)), sep=',')
