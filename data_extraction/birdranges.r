# Range sizes from the Birds of the World dataset. Calculated them in QGIS, now load into R
fp <- 'C:/Users/Q/Dropbox/projects/verts/'
botw <- read.csv(file.path(fp,'botw_table.csv'), stringsAsFactors = FALSE)

# Should be able to get the migrant status out of this too.
botw <- mutate(botw, taxon = gsub(' ', '_', SCINAME))

library(dplyr)

rangedat <- function(x) {
  migrant_status <- 'none'
  if (any(x$SEASONAL %in% 2:3)) migrant_status <- 'partial'
  if (!any(x$SEASONAL == 1)) migrant_status <- 'obligate'
  data.frame(migrant_status = migrant_status, range_size = sum(x$Area_km2[x$SEASONAL %in% 1:2]))
}

botw_ranges_all <- botw %>%
  filter(PRESENCE == 1) %>% # Only where species is known present.
  group_by(taxon) %>%
  do(rangedat(.))

write.csv(botw_ranges_all, 'C:/Users/Q/Dropbox/projects/nasabiodiv/birdranges.csv', row.names = FALSE)

# Load the previous bird trait matrix and add migrant status if possible
brange <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/birdranges.csv', stringsAsFactors = FALSE)
brange$taxon <- gsub('_', ' ', brange$taxon)
bt <- read.csv('C:/Users/Q/Dropbox/projects/aquaxterra/birdtraitmerged.csv', stringsAsFactors = FALSE)

library(dplyr)
#bt <- left_join(bt, brange %>% rename(Latin_Name_clean = taxon))
#bt <- left_join(bt, brange %>% rename(Latin_Name_synonym = taxon), by = 'Latin_Name_synonym')
#bt <- left_join(bt, brange %>% rename(Latin_Name_synonym2 = taxon), by = 'Latin_Name_synonym2')

#mstatus <- pmin(bt$migrant_status, bt$migrant_status.x, bt$migrant_status.y, na.rm = T)

# Fish and Wildlife's list of migratory birds
migrants <- read.csv('C:/Users/Q/Documents/GitHub/nasabio/specieslists/migratorybirds.csv', stringsAsFactors = F)

# Get rid of asterisks in the Latin name strings, and spaces at the end. 
trim.trailing <- function (x) sub("\\s+$", "", x)
migrants$Latin.name <- sub('\\*', '', migrants$Latin.name)
migrants$Latin.name <- trim.trailing(migrants$Latin.name)

migrants <- subset(migrants, Latin.name != '')
migrants$Latin.name %in% bt$Latin_Name_clean
migrants$Latin.name[!migrants$Latin.name %in% bt$Latin_Name_clean]
migrants$Latin.name[!migrants$Latin.name %in% bt$Latin_Name_synonym]

ismigrant <- bt$Latin_Name_clean %in% migrants$Latin.name | bt$Latin_Name_synonym %in% migrants$Latin.name | bt$Latin_Name_synonym2 %in% migrants$Latin.name

bt$migrant_status <- ismigrant
write.csv(bt, 'C:/Users/Q/Dropbox/projects/aquaxterra/birdtraitmerged.csv', row.names = FALSE)
