### TOKEN FOR ACCESSING IUCN API ###########################################
### This is only authorized to be used by people from MSU ##################
token <- '3b3db1c753a7cad0616e16146363ec88eead97d185cb3304070d7343123516fd'
############################################################################


# Names of all threatened species in CO/EC ------------------------------

library(httr)
library(jsonlite)
library(dplyr)

api_url <- 'apiv3.iucnredlist.org/api/v3'

# Get all species in Colombia and Ecuador

co_spp <- GET(url = paste0(api_url, '/country/getspecies/CO?token=', token))
co_spp <- fromJSON(content(co_spp, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)$result

ec_spp <- GET(url = paste0(api_url, '/country/getspecies/EC?token=', token))
ec_spp <- fromJSON(content(ec_spp, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)$result

# Break down by species group if possible ---------------------------------

# Get codes for the different species groups
grp_codes <- GET(url = paste0(api_url, '/comp-group/list?token=', token))
fromJSON(content(grp_codes, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)

# Get species list of all birds across all countries
all_birds <- GET(url = paste0(api_url, '/comp-group/getspecies/birds?token=', token))
all_birds <- fromJSON(content(all_birds, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)$result

# Get species list of all mammals across all countries
all_mammals <- GET(url = paste0(api_url, '/comp-group/getspecies/mammals?token=', token))
all_mammals <- fromJSON(content(all_mammals, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)$result


# Mammals and birds occurring in CO and EC --------------------------------

co_ec_IDs <- union(co_spp$taxonid, ec_spp$taxonid)
all_spp <- rbind(all_mammals, all_birds)
co_ec_spp <- subset(all_spp, taxonid %in% co_ec_IDs)


# Get habitat for each species --------------------------------------------

# Include error-checking code in the function so that it will return NA for everything if there's no data.

get_habitat <- function(ID) {
  habitat_by_sp <- GET(url = paste0(api_url, '/habitats/species/id/', ID,'?token=', token))
  habitat_by_sp <- fromJSON(content(habitat_by_sp, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)$result
  if (class(habitat_by_sp) == 'data.frame') {
    data.frame(taxonid=ID, habitat_by_sp)
  } else {
    data.frame(taxonid=ID)
  }
}

# Apply get_habitat to each row of co_ec_spp
all_habitat <- co_ec_spp %>%
  rowwise %>%
  do(get_habitat(.$taxonid))

write.csv(all_habitat, 'habitat_by_species.csv', row.names = FALSE)

# For testing purposes, run on 1st 50 rows
# all_habitat <- co_ec_spp[1:50,] %>%
#   rowwise %>%
#   do(get_habitat(.$taxonid))
