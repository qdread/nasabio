### TOKEN FOR ACCESSING IUCN API ###########################################
### This is only authorized to be used by people from MSU ##################
token <- '3b3db1c753a7cad0616e16146363ec88eead97d185cb3304070d7343123516fd'
############################################################################


# Names of all threatened species in the USA ------------------------------



library(httr)
library(jsonlite)
library(dplyr)

api_url <- 'apiv3.iucnredlist.org/api/v3'

# Get all species in the United States.

us_spp <- GET(url = paste0(api_url, '/country/getspecies/US?token=', token))
us_spp <- fromJSON(content(us_spp, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)$result

# As of 10 May 2018, this gives 10004 species (including all taxa and least concern)

# How many non least-concern species are there in the USA?
table(us_spp$category) # A little over 3000 are in some way threatened.
us_threatened <- subset(us_spp, !category %in% c('LC'))


# Break down by species group if possible ---------------------------------



# Get codes for the different species groups
grp_codes <- GET(url = paste0(api_url, '/comp-group/list?token=', token))
fromJSON(content(grp_codes, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)

# Relevant codes: "birds" "conifers" "magnolias" "mangrove_plants"

# Get species list of all birds across all countries
all_birds <- GET(url = paste0(api_url, '/comp-group/getspecies/birds?token=', token))
all_birds <- fromJSON(content(all_birds, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)$result

# Get all the ones that could possibly be trees
# There is not really any way to get the full list of plants.
conifers <- GET(url = paste0(api_url, '/comp-group/getspecies/conifers?token=', token))
conifers <- fromJSON(content(conifers, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)$result
magnolias <- GET(url = paste0(api_url, '/comp-group/getspecies/magnolias?token=', token))
magnolias <- fromJSON(content(magnolias, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)$result



