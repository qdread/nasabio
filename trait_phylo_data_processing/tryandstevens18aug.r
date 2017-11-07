# Combine FIA traits from the new TRY data pull with older traits from Stevens' trait dataset.
# Impute any missing values so that we can use it with the full FIA dataset, no missing species.

# Load TRY

trait_try <- read.csv('C:/Users/Q/google_drive/NASABiodiversityWG/Trait_Data/fia_try_17aug/try_trait_byspecies.csv', stringsAsFactors = FALSE)

# Process SLA related traits
slanames <- grep('SLA', names(trait_try), value=T)
slatraits <- trait_try[,slanames]
# Remove outliers.
slatraits[slatraits > 200] <- NA

SLA <- apply(slatraits, 1, mean, na.rm=TRUE)

library(dplyr)
trait_try_use <- trait_try %>% select(AccSpeciesName, Bark.thickness, Plant.lifespan.years, Root.rooting.depth, Seed.dry.mass, Stem.dry.mass.per.stem.fresh.volume..stem.specific.density..SSD..wood.density.) %>%
  rename(Scientific_Name = AccSpeciesName, Plant.lifespan=Plant.lifespan.years, Rooting.depth=Root.rooting.depth, SSD=Stem.dry.mass.per.stem.fresh.volume..stem.specific.density..SSD..wood.density.) %>%
  cbind(SLA) %>%
  mutate(Scientific_Name = gsub('\\ ', '_', Scientific_Name))

# Load Stevens

library(XLConnect)
library(dplyr)
trait_stevens <- readWorksheetFromFile('C:/Users/Q/google_drive/NASABiodiversityWG/Trait_Data/Traits_Stevens_FIA.xlsx', sheet = 'Master')

# Convert to numerics where needed.
numeric_cols <- c(2,5:ncol(trait_stevens))
trait_stevens[, numeric_cols] <- lapply(trait_stevens[,numeric_cols], as.numeric)

# Get rid of subspecies and genus-level trait values
trait_stevens <- filter(trait_stevens, Scientific_Name != 'Tree_unknown', Subsp %in% c(0,0.5), Generic != 1)



# Phylogeny linked with these
library(ape)
allfia_phylogeny <- read.tree('C:/Users/Q/Dropbox/projects/nasabiodiv/allfiaphylogeny/tree_all_final_031716.txt')

# How many of Stevens scientific names are in there.
spmatch <- trait_stevens$Scientific_Name %in% allfia_phylogeny$tip.label
table(spmatch)
trait_stevens$Scientific_Name[!spmatch]

# Correct names that are typos or obsolete names in trait_stevens
name_correction <- c('Chamaecyparis_lawsonia' = 'Chamaecyparis_lawsoniana',
                     'Aesculus_octandra' = 'Aesculus_flava',
                     "Carya_illinoensis" = "Carya_illinoinensis" ,
                     'Carya_tomentosa' = 'Carya_alba',
                     'Populus_trichocarpa' = "Populus_balsamifera_trichocarpa",
                     "Quercus_chrysolepsis" = 'Quercus_chrysolepis',
                     'Quercus_nutallii' = 'Quercus_texana',
                     'Quercus_wislizenii' = 'Quercus_wislizeni',
                     'Tilia_heterophylla' = 'Tilia_americana_heterophylla',
                     'Ulmus_pumilia' = 'Ulmus_pumila')

trait_stevens$Scientific_Name[match(names(name_correction), trait_stevens$Scientific_Name)] <- name_correction

# 3 still don't match
spmatch <- trait_stevens$Scientific_Name %in% allfia_phylogeny$tip.label
table(spmatch)
trait_stevens <- filter(trait_stevens, spmatch)

# Get the traits we want only.
trait_stevens_use <- trait_stevens %>% select(Scientific_Name, Bark.thickness, SLA, Plant.lifespan, Seed.dry.mass, SSD)


# Combine TRY and Stevens
table(trait_try_use$Scientific_Name %in% trait_stevens_use$Scientific_Name)

trait_trystevens_use <- full_join(trait_try_use, trait_stevens_use, by = "Scientific_Name")

# Use newer TRY value for other traits. If NA, then use Stevens value.
trait_trystevens_use$Plant.lifespan.x[is.na(trait_trystevens_use$Plant.lifespan.x)] <- trait_trystevens_use$Plant.lifespan.y[is.na(trait_trystevens_use$Plant.lifespan.x)]
trait_trystevens_use$SSD.x[is.na(trait_trystevens_use$SSD.x)] <- trait_trystevens_use$SSD.y[is.na(trait_trystevens_use$SSD.x)]
trait_trystevens_use$SLA.x[is.na(trait_trystevens_use$SLA.x)] <- trait_trystevens_use$SLA.y[is.na(trait_trystevens_use$SLA.x)]
trait_trystevens_use$Seed.dry.mass.x[is.na(trait_trystevens_use$Seed.dry.mass.x)] <- trait_trystevens_use$Seed.dry.mass.y[is.na(trait_trystevens_use$Seed.dry.mass.x)]


trait_all_use <- with(trait_trystevens_use, data.frame(Scientific_Name=Scientific_Name,
                                                       Bark.thickness = Bark.thickness.y,
                                                       SLA = SLA.x,
                                                       SSD = SSD.x, 
                                                       Seed.dry.mass = Seed.dry.mass.x,
                                                       Rooting.depth = Rooting.depth,
                                                       Plant.lifespan = Plant.lifespan.x))

########################################################
# Get rid of everything that isn't in the western plots.

pnwspp <- c("Abies_amabilis","Abies_concolor","Abies_grandis","Abies_lasiocarpa","Abies_magnifica","Abies_procera","Chamaecyparis_lawsoniana","Chamaecyparis_nootkatensis","Cupressus_bakeri","Cupressus_sargentii","Juniperus_californica","Juniperus_occidentalis","Juniperus_osteosperma","Juniperus_scopulorum","Larix_lyallii","Larix_occidentalis","Calocedrus_decurrens","Picea_breweriana","Picea_engelmannii","Picea_glauca","Picea_mariana","Picea_sitchensis","Pinus_albicaulis","Pinus_attenuata","Pinus_balfouriana","Pinus_contorta","Pinus_coulteri","Pinus_flexilis","Pinus_jeffreyi","Pinus_lambertiana","Pinus_monticola","Pinus_muricata","Pinus_ponderosa","Pinus_radiata","Pinus_sabiniana","Pinus_sylvestris","Pinus_monophylla","Pinus_washoensis","Pinus_longaeva","Pseudotsuga_macrocarpa","Pseudotsuga_menziesii","Sequoia_sempervirens","Sequoiadendron_giganteum","Taxus_brevifolia","Thuja_plicata","Torreya_californica","Tsuga_heterophylla","Tsuga_mertensiana","Acer_macrophyllum","Acer_negundo","Acer_platanoides","Acer_glabrum","Aesculus_sp","Aesculus_californica","Ailanthus_altissima","Alnus_rubra","Alnus_rhombifolia","Arbutus_menziesii","Betula_occidentalis","Betula_papyrifera","Chrysolepis_chrysophylla","Cercocarpus_ledifolius","Cornus_nuttallii","Eucalyptus_globulus","Fraxinus_sp","Fraxinus_latifolia","Juglans_hindsii","Juglans_californica","Liquidambar_styraciflua","Lithocarpus_densiflorus","Malus_fusca","Platanus_racemosa","Populus_tremuloides","Populus_balsamifera_trichocarpa","Populus_fremontii","Prosopis_glandulosa","Prosopis_pubescens","Prunus_virginiana","Prunus_emarginata","Prunus_avium","Quercus_agrifolia","Quercus_chrysolepis","Quercus_douglasii","Quercus_engelmannii","Quercus_garryana","Quercus_kelloggii","Quercus_lobata","Quercus_muehlenbergii","Quercus_wislizeni","Robinia_pseudoacacia","Salix_nigra","Salix_pyrifolia","Salix_alba","Umbellularia_californica","Olneya_tesota","Elaeagnus_angustifolia")

pnwspp[!pnwspp %in% trait_all_use$Scientific_Name]

# Create mean values for the Aesculus_sp and Fraxinus_sp
aesculus_mean <- trait_all_use %>% filter(grepl('Aesculus', Scientific_Name)) %>%
  summarize_if(is.numeric, mean, na.rm=TRUE)
fraxinus_mean <-  trait_all_use %>% filter(grepl('Fraxinus', Scientific_Name)) %>%
  summarize_if(is.numeric, mean, na.rm=TRUE)

trait_all_use <- rbind(trait_all_use, cbind(Scientific_Name=c('Aesculus_sp','Fraxinus_sp'), rbind(aesculus_mean, fraxinus_mean)))

pnwspp[!pnwspp %in% trait_all_use$Scientific_Name]

trait_all_use <- left_join(data.frame(Scientific_Name = pnwspp), trait_all_use)

########################################################
# Run phylogenetic imputation on the remaining values.

library(Rphylopars)
test_tree <- drop.tip(allfia_phylogeny, tip = allfia_phylogeny$tip.label[!allfia_phylogeny$tip.label %in% trait_all_use$Scientific_Name])

traits_goodspp <- trait_all_use[!trait_all_use$Scientific_Name %in% c('Aesculus_sp','Fraxinus_sp'),]
dimnames(traits_goodspp)[[1]] <- traits_goodspp$Scientific_Name

# Must rename first column of trait matrix to 'species' to get imputation to work.
# Also log-transform the lifespan value.

set.seed(313)
phyimp <- phylopars(trait_data = traits_goodspp %>%
                      rename(species = Scientific_Name) %>% 
                      mutate(Plant.lifespan = log(Plant.lifespan),
                             Seed.dry.mass = log(Seed.dry.mass)), 
                    tree = test_tree,
                    model = 'OU')
apply(phyimp$anc_recon, 2, min) # check for validity

# Back-transform and combine everything together.
traits_imputed <- phyimp$anc_recon[1:94,] %>%
  as.data.frame %>%
  mutate(Plant.lifespan = exp(Plant.lifespan),
         Seed.dry.mass = exp(Seed.dry.mass))

dimnames(traits_imputed)[[1]] <- dimnames(phyimp$anc_recon)[[1]][1:94]

traits_imputed <- rbind(traits_imputed, traits_imputed[grep('Aesculus|Fraxinus', dimnames(traits_imputed)[[1]]),])  
dimnames(traits_imputed)[[1]][95:96] <- c('Aesculus_sp', 'Fraxinus_sp')

traits_imputed <- cbind(Scientific_Name = dimnames(traits_imputed)[[1]], traits_imputed)
write.csv(traits_imputed, file = 'C:/Users/Q/google_drive/NASABiodiversityWG/Trait_Data/traits_imputed_22aug.csv', row.names = FALSE)
