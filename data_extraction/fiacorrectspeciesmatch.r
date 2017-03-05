# Match correct names
fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/'
load(file.path(fp, 'fiamatrices.RData'))

fiasums_plot <- fiapnw %>%
  filter(STATUSCD == 1) %>%
  mutate(ba = pi * (DIA/200)^2) %>% # basal area in m2.
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT, MEASYEAR, SPCD) %>%
  summarize(basalarea = sum(ba),
            n = n())

sppids <- sort(unique(fiasums_plot$SPCD))

fiaplotlist <- fiasums_plot %>% do(x = sapply(sppids, function(z) .$basalarea[.$SPCD == z]))
fiaplotmat <- do.call('rbind', fiaplotlist$x)

# Try a different way

area_by_sp <- function(dat, sppids) {
  areas <- numeric(length(sppids))
  for (i in 1:length(sppids)) {
    areas[i] <- sum(dat$basalarea[dat$SPCD == sppids[i]])
  }
  areas
}

fiaplotlist2 <- fiasums_plot %>% do(x = area_by_sp(., sppids))
fiaplotmat2 <- do.call('rbind', fiaplotlist2$x)

# Put correct dimnames on fiaplotmat2.
idx <- match(sppids, fiataxa$FIA.Code)
fiataxa$sciname <- paste(gsub(' ', '', fiataxa$Genus), gsub(' ', '', fiataxa$Species), sep = '_')
fiataxa$sciname[fiataxa$sciname == 'Chrysolepis_chrysophyllavar.chrysophylla'] <- 'Chrysolepis_chrysophylla'
fiataxa$sciname[fiataxa$sciname == 'Populus_balsamiferassp.Trichocarpa'] <- 'Populus_trichocarpa'

dimnames(fiaplotmat2)[[2]] <- fiataxa$sciname[idx]
