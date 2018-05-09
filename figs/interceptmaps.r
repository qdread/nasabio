
 # Intercept maps. 
	
arrangeMaps <- function(x, fpfig, region_name, titles, raw_names, div_type = 'incidence', rad = radius) {
  geo_name <- 'intercept' # For file name by geo variable
  x$bio_title <- titles[match(x$rv, raw_names)] # Short title by bio variable
  x <- x[match(titles, x$bio_title),] # Put bio variables in correct order
  png(file.path(fpfig, paste0(region_name, '_', rad, 'k_', div_type, '_', geo_name, '.png')), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(x$maps, x$bio_title, function(p, name) ggplotGrob(p + ggtitle(name) + tw)), nrow = 3)
  dev.off()
  return('i just made a map :-)')
}

maps_bbs_huc <- coef_all %>%
    filter(taxon == 'bbs', effect == 'random', ecoregion == 'HUC4', parameter == 'Intercept') %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(HUC4 = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, huc4, states))

maps_bbs_huc %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpbbs, region_name = 'HUC4', titles = bio_titles, raw_names = bio_names))


  maps_bbs_bcr <- coef_all %>%
    filter(taxon == 'bbs', effect == 'random', ecoregion == 'BCR', parameter == 'Intercept') %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(BCRNAME = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, bcr, states))  
  maps_bbs_bcr %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpbbs, region_name = 'BCR', titles = bio_titles, raw_names = bio_names))
  
  


  maps_bbs_tnc <- coef_all %>%
    filter(taxon == 'bbs', effect == 'random', ecoregion == 'TNC', parameter == 'Intercept') %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(ECODE_NAME = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, tnc, states))
  maps_bbs_tnc %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpbbs, region_name = 'TNC', titles = bio_titles, raw_names = bio_names))
  
  

  maps_fia_huc_incid <- coef_all %>%
    filter(taxon == 'fia', effect == 'random', ecoregion == 'HUC4', parameter == 'Intercept', rv %in% fia_bio_names_incid) %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(HUC4 = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, huc4, states))
  maps_fia_huc_incid %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfia, region_name = 'HUC4', titles = bio_titles_incidence, raw_names = fia_bio_names_incid))
  



  maps_fia_bcr_incid <- coef_all %>%
    filter(taxon == 'fia', effect == 'random', ecoregion == 'BCR', parameter == 'Intercept', rv %in% fia_bio_names_incid) %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(BCRNAME = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, bcr, states))  
  maps_fia_bcr_incid %>%
    group_by(parameter) %>%
	do(foo = arrangeMaps(., fpfig = fpfia, region_name = 'BCR', titles = bio_titles_incidence, raw_names = fia_bio_names_incid))

  maps_fia_tnc_incid <- coef_all %>%
    filter(taxon == 'fia', effect == 'random', ecoregion == 'TNC', parameter == 'Intercept', rv %in% fia_bio_names_incid) %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(ECODE_NAME = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, tnc, states))
  maps_fia_tnc_incid %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfia, region_name = 'TNC', titles = bio_titles_incidence, raw_names = fia_bio_names_incid))


  maps_fia_huc_abund <- coef_all %>%
    filter(taxon == 'fia', effect == 'random', ecoregion == 'HUC4', parameter == 'Intercept', rv %in% fia_bio_names_abund) %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(HUC4 = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, huc4, states))
  maps_fia_huc_abund %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfia, region_name = 'HUC4', titles = bio_titles_abundance, raw_names = fia_bio_names_abund, div_type = 'abundance'))


  maps_fia_bcr_abund <- coef_all %>%
    filter(taxon == 'fia', effect == 'random', ecoregion == 'BCR', parameter == 'Intercept', rv %in% fia_bio_names_abund) %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(BCRNAME = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, bcr, states))  
  maps_fia_bcr_abund %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfia, region_name = 'BCR', titles = bio_titles_abundance, raw_names = fia_bio_names_abund, div_type = 'abundance'))



  maps_fia_tnc_abund <- coef_all %>%
    filter(taxon == 'fia', effect == 'random', ecoregion == 'TNC', parameter == 'Intercept', rv %in% fia_bio_names_abund) %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(ECODE_NAME = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, tnc, states))
  maps_fia_tnc_abund %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfia, region_name = 'TNC', titles = bio_titles_abundance, raw_names = fia_bio_names_abund, div_type = 'abundance'))
