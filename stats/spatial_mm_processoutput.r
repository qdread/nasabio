# Summarize all spatial mixed models and combine all their coefficients into one data frame for plotting.
# QDR/NASABIOXGEO/26 Apr 2018

respnames_fia <- c("alpha_richness", "alpha_effspn", "alpha_phy_pa",
                   "alpha_phy", "alpha_func_pa", "alpha_func", "beta_td_sorensen_pa",
                   "beta_td_sorensen", "beta_phy_pa", "beta_phy", "beta_func_pa",
                   "beta_func", "gamma_richness", "gamma_effspn", "gamma_phy_pa",
                   "gamma_phy", "gamma_func_pa", "gamma_func")
respnames_bbs <- c("alpha_richness", "alpha_phy_pa", "alpha_func_pa",
                   "beta_td_sorensen_pa", "beta_phy_pa", "beta_func_pa", "gamma_richness",
                   "gamma_phy_pa", "gamma_func_pa")


task_table <- data.frame(taxon = rep(c('fia','bbs'), c(length(respnames_fia), length(respnames_bbs))),
                         rv = c(respnames_fia, respnames_bbs),
                         ecoregion = rep(c('HUC4','BCR','TNC'), each = length(respnames_fia) + length(respnames_bbs)),
                         stringsAsFactors = FALSE)

fp <- '/mnt/research/nasabio/temp/spammfit'
library(brms)
library(purrr)
library(reshape2)

model_coef <- list()
model_summ <- list()

n_fits <- length(dir(fp, pattern = 'fit')) # 81

# Added April 30: For the ones that didn't converge the first time, load the model fit from a different directory
fp2 <- '/mnt/research/nasabio/temp/spammfit_moreiter'

# Run summary on each fit.
for (i in 1:n_fits) {
	print(i)
  if (!(file.exists(file.path(fp2, paste0('fit', i, '.RData'))))) {
	  load(file.path(fp, paste0('fit', i, '.RData')))
  } else {
    load(file.path(fp2, paste0('fit', i, '.RData')))
  }
	
	model_summ[[i]] <- summary(fit$model, waic = FALSE, loo = FALSE, R2 = TRUE) 
	# Will need to go to parallel if we want to calculate ICs, as it takes way too long.
	model_coef[[i]] <- fit$coef
}

# Reshape coefficient data frame to slightly wider form, and then add identifying columns.
model_coef <- map2(model_coef, 1:n_fits, function(x, y) {
	coef_cast <- dcast(x, effect + region + parameter ~ stat)
	names(coef_cast) <- gsub('2.5%ile', 'q025', names(coef_cast))
	names(coef_cast) <- gsub('97.5%ile', 'q975', names(coef_cast))
	cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], coef_cast)
})

model_coef <- do.call(rbind, model_coef)

# Added 27 Apr: Reduce size of summary to just show the coefficients' convergence stats.
get_pars <- function(x, y) {
  par_df <- as.data.frame(rbind(x$fixed, x$spec_pars, x$cor_pars, x$random$region))
  cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], parameter = row.names(par_df), par_df)
}

model_r2s <- map_dbl(model_summ, 'R2')
write.csv(cbind(task_table, R2 = model_r2s), '/mnt/research/nasabio/data/fia/spatial_r2s.csv', row.names = FALSE)

model_summ <- map2(model_summ, 1:n_fits, get_pars)
model_summ <- do.call(rbind, model_summ)

model_coef_fia <- subset(model_coef, taxon == 'fia')
model_coef_bbs <- subset(model_coef, taxon == 'bbs')

model_summ_fia <-  subset(model_summ, taxon == 'fia')
model_summ_bbs <-  subset(model_summ, taxon == 'bbs')

write.csv(model_coef_fia, '/mnt/research/nasabio/data/fia/spatial_coef_fia.csv', row.names = FALSE)
write.csv(model_coef_bbs, '/mnt/research/nasabio/data/bbs/spatial_coef_bbs.csv', row.names = FALSE)

write.csv(model_summ_fia, '/mnt/research/nasabio/data/fia/spatial_summ_fia.csv', row.names = FALSE)
write.csv(model_summ_bbs, '/mnt/research/nasabio/data/bbs/spatial_summ_bbs.csv', row.names = FALSE)
