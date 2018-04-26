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

# Run summary on each fit.
for (i in 1:81) {
	print(i)
	load(file.path(fp, paste0('fit', i, '.RData')))
	model_summ[[i]] <- summary(fit$model)
	model_coef[[i]] <- fit$coef
}

# Reshape coefficient data frame to slightly wider form, and then add identifying columns.
model_coef <- map2(model_coef, 1:81, function(x, y) {
	coef_cast <- dcast(x, effect + region + parameter ~ stat)
	names(coef_cast) <- gsub('2.5%ile', 'q025', names(coef_cast))
	names(coef_cast) <- gsub('97.5%ile', 'q975', names(coef_cast))
	cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], coef_cast)
})

model_coef <- do.call(rbind, model_coef)

model_coef_fia <- subset(model_coef, taxon == 'fia')
model_coef_bbs <- subset(model_coef, taxon == 'bbs')

model_summ_fia <- model_summ[which(task_table$taxon == 'fia')]
model_summ_bbs <- model_summ[which(task_table$taxon == 'bbs')]

write.csv(model_coef_fia, '/mnt/research/nasabio/data/fia/spatial_coef_fia.csv', row.names = FALSE)
write.csv(model_coef_bbs, '/mnt/research/nasabio/data/bbs/spatial_coef_bbs.csv', row.names = FALSE)

save(model_summ_fia, file = '/mnt/research/nasabio/data/fia/spatial_summ_fia.RData')
save(model_summ_bbs, file = '/mnt/research/nasabio/data/bbs/spatial_summ_bbs.RData')