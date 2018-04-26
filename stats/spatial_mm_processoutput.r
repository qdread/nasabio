# Summarize all spatial mixed models and combine all their coefficients into one data frame for plotting.
# QDR/NASABIOXGEO/26 Apr 2018

task_table <- data.frame(taxon = rep(c('fia','bbs'), c(length(respnames_fia), length(respnames_bbs))),
                         rv = c(respnames_fia, respnames_bbs),
                         ecoregion = rep(c('HUC4','BCR','TNC'), each = length(respnames_fia) + length(respnames_bbs)),
                         stringsAsFactors = FALSE)

fp <- '/mnt/research/nasabio/temp/spammfit'
library(brms)
library(purrr)
library(reshape2)

# Run summary on each fit.
model_summ <- map(1:81, function(i) {
	load(file.path(fp, paste0('fit', i, '.RData'))
	summ <- summary(fit$model)
	list(summ = summ, coef = fit$coef)
}

model_coef <- map(model_summ, 'coef')
model_summ <- map(model_summ, 'summ')

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

model_summ_fia <- model_summ[which(task_table$taxon == 'fia']
model_summ_bbs <- model_summ[which(task_table$taxon == 'bbs']

write.csv(model_coef_fia, '/mnt/research/nasabio/data/fia/spatial_coef_fia.csv', row.names = FALSE)
write.csv(model_coef_bbs, '/mnt/research/nasabio/data/bbs/spatial_coef_bbs.csv', row.names = FALSE)

save(model_summ_fia, file = '/mnt/research/nasabio/data/fia/spatial_summ_fia.RData')
save(model_summ_bbs, file = '/mnt/research/nasabio/data/bbs/spatial_summ_bbs.RData')