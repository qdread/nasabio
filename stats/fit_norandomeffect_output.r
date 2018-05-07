# Process output of models fit with no random effect
# QDR 07 May 2018

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
                         stringsAsFactors = FALSE)

library(purrr)

# Run summary on each fit.
summarize_model <- function(fit) {

	model_summ <- summary(fit$model, waic = FALSE, loo = FALSE, R2 = TRUE) 
	model_coef <- fit$coef
	model_pred <- cbind(observed = fit$model$data[,1], predict(fit$model))

	return(list(model_summ=model_summ, model_coef=model_coef, model_pred=model_pred))
	
}

load('/mnt/research/nasabio/temp/fits_nospatial_50k.RData')
summary_50k <- map(fits, summarize_model)

load('/mnt/research/nasabio/temp/fits_nospatial_100k.RData')
summary_100k <- map(fits, summarize_model)

get_pars <- function(summ) {
	par_df <- as.data.frame(summ$model_summ$fixed)
}

model_pars_50k <- map(summary_50k, get_pars)
model_pars_100k <- map(summary_100k, get_pars)

# Check R-hats of the parameters
map_int(model_pars_50k, function(x) sum(x$Rhat > 1.1)) # All converge
map_int(model_pars_100k, function(x) sum(x$Rhat > 1.1)) # All converge

model_r2s_50k <- map_dbl(summary_50k, function(x) x$model_summ$R2)
model_r2s_100k <- map_dbl(summary_100k, function(x) x$model_summ$R2)
write.csv(cbind(task_table, radius = rep(c(50,100),each=nrow(task_table)), R2 = c(model_r2s_50k, model_r2s_100k)), '/mnt/research/nasabio/data/modelfits/nonspatial_r2s.csv', row.names = FALSE)

model_pred_50k <- map2_dfr(summary_50k, 1:nrow(task_table), function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], radius=50, as.data.frame(x$model_pred)))
model_pred_100k <- map2_dfr(summary_100k, 1:nrow(task_table), function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], radius=100, as.data.frame(x$model_pred)))


# Reshape coefficient data frame to slightly wider form, and then add identifying columns.
library(reshape2)
model_coef_50k <- map2(summary_50k, 1:nrow(task_table), function(x, y) {
	coef_cast <- dcast(x$model_coef, parameter ~ stat)
	names(coef_cast) <- gsub('2.5%ile', 'q025', names(coef_cast))
	names(coef_cast) <- gsub('97.5%ile', 'q975', names(coef_cast))
	cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], radius=50, coef_cast)
})
model_coef_100k <- map2(summary_100k, 1:nrow(task_table), function(x, y) {
	coef_cast <- dcast(x$model_coef, parameter ~ stat)
	names(coef_cast) <- gsub('2.5%ile', 'q025', names(coef_cast))
	names(coef_cast) <- gsub('97.5%ile', 'q975', names(coef_cast))
	cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], radius=100, coef_cast)
})

write.csv(do.call(rbind, model_coef_50k), '/mnt/research/nasabio/data/modelfits/nonspatial_coef_50k.csv', row.names = FALSE)
write.csv(do.call(rbind, model_coef_100k), '/mnt/research/nasabio/data/modelfits/nonspatial_coef_100k.csv', row.names = FALSE)

write.csv(model_pars_50k, '/mnt/research/nasabio/data/modelfits/nonspatial_summ_50k.csv', row.names = FALSE)
write.csv(model_pars_100k, '/mnt/research/nasabio/data/modelfits/nonspatial_summ_100k.csv', row.names = FALSE)

write.csv(model_pred_50k, '/mnt/research/nasabio/data/modelfits/nonspatial_pred_50k.csv', row.names = FALSE)
write.csv(model_pred_100k, '/mnt/research/nasabio/data/modelfits/nonspatial_pred_100k.csv', row.names = FALSE)