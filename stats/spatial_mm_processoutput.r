# Summarize all spatial mixed models and combine all their coefficients into one data frame for plotting.
# QDR/NASABIOXGEO/26 Apr 2018

# Edited 03 May: Also extract the fit statistics from the k-fold cross validation output
# Edited 01 May: Add "predict" step so we can get RMSE.

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
fpkf <- '/mnt/research/nasabio/temp/spammkfold'
library(brms)
library(purrr)
library(reshape2)

model_coef <- list()
model_summ <- list()
model_pred <- list()

n_fits <- nrow(task_table) # 81

# Run summary on each fit.
for (i in 1:n_fits) {
	print(i)
  if (!(file.exists(file.path(fp, paste0('fit', i, '.RData'))))) {
	  model_coef[[i]] <- model_summ[[i]] <- model_pred[[i]] <- 'No data'
  } else {
    load(file.path(fp, paste0('fit', i, '.RData')))
  	
	model_summ[[i]] <- summary(fit$model, waic = FALSE, loo = FALSE, R2 = TRUE) 
	model_coef[[i]] <- fit$coef
	model_pred[[i]] <- cbind(observed = fit$model$data[,1], predict(fit$model))
  }
}

# Reshape coefficient data frame to slightly wider form, and then add identifying columns.
model_coef <- map2(model_coef, 1:n_fits, function(x, y) {
	coef_cast <- dcast(x, effect + region + parameter ~ stat)
	names(coef_cast) <- gsub('2.5%ile', 'q025', names(coef_cast))
	names(coef_cast) <- gsub('97.5%ile', 'q975', names(coef_cast))
	cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], coef_cast)
})

model_coef <- do.call(rbind, model_coef)

model_r2s <- map_dbl(model_summ, 'R2')
write.csv(cbind(task_table, R2 = model_r2s), '/mnt/research/nasabio/data/fia/spatial_r2s.csv', row.names = FALSE)

# Added 27 Apr: Reduce size of summary to just show the coefficients' convergence stats.
get_pars <- function(x, y) {
  par_df <- as.data.frame(rbind(x$fixed, x$spec_pars, x$cor_pars, x$random$region))
  cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], parameter = row.names(par_df), par_df)
}


model_summ <- map2(model_summ, 1:n_fits, get_pars)

# Diagnostic step: Check R-hats of the parameters
did_not_converge <- map(model_summ, function(x) x$parameter[x$Rhat > 1.1])
which(map_int(did_not_converge, length) > 0)

model_summ <- do.call(rbind, model_summ)

model_pred <- map2_dfr(model_pred, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], as.data.frame(x)))


model_coef_fia <- subset(model_coef, taxon == 'fia')
model_coef_bbs <- subset(model_coef, taxon == 'bbs')

model_summ_fia <-  subset(model_summ, taxon == 'fia')
model_summ_bbs <-  subset(model_summ, taxon == 'bbs')

model_pred_fia <-  subset(model_pred, taxon == 'fia')
model_pred_bbs <-  subset(model_pred, taxon == 'bbs')

write.csv(model_coef_fia, '/mnt/research/nasabio/data/fia/spatial_coef_fia.csv', row.names = FALSE)
write.csv(model_coef_bbs, '/mnt/research/nasabio/data/bbs/spatial_coef_bbs.csv', row.names = FALSE)

write.csv(model_summ_fia, '/mnt/research/nasabio/data/fia/spatial_summ_fia.csv', row.names = FALSE)
write.csv(model_summ_bbs, '/mnt/research/nasabio/data/bbs/spatial_summ_bbs.csv', row.names = FALSE)

write.csv(model_pred_fia, '/mnt/research/nasabio/data/fia/spatial_pred_fia.csv', row.names = FALSE)
write.csv(model_pred_bbs, '/mnt/research/nasabio/data/bbs/spatial_pred_bbs.csv', row.names = FALSE)

# Added 03 May: Extract k-fold results.
model_kfold <- list()

# Get k-fold output for each fit
for (i in 1:n_fits) {
	print(i)
	if (!(file.exists(file.path(fpkf, paste0('kfold_', i, '.RData'))))) {
	  model_kfold[[i]] <- 'No data'
    } else {
    load(file.path(fpkf, paste0('kfold_', i, '.RData')))
  	
	model_kfold[[i]] <- kf
  }
}

model_kfold_pred <- map2_dfr(model_kfold, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], x$oos_pred))
model_kfold_stats <- map_dfr(model_kfold,function(x) {
	res <- data.frame(rmse_total=NA, rmse1=NA, rmse2=NA, rmse3=NA, rmse4=NA, rmse5=NA, kfoldic=NA, kfoldic_se=NA)
	if (length(x) > 1) res <- data.frame(rmse_total = x$rmse_total, 
									   rmse1 = x$rmse_fold[1], 
									   rmse2 = x$rmse_fold[2],
									   rmse3 = x$rmse_fold[3],
									   rmse4 = x$rmse_fold[4],
									   rmse5 = x$rmse_fold[5],
									   kfoldic = x$kfold_estimates['kfoldic','Estimate'],
									   kfoldic_se = x$kfold_estimates['kfoldic','SE'])
	return(res)								   
									   })

write.csv(cbind(task_table, model_kfold_stats), '/mnt/research/nasabio/data/fia/spatial_kfold_stats.csv', row.names = FALSE)
															   
model_kfold_pred_fia <-  subset(model_kfold_pred, taxon == 'fia')
model_kfold_pred_bbs <-  subset(model_kfold_pred, taxon == 'bbs')

write.csv(model_kfold_pred_fia, '/mnt/research/nasabio/data/fia/spatial_kfoldpred_fia.csv', row.names = FALSE)
write.csv(model_kfold_pred_bbs, '/mnt/research/nasabio/data/bbs/spatial_kfoldpred_bbs.csv', row.names = FALSE)
																			   