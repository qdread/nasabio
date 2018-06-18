# Summarize spatial multivariate mixed model output
# QDR/NASABIOXGEO/25 May 2018

# Modified 14 June: include newer "null" and subset models
# Modified 4 June: get RMSE for each iteration so we can put a credible interval on it.

task_table <- data.frame(taxon = rep(c('fia','bbs'), each = 3),
                         rv = c('alpha', 'beta', 'gamma'),
                         ecoregion = 'TNC',
                         model = rep(c('full','climate','space'),each=6),
                         stringsAsFactors = FALSE)

n_fits <- nrow(task_table)

fp <- '/mnt/research/nasabio/temp/mvspam'

library(brms)
library(purrr)
library(reshape2)
library(tidyr)
library(dplyr)
library(foreach)
library(doParallel)

registerDoParallel(cores = n_fits)

# Inspect summaries.
model_summ <- foreach (i = 1:n_fits) %dopar% {
  load(file.path(fp, paste0('fit', i, '.RData')))
  summary(fit$model)
}

# Check R-hat from summary.
model_rhat <- map(model_summ, function(x) {
  pars <- with(x, rbind(fixed, spec_pars, cor_pars, random[[1]]))
  pars[order(pars[,'Rhat'], decreasing = TRUE), ]
})

map(model_rhat, function(x) x[x[,'Rhat'] >= 1.1, ])

model_stats <- foreach (i = 1:n_fits) %dopar% {
  load(file.path(fp, paste0('fit', i, '.RData')))
  
  # Get the correct variable names and apply them where needed.
  raw_resp_names <- fit$model$formula$responses
  resp_idx <- match(raw_resp_names, gsub('_', '', names(fit$model$data)))
  resp_names <- names(fit$model$data)[resp_idx]
  
  fit$coef <- fit$coef %>%
    separate(parameter, into = c('response', 'parameter'), extra = 'merge') %>%
    mutate(response = resp_names[match(response, raw_resp_names)])
  
  model_coef <- fit$coef # Includes fixed, random, and coefficient.
  model_pred <- predict(fit$model) # Returns raw array: n data x 4 stats x n response variables.
  # The predict() call takes a long time (~20 min or so in some cases)
  
  dimnames(model_pred)[[3]] <- resp_names
  model_pred <- melt(model_pred, varnames=c('idx','stat','response'))
  
  # Join predicted with observed values
  model_obs <- melt(cbind(idx = 1:nrow(fit$model$data), fit$model$data[, resp_names]), id.vars = 1, value.name = 'observed', variable.name = 'response')
  model_pred <- dcast(model_pred, idx + response ~ stat) %>%
    left_join(model_obs)
  
  # Here, do the RMSE for the model.
  # Prediction raw values. 
  pred_raw <- predict(fit$model, summary = FALSE)
  dimnames(pred_raw)[[3]] <- resp_names
  
  # Observed raw values
  obs_raw <- fit$model$data[, resp_names]
  
  # Get RMSE for each iteration and their quantiles
  rmse_quantiles <- sweep(pred_raw, 2:3, as.matrix(obs_raw), FUN = '-') %>% # Subtract predicted - observed
    melt(varnames = c('iter', 'idx', 'response')) %>%
    group_by(response, iter) %>%
    summarize(RMSE = sqrt(mean(value^2))) %>%
    ungroup %>% group_by(response) %>%
    summarize(RMSE_mean = sqrt(mean(RMSE^2)), 
              RMSE_q025 = quantile(RMSE, probs = 0.025), 
              RMSE_q975 = quantile(RMSE, probs = 0.975))
  
  # Generate ranges of observed data and divide this by the RMSE values to get the relative RMSE values
  model_rmse <- model_pred %>%
    group_by(response) %>%
    summarize(range_obs = diff(range(observed))) %>%
    left_join(rmse_quantiles) %>%
    mutate_at(vars(starts_with('RMSE')), funs(relative = ./range_obs))
  
  # Bayesian R-squared
  model_r2 <- cbind(task_table[i, ], response = resp_names, bayes_R2(fit$model))
  
  # Addition 07 June: calculate coefficient estimates for all iterations, then get the SD of them (with credible interval)
  # This will be used to say which relationships vary more spatially.
  co_raw <- coef(fit$model, summary = FALSE)
  co_sds <- apply(co_raw[[1]], c(1,3), sd)
  co_mean <- apply(co_raw[[1]], c(1,3), mean)
  co_cv <- co_sds/abs(co_mean)
  
  co_cv_quant <- as.data.frame(t(apply(co_cv, 2, quantile, probs = c(0.5, 0.025, 0.975))))
  names(co_cv_quant) <- c('cv', 'cv_q025', 'cv_q975')
  
  co_cv_quant <- cbind(name = row.names(co_cv_quant), co_cv_quant)
  
  co_cv_quant <- co_cv_quant %>%
    separate(name, into = c('response', 'parameter'), extra = 'merge') %>%
    mutate(response = resp_names[match(response, raw_resp_names)])
  
  
  list(coef = model_coef, pred = model_pred, rmse = model_rmse, r2 = model_r2, coef_variation = co_cv_quant)

}

model_coef <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], model = task_table$model[y], as.data.frame(x$coef)))
model_pred <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], model = task_table$model[y], as.data.frame(x$pred)))
model_rmse <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], model = task_table$model[y], as.data.frame(x$rmse)))
model_r2 <- map_dfr(model_stats, 'r2')
model_coef_variation <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], model = task_table$model[y], x$coef_variation))


write.csv(model_coef, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_coef.csv', row.names = FALSE)
write.csv(model_pred, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_pred.csv', row.names = FALSE)
write.csv(model_rmse, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_rmse.csv', row.names = FALSE)
write.csv(model_r2, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_r2.csv', row.names = FALSE)
write.csv(model_coef_variation, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_coef_variation.csv', row.names = FALSE)

# K-fold output -----------------------------------------------------------

model_kfold_pred <- list()
model_kfold_stats <- list()

# For each fit, load each k-fold subset (5) and combine the outputs.

for (i in 1:n_fits) {
  folds_i <- map(1:5, function(j) {
    load(file.path(fp, paste0('kfold_', i, '_', j, '.RData')))
    kf
  })
  # Combine all the predicted values into a single data frame.
  model_kfold_pred[[i]] <- map_dfr(folds_i, 'oos_pred')
  
  # Calculate RMSE for each fold and each iteration
  RMSE_all <- model_kfold_pred[[i]] %>%
    group_by(fold, iter, response) %>%
    summarize(RMSE = sqrt(mean(value^2)))
  RMSE_quantiles <- RMSE_all %>%
    ungroup %>%
    group_by(response) %>%
    summarize(kfold_RMSE_mean = sqrt(mean(RMSE^2)), 
              kfold_RMSE_q025 = quantile(RMSE, probs = 0.025), 
              kfold_RMSE_q975 = quantile(RMSE, probs = 0.975))
  RMSE_quantiles_byfold <- RMSE_all %>%
    ungroup %>%
    group_by(fold, response) %>%
    summarize(kfold_RMSE_mean = sqrt(mean(RMSE^2)), 
              kfold_RMSE_q025 = quantile(RMSE, probs = 0.025), 
              kfold_RMSE_q975 = quantile(RMSE, probs = 0.975))
  
  
  # Combine all the k-fold ICs and RMSEs into a single object.
  model_kfold_stats[[i]] <- bind_rows(data.frame(fold = NA, RMSE_quantiles),
                                  data.frame(RMSE_quantiles_byfold, 
                                  map_dfr(folds_i, function(x) data.frame(kfoldic = rep(x$kfold_estimates['kfoldic','Estimate'],3),
                                                                          kfoldic_se = rep(x$kfold_estimates['kfoldic','SE'],3)))))
    
}

model_kfold_stats <- map2_dfr(model_kfold_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], model = task_table$model[y], as.data.frame(x)))

# model_kfold_pred <- map2_dfr(model_kfold_pred, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], as.data.frame(x))) %>%
#   arrange(taxon, rv, yvar, idx)

write.csv(model_kfold_stats, '/mnt/research/nasabio/data/modelfits/multivariate_kfold_rmse.csv', row.names = FALSE)
# write.csv(model_kfold_pred, '/mnt/research/nasabio/data/modelfits/multivariate_kfold_pred.csv', row.names = FALSE)
