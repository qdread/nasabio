# Check convergence of models for scaling analysis
# If they converged, get coefficients, R-squared, RMSE, and information criteria.
# QDR NASAbioXgeo 18 June 2018

# Table of factors for each model
task_table <- expand.grid(taxon = c('fia','bbs'),
                          rv = c('alpha_richness', 'beta_td_sorensen_pa'),
                          random_effect = c('spatial', 'nonspatial'),
                          scale = c(5, 20, 100),
                          predictor_group = c('full', 'climate', 'geo', 'human'),
                          stringsAsFactors = FALSE)
task_table <- merge(task_table, data.frame(taxon = c('fia','bbs'), rv = c('alpha_richness','alpha_richness','beta_td_sorensen_pa','beta_td_sorensen_pa'), random_effect = rep(c('spatial','nonspatial'),each=4), predictor_group = 'null'), all = TRUE)

n_fits <- nrow(task_table)

fp <- '/mnt/research/nasabio/temp/scalingfits'

library(brms)
library(purrr)
library(reshape2)
library(tidyr)
library(dplyr)
library(foreach)
library(doParallel)

registerDoParallel(cores = 20)


# Check convergence -------------------------------------------------------

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

# Show all parameters that have Rhat of >1.1 (indicating it did not converge)
map(model_rhat, function(x) x[x[,'Rhat'] >= 1.1, ])
failed <- map_lgl(model_rhat, function(x) any(x[,'Rhat'] >= 1.1))
cbind(task_table, failed)


# Get model fit stats and predictions -------------------------------------

model_stats <- foreach (i = 1:n_fits) %dopar% {
  load(file.path(fp, paste0('fit', i, '.RData'))) # Load model
  
  model_pred <- as.data.frame(predict(fit$model)) # Returns matrix of predicted values: n rows x 4 columns (summary stats)
  names(model_pred) <- c('predicted', 'predicted_error', 'predicted_q025', 'predicted_q975')
  # Join predicted with observed values
  model_obs <- data.frame(idx = 1:nrow(fit$model$data), observed = as.numeric(fit$model$data[, 1]))
  model_pred <- cbind(task_table[i,], model_obs, model_pred)
  
  # Here, do the RMSE for the model.
  # Prediction raw values. 
  pred_raw <- predict(fit$model, summary = FALSE) # Returns n samples x n data points matrix.
  # Observed raw values
  obs_raw <- fit$model$data[, 1]
  
  # Get RMSE for each iteration and their quantiles
  rmse_quantiles <- sweep(pred_raw, 2, obs_raw, FUN = '-') %>% # Subtract predicted - observed
    melt(varnames = c('iter', 'idx')) %>%
    group_by(iter) %>%
    summarize(RMSE = sqrt(mean(value^2))) %>%
    summarize(RMSE_mean = sqrt(mean(RMSE^2)), 
              RMSE_q025 = quantile(RMSE, probs = 0.025), 
              RMSE_q975 = quantile(RMSE, probs = 0.975))
  
  # Generate range of observed data and divide this by the RMSE values to get the relative RMSE values
  rmse_quantiles <- rmse_quantiles %>%
    mutate(range_obs = diff(range(model_pred$observed))) %>%
    mutate_at(vars(starts_with('RMSE')), funs(relative = ./range_obs))
  model_rmse <- cbind(task_table[i,], rmse_quantiles)
  
  # Bayesian R-squared
  model_r2 <- bayes_R2(fit$model)
  names(model_r2) <- c('r2', 'r2_error', 'r2_q025', 'r2_q975')
  model_r2 <- cbind(task_table[i,], model_r2)
  
  # LOO information criterion
  # model_loo <- LOO(fit$model) # Takes ~5 minutes to run even if model isn't refit.
  # model_loo_est <- as.numeric(model_loo$estimates)
  # names(model_loo_est) <- c('elpd_LOO', 'p_LOO', 'LOOIC', 'elpd_LOO_se', 'p_LOO_se', 'LOOIC_se')
  # model_loo_est <- cbind(task_table[i,], t(model_loo_est))
  
  list(coef = cbind(task_table[i,], fit$coef), pred = model_pred, rmse = model_rmse, r2 = model_r2)
  
}

# After the fit stats have run for each model separately, combine all the output into data frames.
model_coef <- map_dfr(model_stats, 'coef')
model_pred <- map_dfr(model_stats, 'pred')
model_rmse <- map_dfr(model_stats, 'rmse')
model_r2 <- map_dfr(model_stats, 'r2')
model_looic <- map_dfr(model_stats, 'loo')
model_fitstats <- Reduce(left_join, list(model_rmse, model_r2, model_looic))

# Write the data frames to csv.
fp_output <- '/mnt/research/nasabio/data/modelfits'
write.csv(model_coef, file.path(fp_output, 'scaling_coef.csv'), row.names = FALSE)
write.csv(model_pred, file.path(fp_output, 'scaling_coef.csv'), row.names = FALSE)
write.csv(model_fitstats, file.path(fp_output, 'scaling_fitstats.csv'), row.names = FALSE)
