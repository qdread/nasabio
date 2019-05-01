# Summarize spatial univariate CAR output

alpha_resp <- c('alpha_richness', 'alpha_phy_pa', 'alpha_func_pa')
beta_resp <- c('beta_td_sorensen_pa', 'beta_phy_pa', 'beta_func_pa')
gamma_resp <- c('gamma_richness', 'gamma_phy_pa', 'gamma_func_pa')

task_table <- expand.grid(taxon = c('fia', 'bbs'),
						  rv = c(alpha_resp, beta_resp, gamma_resp),
						  model = c('full','climate','space', 'geo'),
						  stringsAsFactors = FALSE)

n_fits <- nrow(task_table)

fp <- '/mnt/research/nasabio/temp/mvspam'

library(purrr)
library(reshape2)
library(tidyr)
library(dplyr)
library(foreach)
library(doParallel)

registerDoParallel(cores = 10) # Change this number as needed.

# Load all fits
model_fits <- foreach (i = 1:n_fits) %dopar% {
  load(file.path(fp, paste0('uvfit', i, '.RData')))
  list(fit = fit, fit_folds = fit_folds)
}

# Check effective sample sizes from summary to see if any need to be run again.
# There are 10000 actual samples. It looks like most of them are above 500 for effective sample size.
model_Neff <- map(model_fits, function(x) {
  pars <- x$fit$summary.results
  pars[order(pars[,'n.effective']), ]
})

map(model_Neff, function(x) x[x[,'n.effective'] < 200, ])

model_stats <- foreach (i = 1:n_fits) %dopar% {
  fit <- model_fits$fit[[i]]
  fit_folds <- model_fits$fit_folds[[i]]
  
  # Get fixed effect estimates, random effect estimates, and region-level coefficient estimates (fixed + random)
  # Manually construct the credible intervals for the estimates.
  
  # To get the coefficients, we have to add the fixed+random effects together, using sweep
  
  
  # Get fitted values that are included in the model fit object.
  
  # Join predicted with observed values to calculate the RMSE.
  
  # Manually calculate the Bayesian R-squared for the model.
  
				   # 4. Plug in dbh (x) to get posterior estimates of linear predictor of production
				  powerlaw_exp_log <- function(x, a, b, c, beta0, beta1) exp(-beta0) * x^beta1 * (-a * x ^ -b + c)
				  powerlaw_log <- function(x, beta0, beta1) exp(-beta0) * x^beta1
				  
				  # Take the log of the fitted values
				  if (prod_model == 'power') {
					prod_fitted <- log(do.call(rbind, pmap(pars, powerlaw_log, x = x)))
				  } else {
					prod_fitted <- log(do.call(rbind, pmap(pars, powerlaw_exp_log, x = x)))
				  }

				  # 5. Get residuals by subtracting log y from linear predictor
				  resids <- -1 * sweep(prod_fitted, 2, log(y))
				  
				  # 6. Calculate variances and ratio
				  pred_var <- apply(prod_fitted, 1, var)
				  resid_var <- apply(resids, 1, var)
				  r2s <- pred_var / (pred_var + resid_var)
				  
				  # Quantiles of rsq
				  r2_quant <- quantile(r2s, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
				  setNames(r2_quant, c('q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
				}
  
  # For the k-fold cross validation, extract the fitted values for only the holdout points in each fold.
  
  # Calculate the k-fold RMSE
  
  # Calculate the k-fold Bayesian R-squared
  
  # Extract any information criteria that were already calculated (for model and k-fold)
  
  # Return: fixed/random/coef, fitted values, RMSE, Bayesian R-squared, information criteria, k-fold RMSE, k-fold Bayesian R-squared, k-fold information criteria
  
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
  
    list(coef = model_coef, pred = model_pred, rmse = model_rmse, r2 = model_r2)

}

model_coef <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], model = task_table$model[y], as.data.frame(x$coef)))
model_pred <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], model = task_table$model[y], as.data.frame(x$pred)))
model_rmse <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], model = task_table$model[y], as.data.frame(x$rmse)))
model_r2 <- map_dfr(model_stats, 'r2')

write.csv(model_coef, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_coef.csv', row.names = FALSE)
write.csv(model_pred, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_pred.csv', row.names = FALSE)
write.csv(model_rmse, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_rmse.csv', row.names = FALSE)
write.csv(model_r2, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_r2.csv', row.names = FALSE)

