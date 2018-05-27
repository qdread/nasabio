# Summarize spatial multivariate mixed model output
# QDR/NASABIOXGEO/25 May 2018

task_table <- data.frame(taxon = rep(c('fia','bbs'), each = 3),
                         rv = c('alpha', 'beta', 'gamma'),
                         ecoregion = 'TNC',
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
  
  dimnames(model_pred)[[3]] <- resp_names
  model_pred <- melt(model_pred, varnames=c('idx','stat','response'))
  
  # Join predicted with observed values
  model_obs <- melt(cbind(idx = 1:nrow(fit$model$data), fit$model$data[, resp_names]), id.vars = 1, value.name = 'observed', variable.name = 'response')
  model_pred <- dcast(model_pred, idx + response ~ stat) %>%
    left_join(model_obs)
  
  # Here, do the RMSE for the model.
  model_rmse <- model_pred %>%
    group_by(response) %>%
    summarize(RMSE = sqrt(mean((observed-Estimate)^2)),
              range_obs = diff(range(observed)),
              rRMSE = RMSE/range_obs)
  
  list(coef = model_coef, pred = model_pred, rmse = model_rmse)

}

model_coef <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], as.data.frame(x$coef)))
model_pred <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], as.data.frame(x$pred)))
model_rmse <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], as.data.frame(x$rmse)))

write.csv(model_coef, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_coef.csv', row.names = FALSE)
write.csv(model_pred, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_pred.csv', row.names = FALSE)
write.csv(model_rmse, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_rmse.csv', row.names = FALSE)


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
  
  # Combine all the k-fold ICs and RMSEs into a single object.
  model_kfold_stats[[i]] <- map_dfr(folds_i, function(x) data.frame(                                                                        rmse_y1 = x$rmse_fold[1], 
                                                                                                                                            rmse_y2 = x$rmse_fold[2], 
                                                                                                                                            rmse_y3 = x$rmse_fold[3], 
                                                                                                                                            kfoldic = x$kfold_estimates['kfoldic','Estimate'], 
                                                                                                                                            kfoldic_se = x$kfold_estimates['kfoldic','SE']))
  
}

model_kfold_stats <- map2_dfr(model_kfold_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], fold = 1:5, as.data.frame(x)))

model_kfold_pred <- map2_dfr(model_kfold_pred, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], as.data.frame(x))) %>%
  arrange(taxon, rv, yvar, idx)

write.csv(model_kfold_stats, '/mnt/research/nasabio/data/modelfits/multivariate_kfold_rmse.csv', row.names = FALSE)
write.csv(model_kfold_pred, '/mnt/research/nasabio/data/modelfits/multivariate_kfold_pred.csv', row.names = FALSE)
