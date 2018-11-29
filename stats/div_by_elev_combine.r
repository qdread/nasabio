# Code to compile the div_by_elev subsamples for the revision.
# 1e5 iterations were run on the cluster. Load them there, calculate summary information, and save to make plots locally.

radii <- c(5, 10, 20, 50, 100)
div_names <- c('alpha_diversity','beta_diversity','gamma_diversity')

# x-values for predicted y-values
xrange <- c(2, 1211)
newx <- round(seq(xrange[1], xrange[2], length.out = 50))

fp <- '/mnt/research/nasabio/temp/pnwfit'

library(dplyr)

pred_list <- list()
r2_list <- list()
slope_list <- list()

fitnames <- dir(fp, pattern='fit_')

for (i in 1:100) {
  load(file.path(fp, fitnames[i]))
  pred_list[[i]] <- pred_val_lm_array
  r2_list[[i]] <- r2_lm_array
  slope_list[[i]] <- coef_array
  print(i)
}

library(abind)
pred_val_lm_array <- do.call('abind', c(pred_list, along = 3)) # very big.
r2_lm_array <- do.call('abind', c(r2_list, along = 3))
coef_array <- do.call('abind', c(slope_list, along=3))

# Convert array to data frame
library(reshape2)

dimnames(r2_lm_array) <- list(div_names, radii, NULL)
dimnames(pred_val_lm_array) <- list(div_names, radii, NULL, NULL)
dimnames(coef_array) <- list(div_names, radii, NULL)
coef_df <- melt(coef_array, varnames = c('diversity_type', 'radius', 'iteration'))
r2_lm_df <- melt(r2_lm_array, varnames = c('diversity_type', 'radius', 'iteration'))
pred_val_lm_df <- melt(pred_val_lm_array, varnames = c('diversity_type', 'radius', 'iteration', 'x'))
pred_val_lm_df$x <- newx[pred_val_lm_df$x]

# Predicted values
pred_val_lm_quant <- pred_val_lm_df %>%
  group_by(diversity_type, radius, x) %>%
  summarize(pred_y = quantile(value, probs = 0.5, na.rm = TRUE),
            pred_y_q025 = quantile(value, probs = 0.025, na.rm = TRUE),
            pred_y_q975 = quantile(value, probs = 0.975, na.rm = TRUE),
            pred_y_q25 = quantile(value, probs = 0.25, na.rm = TRUE),
            pred_y_q75 = quantile(value, probs = 0.75, na.rm = TRUE),
            pred_y_mean = mean(value, na.rm = TRUE)) %>%
  arrange(diversity_type, radius, x)

# R2s
r2_lm_quant <-r2_lm_df %>%
  group_by(diversity_type, radius) %>%
  summarize(r2 = quantile(value, probs = 0.5, na.rm = TRUE),
            r2_q025 = quantile(value, probs = 0.025, na.rm = TRUE),
            r2_q975 = quantile(value, probs = 0.975, na.rm = TRUE),
            r2_q25 = quantile(value, probs = 0.25, na.rm = TRUE),
            r2_q75 = quantile(value, probs = 0.75, na.rm = TRUE),
            r2_mean = mean(value, na.rm = TRUE)) %>%
  arrange(diversity_type, radius)

# Slopes or coefficients
coef_quant <- coef_df %>%
  group_by(diversity_type, radius) %>%
  summarize(coef = quantile(value, probs = 0.5, na.rm = TRUE),
            coef_q025 = quantile(value, probs = 0.025, na.rm = TRUE),
            coef_q975 = quantile(value, probs = 0.975, na.rm = TRUE),
            coef_q25 = quantile(value, probs = 0.25, na.rm = TRUE),
            coef_q75 = quantile(value, probs = 0.75, na.rm = TRUE),
            coef_mean = mean(value, na.rm = TRUE)) %>%
  arrange(diversity_type, radius)


save(pred_val_lm_quant, r2_lm_quant, coef_quant, file = file.path(fp, 'fiafitplotdat_pnw.R'))
