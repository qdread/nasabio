# Extract sds of coefficients to see which ones vary more spatially.
# This is the correct way of doing it, better than the previous overly complex method.

task_table <- data.frame(taxon = rep(c('fia','bbs'), each = 3),
                         rv = c('alpha', 'beta', 'gamma'),
                         ecoregion = 'TNC',
                         model = rep(c('full','climate','space','geo'),each=6),
                         stringsAsFactors = FALSE)

n_fits <- nrow(task_table)

fp <- '/mnt/research/nasabio/temp/mvspam' 

library(brms)
library(purrr)
library(dplyr)
library(foreach)
library(doParallel)

registerDoParallel(cores = 2)

# Inspect summaries.
# Only need to do the first 6, full models.
model_sds <- foreach (i = 1:6) %dopar% {
  load(file.path(fp, paste0('fit', i, '.RData')))
  sds <- summary(fit$model)$random$region[,c(1,3,4)]
  dimnames(sds)[[2]] <-  c('Estimate', 'q025', 'q975')
  sds
}

# Parse the response and parameter names out.
parse_names <- function(ns) {
  ns <- map_chr(strsplit(ns, '\\(|\\)'), 2) # get rid of parentheses around name
  response <- unlist(regmatches(ns, gregexpr('^[^_]*', ns)))
  parameter <- unlist(regmatches(ns, gregexpr('_.*$', ns)))
  parameter <- substr(parameter, 2, nchar(parameter))
  return(data.frame(response = response, parameter = parameter))
}

model_sds <- map(model_sds, function(x) cbind(parse_names(dimnames(x)[[1]]), x))

model_sds <- cbind(taxon = rep(task_table$taxon[1:6], map_int(model_sds, nrow)),
                   do.call(rbind, model_sds))

write.csv(model_sds, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_coef_variation_corrected.csv', row.names = FALSE)
