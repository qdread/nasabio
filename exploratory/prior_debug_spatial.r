pred_df = geodat 
resp_df = biodat 
pred_vars = prednames
resp_var = rv
id_var = siteid
region_var = ecoregion
distribution = distrib
adj_matrix = eco_mat
n_chains = NC
n_iter = NI
n_warmup = NW
require(dplyr)
  require(brms)
  require(reshape2)
  # Build formula and data
  id_df <- pred_df[, c(id_var, region_var)]
  resp_df <- resp_df[, c(id_var, resp_var)]
  pred_df <- pred_df[, c(id_var, pred_vars)]
  pred_df[,-1] <- scale(pred_df[,-1])
  pred_var_names <- names(pred_df)[-1]
  names(id_df)[2] <- 'region' # Make sure the name of the region is called region so that the random effects will specify correctly.
  region_name <- 'region'
  fixed_effects <- paste(pred_var_names, collapse = '+')
  random_effects <- paste(c(paste('(1|', region_name, ')', sep = ''), paste('(', pred_var_names, ' - 1|', region_name, ')', sep = '')), collapse = '+')
  formula_string <- paste(names(resp_df)[2], '~', fixed_effects, '+', random_effects)
  dat <- Reduce(left_join, list(id_df, resp_df, pred_df)) %>% filter(complete.cases(.))
  
mmprior <- get_prior(formula_string, data = dat, family = distribution, autocor = cor_car(adj_matrix, formula = ~ 1|region, type = 'esicar'))  

# Insert priors into the existing prior df for b and sd classes

make_stancode(formula_string, data = dat, family = distribution, autocor = cor_car(adj_matrix, formula = ~ 1|region, type = 'esicar'))  

added_priors <-  c(prior('normal(0,10)', class = 'b'),
				   prior('lognormal(1,1)', class = 'sd', group = 'region'),
				   prior('student_t(3,0,3)', class = 'sdcar')
				 )
				 
make_stancode(formula_string, data = dat, family = distribution, autocor = cor_car(adj_matrix, formula = ~ 1|region, type = 'esicar'), prior = added_priors)  