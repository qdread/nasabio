# Extract convergence diagnostics from summaries

get_pars <- function(x) {
  rbind(x$fixed, x$spec_pars, x$cor_pars, x$random$region) 
}

check_rhat <- function(x) {
  row.names(x)[x[,'Rhat'] > 1.1]
}


library(purrr)

load('/mnt/research/nasabio/data/bbs/spatial_summ_bbs.RData')
bbs_pars <- map(model_summ_bbs, get_pars)
# See if Rhats are greater than 1.1
high_rhat_bbs <- map(bbs_pars, check_rhat)

load('/mnt/research/nasabio/data/fia/spatial_summ_fia.RData')
fia_pars <- map(model_summ_fia, get_pars)
# See if Rhats are greater than 1.1
high_rhat_fia <- map(fia_pars, check_rhat)

notconverge <- which( map_int(c(high_rhat_fia, high_rhat_bbs), length) > 0 )

paste(notconverge, collapse = ',')

# Models that did not converge as of 27 Apr 2018:
# 3,13,25,38,40,43,51,55,60,61,64,73,75,78
