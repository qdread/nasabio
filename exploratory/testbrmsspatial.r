# Exploratory analysis: test conditional autoregressive models in BRMS
# Using brms 2.2.0 (older version did not have CAR implemented)

library(brms) # The likes of P. Bürkner does things!
library(rstan)
options(mc.cores = 3)
rstan_options(auto_write = TRUE)

# Run example from help documentation

## Not run: 
# generate some spatial data
east <- north <- 1:10
Grid <- expand.grid(east, north)
K <- nrow(Grid)

# set up distance and neighbourhood matrices
distance <- as.matrix(dist(Grid))
W <- array(0, c(K, K))
W[distance == 1] <- 1 	

# generate the covariates and response data
x1 <- rnorm(K)
x2 <- rnorm(K)
theta <- rnorm(K, sd = 0.05)
phi <- rmulti_normal(
  1, mu = rep(0, K), Sigma = 0.4 * exp(-0.1 * distance)
)
eta <- x1 + x2 + phi
prob <- exp(eta) / (1 + exp(eta))
size <- rep(50, K)
y <- rbinom(n = K, size = size, prob = prob)
dat <- data.frame(y, size, x1, x2)

# fit a CAR model
fit <- brm(y | trials(size) ~ x1 + x2, data = dat, 
           family = binomial(), autocor = cor_car(W)) 
summary(fit)

## End(Not run)

# generate some fake data

set.seed(2)
x1 <- sample(100)
x2 <- sample(100)
y <- 2*x1 + -1*x2 + rnorm(100, 0, 10)
region <- rep(letters[1:4], each = 25)

y[region == 'a'] <- y[region == 'a'] + 5 + 2 * x1[region == 'a']
y[region == 'b'] <- y[region == 'b'] + 25 + 5 * x1[region == 'b']

fakedat <- data.frame(y,x1,x2,region)

# generate adjacency matrix.
adjmat <- cbind(c(1,1,0,0),
                c(1,1,1,1),
                c(0,1,1,1),
                c(0,1,1,1))
dimnames(adjmat) <- list(letters[1:4], letters[1:4])

# Fit in brms

fakecar <- brm(y ~ x1 + x2 + (1|region) + (x1 - 1|region) + (x2 - 1|region), 
    data = fakedat, family = gaussian, autocor = cor_car(adjmat, formula = ~ 1|region, type = 'esicar'),
    chains = 2, iter = 2000, warmup = 1000)

summary(fakecar)
fakepars <- fakecar$fit@sim$samples[[1]] # Chain 1 of 2.

# Plot.
# Parameter names
par_names <- fakecar$fit@model_pars
stanplot(fakecar, pars = 'b', type = 'trace')

# Code to extract parameter estimates and CIs from fit object in a tidy format.
fit_sum <- summary(fakecar)
fit_sum$fixed
fit_sum$random
fit_sum$spec_pars
fit_sum$cor_pars
fre <- ranef(fakecar)
ffe <- fixef(fakecar)

# Random effects come out as a list containing a 3D array; fixed effects come out as a 2D array

# Flatten random effects into 2D array and combine with fixed effects.
ffe_df <- cbind(effect = 'fixed', as.data.frame.table(ffe, responseName = 'value'))
fre_df <- cbind(effect = 'random', as.data.frame.table(fre$region, responseName = 'value'))

rbind(ffe_df, fre_df)

library(reshape2)

ffe <- melt(ffe, varnames = c('parameter', 'stat'))
fre <- melt(fre$region, varnames = c('region', 'parameter', 'stat'))
rbind(cbind(effect = 'fixed', region = NA, ffe),
      cbind(effect = 'random', fre))

### try to fit model here.
load('C:/Users/Q/Dropbox/projects/nasabiodiv/fia_spatial_mm_dat.RData')
# Do this on cluster instead.
load('/mnt/research/nasabio/temp/fia_spatial_mm_dat.RData')

prednames <- c('elevation_5k_100_sd', 'bio1_5k_100_mean', 'geological_age_5k_100_diversity', 'soil_type_5k_100_diversity', 'bio12_5k_100_mean', 'bio12_5k_100_sd', 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean')

# This matching may not be necessary.
matchbcr <- match(fiageo$BCR, dimnames(bcr_bin)[[1]])
fiageo$BCR <- matchbcr
dimnames(bcr_bin) <- list(1:nrow(bcr_bin), 1:nrow(bcr_bin))

# Which ones have a BCR, HUC4, and TNC?
# Most of them. We can just get rid of the few that do not.

fiageo <- filter(fiageo, BCR %in% dimnames(bcr_bin)[[1]] & TNC %in% dimnames(tnc_bin)[[1]] & HUC4 %in% dimnames(huc_bin)[[1]])

test_mm1 <- fit_spatial_mm(fiageo, fiabio, prednames, resp_var = 'alpha_richness', id_var = 'PLT_CN', region_var = 'BCR', adj_matrix = bcr_bin, distribution = 'gaussian', n_chains = 2, n_iter = 500, n_warmup = 400) # This at least starts sampling, though it is taking a very long time even for the small number of iterations.

test_mm2 <- fit_spatial_mm(fiageo, fiabio, prednames, resp_var = 'beta_td_sorensen_pa', id_var = 'PLT_CN', region_var = 'BCR', adj_matrix = bcr_bin, distribution = 'beta', n_chains = 2, n_iter = 500, n_warmup = 400) 

# Do this on cluster instead.
load('/mnt/research/nasabio/temp/fia_spatial_mm_dat.RData')

# Some failed attempts below this line ------------------------------------



library(spaMM)

fakefit <- HLCor(y ~ x1 + x2 + adjacency(1|region), data = fakedat, family = gaussian(), adjMatrix = adjmat)

### lmer(DV ~ IV1 + IV2 + (1|Subject) + (IV1 - 1| Subject) + (IV2 - 1| Subject)) 

fakefit_randomslopes <- HLCor(y ~ x1 + x2 + adjacency(1|region) + adjacency(x1 - 1|region) + adjacency(x2 - 1|region), data = fakedat, family = gaussian(), adjMatrix = adjmat)
fakefit_onerandomslope <- HLCor(y ~ x1 + adjacency(x1 - 1|region), data = fakedat, family = gaussian(), adjMatrix = adjmat)

fakefit_randomslopes <- corrHLfit(y ~ x1 + x2 + adjacency(1|region) + adjacency(x1 - 1|region) + adjacency(x2 - 1|region), data = fakedat, family = gaussian(), adjMatrix = adjmat)


library(CARBayes)

adjmat <- cbind(c(0,1,0,0),
                c(1,0,1,1),
                c(0,1,0,1),
                c(0,1,1,0))
dimnames(adjmat) <- list(letters[1:4], letters[1:4])

model <- S.CARmultilevel(formula=y ~ x1 + x2, family="gaussian", ind.area=as.numeric(fakedat$region),
                         ind.re=NULL, W=adjmat, burnin=20000, n.sample=100000, data = fakedat)

# try HSAR by Dong et al.
# Use example they have

library(HSAR)
library(spdep)

data(landprice)
data("Beijingdistricts")
data(landSPDF)

# Points nested within regions (districts)!
plot(Beijingdistricts,border="green")
plot(landSPDF,add=TRUE,col="red",pch=16,cex=0.8)

model.data <- landprice[order(landprice$district.id),] # Sort increasing order of region ID.
head(model.data,50)

# the number of individuals within each neighbourhood
MM <- as.data.frame(table(model.data$district.id))
# the total number of neighbourhood, 100
Utotal <- dim(MM)[1]
Unum <- MM[,2]
Uid <- rep(c(1:Utotal),Unum)

# Which region is each one in (Z matrix)
n <- nrow(model.data)
Delta <- matrix(0,nrow=n,ncol=Utotal)
for(i in 1:Utotal) {
  Delta[Uid==i,i] <- 1
}
Delta <- as(Delta,"dgCMatrix") # Coerce to sparse.

nb.list <- poly2nb(Beijingdistricts)
mat.list <- nb2mat(nb.list,style="W") # Row normalized distance matrix
M <- as(mat.list, 'dgCMatrix') # Coerce to sparse

# extract the land parcel level spatial weights matrix
nb.25 <- dnearneigh(landSPDF,0,2500)
# to a weights matrix
dist.25 <- nbdists(nb.25,landSPDF)
dist.25 <- lapply(dist.25,function(x) exp(-0.5 * (x / 2500)^2))
mat.25 <- nb2mat(nb.25,glist=dist.25,style="W")
W <- as(mat.25,"dgCMatrix")
## run the hsar() function
res.formula <- lnprice ~ lnarea + lndcbd + dsubway + dpark + dele +
  popden + crimerate + as.factor(year)
betas= coef(lm(formula=res.formula,data=landprice)) # Coeffs from fixed effects


pars=list( rho = 0.5,lambda = 0.5, sigma2e = 2.0, sigma2u = 2.0, betas = betas )
## Not run:
res <- hsar(res.formula,data=model.data,W=W,M=M,Delta=Delta,
            burnin=5000, Nsim=10000, thinning = 1, parameters.start=pars)
summary(res)

library(classInt)
library(RColorBrewer)
x <- as.numeric(res$Mus)
breaks <- classIntervals(x,4,"fisher")$brks
groups <- cut(x,breaks,include.lowest=TRUE,labels=FALSE)
palette <- brewer.pal(4, "Blues")
plot(Beijingdistricts,col=palette[groups],border="grey")