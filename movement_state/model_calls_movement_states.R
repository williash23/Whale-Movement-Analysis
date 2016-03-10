# Sara Williams
# 3/9/2016
# Models movement states - script to run models.
#   Formualted after code from Morales et al. 2004.
#   Run in JAGS.
################################################################################

#  Load packages
library(rjags)
library(mcmcplots)
#  Load "glm" module for JAGS
load.module("glm")
################################################################################

#   MCMC settings
nc <- 3
ni <- 1000
nb <- 200
nt <- 2
na <- 200
################################################################################

#  Run "single" model
#   Bundle data
jags.dat <- list(npts = npts_1, l = l_1, theta = theta_1, nind = nind_1)
 
#   Inits function
inits <- function(){list(v0 = runif(1, 0.01,  5), 
                                         lambda0 = runif(1, 0.01, 5), 
                                         rho0 = runif(1, 0.01, 1), 
                                         mu0 = runif(1, -pi, pi))
                                         }

#   Parameters to monitor
params <- c("mean.v","mean.lambda", "mean.mu", "mean.rho", "scale")

out_single <- jags.model(data = jags.dat,
                                            file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_state/models/single.txt", 
                                            inits = inits, 
                                            n.chains = nc, 
                                            n.adapt = na)

update(out_single, n.iter = nb)

single_fit <- coda.samples(out_single,
                                              variable.names= params, 
                                              n.iter = ni, 
                                              thin = nt)

mcmcplot(single_fit)
summary(single_fit)
################################################################################

#  Run "double" model
#   Bundle data
jags.dat <- list(npts = npts_1, l = l_1, theta = thet_1a, ID = ID_1, nind = nind_1)

#   Inits function
inits <- function(){list(v = runif(2, 0.01,  5), 
                                         lambda = c(NA, runif(1, 0.01, 5)), 
                                         eps=runif(1, 0.01, 5),
                                         rho = runif(2, 0.01, 1), 
                                         mu = runif(2, -pi, pi),
                                         beta0 = runif(1, -5, 5))
                                         }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "scale", "beta0", "mean.alpha", "prob1", "prob2")

out_double<- jags.model(data = jags.dat,
                                             file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_state/models/double.txt", 
                                             inits = inits, 
                                             n.chains = nc, 
                                             n.adapt = na)

update(out_double, n.iter = nb)

double_fit <- coda.samples(out_double,
                                                variable.names= params, 
                                                n.iter = ni, 
                                                thin = nt)

mcmcplot(double_fit)
summary(double_fit)
################################################################################

#  Run "double covariate" model
#   Bundle data
jags.dat <- list(npts = npts_1, l = l_1, theta = theta_1, ID = ID_1, nind = nind_1, sst = sst_1)

#   Inits function
inits <- function(){list(v = runif(2, 0.01,  5), 
                                         lambda = c(NA, runif(1, 0.01, 5)), 
                                         eps=runif(1, 0.01, 5),
                                         rho = runif(2, 0.01, 1), 
                                         mu = runif(2, -pi, pi),
                                         beta0 = runif(1, -5, 5),
                                         beta1 = runif(1, -5, 5))
                                         }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "scale", "beta0", "beta1")

out_double_cov <- jags.model(data = jags.dat,
                                                     file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_state/models/double_cov.txt", 
                                                     inits = inits, 
                                                     n.chains = nc, 
                                                     n.adapt = na)

update(out_double_cov, n.iter = nb)

double_cov_fit <- coda.samples(out_double_cov,
                                                        variable.names= params, 
                                                        n.iter = ni, 
                                                        thin = nt)

mcmcplot(double_cov_fit)
summary(double_cov_fit)
################################################################################

#  Run "double switch" model
#   Bundle data
nstate <- 2
jags.dat <- list(npts = npts, l = l, theta = theta, ID = ID, nstate = nstate, nind = nind)

#   Inits function
inits <- function(){list(v = runif(2, 0.01,  5), 
                                         lambda = c(NA, runif(1, 0.01, 5)), 
                                         eps=runif(1, 0.01, 5),
                                         rho = runif(2, 0.01, 1), 
                                         mu = runif(2, -pi, pi),
                                         phi = c(runif(1, 0.01, 1), NA),
                                         beta0 = runif(2, -5, 5))
                                         }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "scale", "beta0")

out_double_sw<- jags.model(data = jags.dat,
                                                    file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_state/models/double_switch.txt", 
                                                    inits = inits, 
                                                    n.chains = nc, 
                                                    n.adapt = na)

update(out_double_sw, n.iter = nb)

double_sw_fit <- coda.samples(out_double_sw,
                                                       variable.names= params, 
                                                       n.iter = ni, 
                                                       thin = nt)

mcmcplot(double_sw_fit)
summary(double_sw_fit)
################################################################################

#  Run "double switch covariate" model
#   Bundle data
nstate <- 2
jags.dat <- list(npts = npts, l = l, theta = theta, ID = ID, nstate = nstate, nind = nind, shipdens = ship_dens)

#   Inits function
inits <- function(){list(v = runif(2, 0.01,  5), 
                                         lambda = c(NA, runif(1, 0.01, 5)), 
                                         eps=runif(1, 0.01, 5),
                                         rho = runif(2, 0.01, 1), 
                                         mu = runif(2, -pi, pi),
                                         phi = c(runif(1, 0.01, 1), NA),
                                         beta0 = runif(2, -5, 5),
                                         beta1 = runif(2, -5, 5))
                                         }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "scale", "beta0", "beta1")

out_double_sw_cov<- jags.model(data = jags.dat,
                                                           file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_state/models/double_switch_cov.txt", 
                                                           inits = inits, 
                                                           n.chains = nc, 
                                                           n.adapt = na)

update(out_double_sw_cov, n.iter = nb)

double_sw_cov_fit <- coda.samples(out_double_sw_cov,
                                                               variable.names= params, 
                                                               n.iter = ni, 
                                                               thin = nt)

mcmcplot(double_sw_cov_fit)
summary(double_sw_cov_fit)