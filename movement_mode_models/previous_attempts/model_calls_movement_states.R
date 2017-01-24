# Sara Williams
# 3/9/2016, last update: 10/28/2016; 1/17/2017
# Models movement states - script to run models.
#   Formualted after code from Morales et al. 2004.
#   Run in JAGS.
################################################################################

#  Load packages
library(rjags)
library(mcmcplots)
library(coda)
#  Load "glm" module for JAGS
load.module("glm")
load.module("dic")
################################################################################

#   MCMC settings
nc <- 3
ni <- 40000
nb <- 10000
nt <- 2
na <- 5000
################################################################################



#  Run "single" model
#   Bundle data
jags.dat <- list(npts = npts_1, l = l, theta = theta)
 
#   Inits function
inits <- function(){list(v = runif(1, 0.01,  5), 
                                       lambda = runif(1, 0.01, 5), 
                                       rho = runif(1, 0.01, 1), 
                                       mu = runif(1, -3.14159265359, 3.14159265359))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho")

#  Initialize model and go through adaptation 
out_single <- jags.model(data = jags.dat,
                                          file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                                          inits = inits, 
                                          n.chains = nc, 
                                          n.adapt = na)

#  Burnin
update(out_single, n.iter = nb)

#  Sample posterior
single_fit <- coda.samples(out_single,
                                            variable.names= params, 
                                            n.iter = ni, 
                                            thin = nt)

#  Calculate Rhat
single_rhat <- gelman.diag(single_fit, multivariate = F)[[1]]
#  Calcualte DIC
single_dic <- dic.samples(out_single, n.iter = 1000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(single_fit)
#  Look at summary output
summary(single_fit)
################################################################################



#  Run "single surface interval" model
#   Bundle data
jags.dat <- list(npts = npts_1_surf, l = l_surf, theta = theta_surf)
 
#   Inits function
inits <- function(){list(v = runif(1, 0.01,  5), 
                                       lambda = runif(1, 0.01, 5), 
                                       rho = runif(1, 0.01, 1), 
                                       mu = runif(1, -3.14159265359, 3.14159265359))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho")

#  Initialize model and go through adaptation 
out_single_surf <- jags.model(data = jags.dat,
                                                     file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                                                     inits = inits, 
                                                     n.chains = nc, 
                                                     n.adapt = na)

#  Burnin
update(out_single_surf, n.iter = nb)

#  Sample posterior
single_fit_surf <- coda.samples(out_single_surf,
                                                       variable.names= params, 
                                                       n.iter = ni, 
                                                       thin = nt)

#  Calculate Rhat
single_surf_rhat <- gelman.diag(single_fit_surf, multivariate = F)[[1]]
#  Calcualte DIC
single_surf_dic <- dic.samples(out_single, n.iter = 1000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(single_fit_surf)
#  Look at summary output
summary(single_fit_surf)
################################################################################



#  Run "single deep dive" model
#   Bundle data
jags.dat <- list(npts = npts_1_dive, l = l_dive, theta = theta_dive)
 
#   Inits function
inits <- function(){list(v = runif(1, 0.01,  5), 
                                       lambda = runif(1, 0.01, 5), 
                                       rho = runif(1, 0.01, 1), 
                                       mu = runif(1, -3.14159265359, 3.14159265359))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho")

#  Initialize model and go through adaptation 
out_single_dive <- jags.model(data = jags.dat,
                                                     file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                                                     inits = inits, 
                                                     n.chains = nc, 
                                                     n.adapt = na)

#  Burnin
update(out_single_dive, n.iter = nb)

#  Sample posterior
single_fit_dive <- coda.samples(out_single_dive,
                                                       variable.names= params, 
                                                       n.iter = ni, 
                                                       thin = nt)

#  Calculate Rhat
single_dive_rhat <- gelman.diag(single_fit_dive, multivariate = F)[[1]]
#  Calcualte DIC
single_dive_dic <- dic.samples(out_single, n.iter = 1000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(single_fit_dive)
#  Look at summary output
summary(single_fit_dive)
################################################################################



#  Run "single close" model
#   Bundle data
jags.dat <- list(npts = npts_1_close, l = l_close, theta = theta_close)
 
#   Inits function
inits <- function(){list(v = runif(1, 0.01,  5), 
                                       lambda = runif(1, 0.01, 5), 
                                       rho = runif(1, 0.01, 1), 
                                       mu = runif(1, -3.14159265359, 3.14159265359))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho")

#  Initialize model and go through adaptation 
out_single_close <- jags.model(data = jags.dat,
                                          file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                                          inits = inits, 
                                          n.chains = nc, 
                                          n.adapt = na)

#  Burnin
update(out_single_close, n.iter = nb)

#  Sample posterior
single_fit_close <- coda.samples(out_single_close,
                                                    variable.names= params, 
                                                    n.iter = ni, 
                                                    thin = nt)

#  Calculate Rhat
single_close_rhat <- gelman.diag(single_fit_close, multivariate = F)[[1]]
#  Calcualte DIC
#single_close_dic <- dic.samples(out_single, n.iter = 1000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(single_fit_close)
#  Look at summary output
summary(single_fit_close)
################################################################################



#  Run "single far" model
#   Bundle data
jags.dat <- list(npts = npts_1_far, l = l_far, theta = theta_far)
 
#   Inits function
inits <- function(){list(v = runif(1, 0.01,  5), 
                                       lambda = runif(1, 0.01, 5), 
                                       rho = runif(1, 0.01, 1), 
                                       mu = runif(1, -3.14159265359, 3.14159265359))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho")

#  Initialize model and go through adaptation 
out_single_far <- jags.model(data = jags.dat,
                                          file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                                          inits = inits, 
                                          n.chains = nc, 
                                          n.adapt = na)

#  Burnin
update(out_single_far, n.iter = nb)

#  Sample posterior
single_fit_far <- coda.samples(out_single_far,
                                                     variable.names= params, 
                                                     n.iter = ni, 
                                                     thin = nt)

#  Calculate Rhat
single_far_rhat <- gelman.diag(single_fit_far, multivariate = F)[[1]]
#  Calcualte DIC
single_far_dic <- dic.samples(out_single, n.iter = 1000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(single_fit_far)
#  Look at summary output
summary(single_fit_far)
################################################################################



#  Run "single side" model
#   Bundle data
jags.dat <- list(npts = npts_1_side, l = l_side, theta = theta_side)
 
#   Inits function
inits <- function(){list(v = runif(1, 0.01,  5), 
                                       lambda = runif(1, 0.01, 5), 
                                       rho = runif(1, 0.01, 1), 
                                       mu = runif(1, -3.14159265359, 3.14159265359))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho")

#  Initialize model and go through adaptation 
out_single_side <- jags.model(data = jags.dat,
                                                      file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                                                      inits = inits, 
                                                      n.chains = nc, 
                                                      n.adapt = na)

#  Burnin
update(out_single_side, n.iter = nb)

#  Sample posterior
single_fit_side <- coda.samples(out_single_side,
                                                       variable.names= params, 
                                                       n.iter = ni, 
                                                       thin = nt)

#  Calculate Rhat
single_side_rhat <- gelman.diag(single_fit_side, multivariate = F)[[1]]
#  Calcualte DIC
single_side_dic <- dic.samples(out_single, n.iter = 1000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(single_fit_side)
#  Look at summary output
summary(single_fit_side)
################################################################################



#  Run "single front" model
#   Bundle data
jags.dat <- list(npts = npts_1_front, l = l_front, theta = theta_front)
 
#   Inits function
inits <- function(){list(v = runif(1, 0.01,  5), 
                                       lambda = runif(1, 0.01, 5), 
                                       rho = runif(1, 0.01, 1), 
                                       mu = runif(1, -3.14159265359, 3.14159265359))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho")

#  Initialize model and go through adaptation 
out_single_front <- jags.model(data = jags.dat,
                                                      file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                                                      inits = inits, 
                                                      n.chains = nc, 
                                                      n.adapt = na)

#  Burnin
update(out_single_front, n.iter = nb)

#  Sample posterior
single_fit_front <- coda.samples(out_single_front,
                                                        variable.names= params, 
                                                        n.iter = ni, 
                                                        thin = nt)

#  Calculate Rhat
single_front_rhat <- gelman.diag(single_fit_front, multivariate = F)[[1]]
#  Calcualte DIC
single_front_dic <- dic.samples(out_single, n.iter = 1000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(single_fit_front)
#  Look at summary output
summary(single_fit_front)
################################################################################



#  Run "double" model
#   Bundle data
jags.dat <- list(l = l_double, theta = theta_double, nind = nind_1, nocc = nocc_1)

#   Inits function
 inits <- function(){list(v = runif(2, 0.01, 5), 
                                        lambda = c(runif(1, 0.01, 5), NA), 
                                        eps=runif(1, 0.01, 5),
                                        rho = runif(2, 0.01, 1), 
                                        mu = runif(2, -3.14159265359, 3.14159265359),
                                        beta0 = runif(1, -5, 5))
                                        }

#   Parameters to monitor
params <- c("v","lambda", "eps", "mu", "rho", "beta0")

#  Run model in rjags.
out_double <- jags.model(data = jags.dat,
                                           file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/double_loop.txt", 
                                           inits = inits, 
                                           n.chains = nc, 
                                           n.adapt = na)
                                           
update(out_double, n.iter = nb)

double_fit <-coda.samples(out_double,
                                               variable.names= params, 
                                               n.iter = ni, 
                                               thin = nt)
                                             
#  Calculate Rhat
double_rhat <- gelman.diag(double_fit, multivariate = F)[[1]]
#  Calcualte DIC
double_dic <- dic.samples(out_double, n.iter = 1000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(double_fit)
#  Look at summary output
summary(double_fit)


# #  Save model objects to use later
# save(single_fit_surf, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/new_categories_output/single_fit_surf.RData")
# save(single_fit_dive, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/new_categories_output/single_fit_dive.RData")
# save(single_fit_close, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/new_categories_output/single_fit_close.RData")
# save(single_fit_far, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/new_categories_output/single_fit_far.RData")
# save(single_fit_side, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/new_categories_output/single_fit_side.RData")
# save(single_fit_front, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/new_categories_output/single_fit_front.RData")
################################################################################























#  Run "single_blow" model
#   Bundle data
jags.dat <- list(npts = npts_1_blow, l = l_single_blow, theta = theta_single_blow)
 
#   Inits function
inits <- function(){list(v = runif(1, 0.01,  5), 
                                       lambda = runif(1, 0.01, 5), 
                                       rho = runif(1, 0.01, 1), 
                                       mu = runif(1, -3.14159265359, 3.14159265359))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho")

#  Initialize model and go through adaptation 
out_single_blow <- jags.model(data = jags.dat,
                                                      file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                                                      inits = inits, 
                                                      n.chains = nc, 
                                                      n.adapt = na)

update(out_single_blow, n.iter = nb)

single_fit_blow <- coda.samples(out_single_blow,
                                                         variable.names= params, 
                                                         n.iter = ni, 
                                                         thin = nt)

################################################################################
#  Run "single_transit" model
#   Bundle data
jags.dat <- list(npts = npts_1_transit, l = l_single_transit, theta = theta_single_transit)
 
#   Inits function
inits <- function(){list(v = runif(1, 0.01,  5), 
                                       lambda = runif(1, 0.01, 5), 
                                       rho = runif(1, 0.01, 1), 
                                       mu = runif(1, -3.14159265359, 3.14159265359))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho")

#  Initialize model and go through adaptation 
out_single_transit <- jags.model(data = jags.dat,
                                                      file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                                                      inits = inits, 
                                                      n.chains = nc, 
                                                      n.adapt = na)

update(out_single_transit, n.iter = nb)

single_fit_transit <- coda.samples(out_single_transit,
                                                         variable.names= params, 
                                                         n.iter = ni, 
                                                         thin = nt)

# #  Calculate Rhat
# single_transit_rhat <- gelman.diag(single_fit_transit, multivariate = F)[[1]]
# #  Calcualte DIC
# single_transit_dic <- dic.samples(out_single_transit, n.iter = 5000, thin = 1, type = "pD")
# #  Look at simple MCMC plots
# mcmcplot(single_fit_transit)
# summary(single_fit_transit)
# ################################################################################

# #  Run "single_station" model
# #   Bundle data
# jags.dat <- list(npts = npts_1_station, l = l_single_station, theta = theta_single_station)
 
# #   Inits function
# inits <- function(){list(v = runif(1, 0.01,  5), 
                                       # lambda = runif(1, 0.01, 5), 
                                       # rho = runif(1, 0.01, 1), 
                                       # mu = runif(1, -3.14159265359, 3.14159265359))
                                       # }

# #   Parameters to monitor
# params <- c("v","lambda", "mu", "rho")

# #  Initialize model and go through adaptation 
# out_single_station <- jags.model(data = jags.dat,
                                                       # file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                                                       # inits = inits, 
                                                       # n.chains = nc, 
                                                       # n.adapt = na)

# update(out_single_station, n.iter = nb)

# single_fit_station <- coda.samples(out_single_station,
                                                         # variable.names= params, 
                                                         # n.iter = ni, 
                                                         # thin = nt)

# #  Calculate Rhat
# single_station_rhat <- gelman.diag(single_fit_station, multivariate = F)[[1]]
# #  Calcualte DIC
# single_station_dic <- dic.samples(out_single_station, n.iter = 5000, thin = 1, type = "pD")
# #  Look at simple MCMC plots
# mcmcplot(single_fit_station)
# summary(single_fit_station)
# ################################################################################

# #  Run "single_cone" model
# #   Bundle data
# jags.dat <- list(npts = npts_1_cone, l = l_single_cone, theta = theta_single_cone)
 
# #   Inits function
# inits <- function(){list(v = runif(1, 0.01,  5), 
                                       # lambda = runif(1, 0.01, 5), 
                                       # rho = runif(1, 0.01, 1), 
                                       # mu = runif(1, -3.14159265359, 3.14159265359))
                                       # }

# #   Parameters to monitor
# params <- c("v","lambda", "mu", "rho")

# out_single_cone <- jags.model(data = jags.dat,
                                                   # file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                                                   # inits = inits, 
                                                   # n.chains = nc, 
                                                   # n.adapt = na)

# update(out_single_cone, n.iter = nb)

# single_fit_cone <- coda.samples(out_single_cone,
                                                     # variable.names= params, 
                                                     # n.iter = ni, 
                                                     # thin = nt)

# #  Calculate Rhat
# single_cone_rhat <- gelman.diag(single_fit_cone, multivariate = F)[[1]]
# #  Calcualte DIC
# single_cone_dic <- dic.samples(out_single_cone, n.iter = 1000, thin = 1, type = "popt")
# #  Look at simple MCMC plots
# mcmcplot(single_fit_cone)
# summary(single_fit_cone)
# ################################################################################
#  Load packages
library(rjags)
library(mcmcplots)
library(coda)
#  Load "glm" module for JAGS
load.module("glm")
load.module("dic")
################################################################################

#   MCMC settings
nc <- 3
ni <- 40000
nb <- 10000
nt <- 2
na <-5000
#######################
#  Run "double" model
#   Bundle data
jags.dat <- list(l = l_double, theta = theta_double, nind = nind_1, nocc = nocc_1)

#   Inits function
 inits <- function(){list(v = runif(2, 0.01, 5), 
                                        lambda = c(NA, runif(1, 0.01, 5)), 
                                        eps=runif(1, 0.01, 5),
                                        rho = runif(2, 0.01, 1), 
                                        mu = runif(2, -3.14159265359, 3.14159265359),
                                        beta0 = runif(1, -5, 5))
                                        }

#   Parameters to monitor
params <- c("v","lambda", "eps", "mu", "rho", "beta0")

#  Run model in rjags.
out_double <- jags.model(data = jags.dat,
                                           file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/double_loop.txt", 
                                           inits = inits, 
                                           n.chains = nc, 
                                           n.adapt = na)
                                           
update(out_double, n.iter = nb)

double_fit <-coda.samples(out_double,
                                               variable.names= params, 
                                               n.iter = ni, 
                                               thin = nt)
                                             
#  Calculate Rhat
double_rhat <- gelman.diag(double_fit, multivariate = F)[[1]]
#  Calcualte DIC
double_dic <- dic.samples(out_double, n.iter = 1000, thin = 1, type = "pD")

mcmcplot(double_fit)
summary(double_fit)


#  Run model in R2jags.
R2jags_double <- jags(data = jags.dat, 
                                     inits = inits, 
                                     parameters.to.save = params, 
                                     model.file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/double_loop.txt", 
                                     n.chains=nc, 
                                     n.iter=ni, 
                                     n.burnin=nb,
                                     n.thin=nt,
                                     DIC=TRUE, 
                                     working.directory=NULL, 
                                     jags.seed = 123,
                                     progress.bar = "text", 
                                     digits=5,
                                     jags.module = c("glm","dic"))
################################################################################

# #  Save model objects to use later
# save(single_fit, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/single_fit.RData")
# save(single_fit_transit, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/single_fit_transit.RData")
# save(single_fit_station, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/single_fit_station.RData")
# save(single_fit_cone, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/single_fit_cone.RData")
# save(double_fit, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/double_fit.RData")

# ################################################################################

#  Run "double covariate" model
#   Bundle data
jags.dat <- list(l = l_double, theta = theta_double, nind = nind_1, nocc = nocc_1, ship = ship_dist_double)

#   Inits function
 inits <- function(){list(v = runif(2, 0.01, 5), 
                                        lambda = c(NA, runif(1, 0.01, 5)), 
                                        eps=runif(1, 0.01, 5),
                                        rho = runif(2, 0.01, 1), 
                                        mu = runif(2, -3.14159265359, 3.14159265359),
                                        beta0 = runif(1, -5, 5))
                                        }

#   Parameters to monitor
params <- c("v","lambda", "eps", "mu", "rho", "scale", "beta0", "beta1")

out_double_cov <- jags.model(data = jags.dat,
                                                   file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/double_cov_loop.txt", 
                                                   inits = inits, 
                                                   n.chains = nc, 
                                                   n.adapt = na)

update(out_double_cov, n.iter = nb)

double_cov_fit <- coda.samples(out_double_cov,
                                                     variable.names= params, 
                                                     n.iter = ni, 
                                                     thin = nt)
#  Calculate Rhat
double_cov_rhat <- gelman.diag(double_cov_fit, multivariate = F)[[1]]
#  Calcualte DIC
double_cov_dic <- dic.samples(out_double_cov, n.iter = 500, thin = 1, type = "pD")

mcmcplot(double_cov_fit)
summary(double_cov_fit)
#################################################################################

#  Run "double switch" model
#   Bundle data
nstate <- 2
jags.dat <- list(l = l, theta = theta, nstate = nstate, nind = nind, nocc = nocc)

#   Inits function
inits <- function(){list(v = runif(2, 0.01, 5), 
                                        lambda = c(runif(1, 0.01, 5), NA), 
                                        eps=runif(1, 0.01, 5),
                                        rho = runif(2, 0.01, 1), 
                                        mu = runif(2, -3.14159265359, 3.14159265359),
                                       phi = c(runif(1, 0.01, 1), NA),
                                       beta0 = runif(2, -5, 5))
                                      }

# #   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "scale", "beta0", "eps", "tau_alpha")

# #  Initialize model and go through adaptation 
out_double_sw <- jags.model(data = jags.dat,
                                                  file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/double_switch_loop.txt", 
                                                  inits = inits, 
                                                  n.chains = nc, 
                                                  n.adapt = na)

#  Burnin
update(out_double_sw, n.iter = nb)

#  Sample posterior
double_sw_fit <- coda.samples(out_double_sw,
                                                    variable.names= params, 
                                                    n.iter = ni, 
                                                    thin = nt)

# #  Calculate Rhat
# double_sw_rhat <- gelman.diag(double_sw_fit, multivariate = F)[[1]]
# #  Calcualte DIC
# double_sw_dic <- dic.samples(out_double_sw, n.iter = 5000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(double_sw_fit)
#  Look at summary output
summary(double_sw_fit)
################################################################################

#  Run "double switch covariate" model
#   Bundle data
nstate <- 2
jags.dat <- list(l = l, theta = theta, nstate = nstate, nind = nind, nocc = nocc, shipdist  = ship_dist)

#   Inits function
inits <- function(){list(v = runif(2, 0.01, 5), 
                                        lambda = c(runif(1, 0.01, 5), NA), 
                                        eps=runif(1, 0.01, 5),
                                        rho = runif(2, 0.01, 1), 
                                        mu = runif(2, -3.14159265359, 3.14159265359),
                                      phi = c(runif(1, 0.01, 1), NA),
                                      beta0 = runif(2, -5, 5),
                                      beta1 = runif(2, -5, 5))
                                      }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "scale", "beta0", "beta1", "eps")

#  Initialize model and go through adaptation 
out_double_sw_cov <- jags.model(data = jags.dat,
                                                         file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/double_switch_cov_loop.txt", 
                                                         inits = inits, 
                                                         n.chains = nc, 
                                                         n.adapt = na)

#  Burnin
update(out_double_sw_cov, n.iter = nb)

#  Sample posterior
double_sw_cov_fit <- coda.samples(out_double_sw_cov,
                                                           variable.names= params, 
                                                           n.iter = ni, 
                                                           thin = nt)
# #  Calculate Rhat
# double_sw_rhat <- gelman.diag(double_sw_fit, multivariate = F)[[1]]
# #  Calcualte DIC
# double_sw_dic <- dic.samples(out_double_sw, n.iter = 5000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(double_sw_cov_fit)
#  Look at summary output
summary(double_sw_cov_fit)
################################################################################

# #  Run "single step only" model
# #   Bundle data
# jags.dat <- list(npts = npts_step, l = l_single_step)
 
# #   Inits function
# inits <- function(){list(v = runif(1, 0.01,  5), 
                                       # lambda = runif(1, 0.01, 5))
                                       # }

# #   Parameters to monitor
# params <- c("v","lambda")

# #  Initialize model and go through adaptation 
# out_single_step <- jags.model(data = jags.dat,
                                                   # file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single_step.txt", 
                                                   # inits = inits, 
                                                   # n.chains = nc, 
                                                   # n.adapt = na)

# #  Burnin
# update(out_single_step, n.iter = nb)

# #  Sample posterior
# single_step_fit <- coda.samples(out_single_step,
                                                     # variable.names= params, 
                                                     # n.iter = ni, 
                                                     # thin = nt)

# #  Calculate Rhat
# single_step_rhat <- gelman.diag(single_step_fit, multivariate = F)[[1]]
# #  Calcualte DIC
# single_step_dic <- dic.samples(out_single_step, n.iter = 5000, thin = 1, type = "pD")
# #  Look at simple MCMC plots
# mcmcplot(single_step_fit)
# #  Look at summary output
# summary(single_step_fit)
# ################################################################################

# #  Run "double step only" model
# #   Bundle data
# jags.dat <- list(l = l_double_step, nind = nind_step, nocc = nocc_step)

# #   Inits function
 # inits <- function(){list(v = runif(2, 0.01, 5), 
                                        # lambda = c(runif(1, 0.01, 5),NA), 
                                        # eps=runif(1, 0.01, 5),
                                        # beta0 = runif(1, -5, 5))
                                        # }

# #   Parameters to monitor
# params <- c("v","lambda", "beta0", "eps")

# out_double_step <- jags.model(data = jags.dat,
                                                    # file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/double_loop_step.txt", 
                                                    # inits = inits, 
                                                    # n.chains = nc, 
                                                    # n.adapt = na)

# update(out_double_step, n.iter = nb)

# double_step_fit <- coda.samples(out_double_step,
                                                       # variable.names= params, 
                                                       # n.iter = ni, 
                                                       # thin = nt)

# #  Calculate Rhat
# double_step_rhat <- gelman.diag(double_step_fit, multivariate = F)[[1]]
# #  Calcualte DIC
# double_step_dic <- dic.samples(out_double_step, n.iter = 5000, thin = 1, type = "pD")

# mcmcplot(double_step_fit)
# summary(double_step_fit)
# ################################################################################