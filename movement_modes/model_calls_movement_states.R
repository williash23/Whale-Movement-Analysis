# Sara Williams
# 3/9/2016
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
ni <- 30000
nb <- 5000
nt <- 2
na <-5000
################################################################################

#  Run "single" model
#   Bundle data
jags.dat <- list(npts = npts_1, l = l_single, theta = theta_single)
 
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
single_dic <- dic.samples(out_single, n.iter = 5000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(single_fit)
#  Look at summary output
summary(single_fit)
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

#  Calculate Rhat
single_transit_rhat <- gelman.diag(single_fit_transit, multivariate = F)[[1]]
#  Calcualte DIC
single_transit_dic <- dic.samples(out_single_transit, n.iter = 5000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(single_fit_transit)
summary(single_fit_transit)
################################################################################

#  Run "single_station" model
#   Bundle data
jags.dat <- list(npts = npts_1_station, l = l_single_station, theta = theta_single_station)
 
#   Inits function
inits <- function(){list(v = runif(1, 0.01,  5), 
                                       lambda = runif(1, 0.01, 5), 
                                       rho = runif(1, 0.01, 1), 
                                       mu = runif(1, -3.14159265359, 3.14159265359))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho")

#  Initialize model and go through adaptation 
out_single_station <- jags.model(data = jags.dat,
                                                       file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                                                       inits = inits, 
                                                       n.chains = nc, 
                                                       n.adapt = na)

update(out_single_station, n.iter = nb)

single_fit_station <- coda.samples(out_single_station,
                                                         variable.names= params, 
                                                         n.iter = ni, 
                                                         thin = nt)

#  Calculate Rhat
single_station_rhat <- gelman.diag(single_fit_station, multivariate = F)[[1]]
#  Calcualte DIC
single_station_dic <- dic.samples(out_single_station, n.iter = 5000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(single_fit_station)
summary(single_fit_station)
################################################################################

#  Run "single_cone" model
#   Bundle data
jags.dat <- list(npts = npts_1_cone, l = l_single_cone, theta = theta_single_cone)
 
#   Inits function
inits <- function(){list(v = runif(1, 0.01,  5), 
                                       lambda = runif(1, 0.01, 5), 
                                       rho = runif(1, 0.01, 1), 
                                       mu = runif(1, -3.14159265359, 3.14159265359))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho")

out_single_cone <- jags.model(data = jags.dat,
                                                   file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                                                   inits = inits, 
                                                   n.chains = nc, 
                                                   n.adapt = na)

update(out_single_cone, n.iter = nb)

single_fit_cone <- coda.samples(out_single_cone,
                                                     variable.names= params, 
                                                     n.iter = ni, 
                                                     thin = nt)

#  Calculate Rhat
single_cone_rhat <- gelman.diag(single_fit_cone, multivariate = F)[[1]]
#  Calcualte DIC
single_cone_dic <- dic.samples(out_single_cone, n.iter = 5000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(single_fit_cone)
summary(single_fit_cone)
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
##
#  Calculate Rhat
double_rhat <- gelman.diag(double_fit, multivariate = F)[[1]]
#  Calcualte DIC
double_dic <- dic.samples(out_double, n.iter = 500, thin = 1, type = "popt")

mcmcplot(double_fit)
summary(double_fit)
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
                                        lambda = c(runif(1, 0.01, 5), NA), 
                                        eps=runif(1, 0.01, 5),
                                        rho = runif(2, 0.01, 1), 
                                        mu = runif(2, -3.14159265359, 3.14159265359),
                                        beta0 = runif(1, -5, 5),
                                        beta1 = runif(1, -5, 5))
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