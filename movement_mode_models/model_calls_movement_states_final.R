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
ni <- 110000
nb <- 10000
nt <- 10
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
                                         file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/single.txt", 
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
                                                  file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/single.txt", 
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
                                                  file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/single.txt", 
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



#  Run "single near" model
#   Bundle data
jags.dat <- list(npts = npts_1_near, l = l_near, theta = theta_near)
 
#   Inits function
inits <- function(){list(v = runif(1, 0.01,  5), 
                                       lambda = runif(1, 0.01, 5), 
                                       rho = runif(1, 0.01, 1), 
                                       mu = runif(1, -3.14159265359, 3.14159265359))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho")

#  Initialize model and go through adaptation 
out_single_near <- jags.model(data = jags.dat,
                                                   file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/single.txt", 
                                                   inits = inits, 
                                                   n.chains = nc, 
                                                   n.adapt = na)

#  Burnin
update(out_single_near, n.iter = nb)

#  Sample posterior
single_fit_near <- coda.samples(out_single_near,
                                                     variable.names= params, 
                                                     n.iter = ni, 
                                                     thin = nt)

#  Calculate Rhat
single_near_rhat <- gelman.diag(single_fit_near, multivariate = F)[[1]]
#  Calcualte DIC
single_near_dic <- dic.samples(out_single_near, n.iter = 1000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(single_fit_near)
#  Look at summary output
summary(single_fit_near)
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
                                                file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/single.txt", 
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
                                                  file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/single.txt", 
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
                                                   file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/single.txt", 
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
                                           file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/double_loop.txt", 
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
################################################################################



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
                                                   file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/double_cov_loop.txt", 
                                                   inits = inits, 
                                                   n.chains = nc, 
                                                   n.adapt = na)

update(out_double_cov, n.iter = nb)

double_cov_fit <- coda.samples(out_double_cov,
                                                     variable.names= params, 
                                                     n.iter = ni, 
                                                     thin = nt)

#  Look at simple MCMC plots
mcmcplot(double_cov_fit)
#  Look at summary output
summary(double_cov_fit)

                                                     
#  Calculate Rhat
double_cov_rhat <- gelman.diag(double_cov_fit, multivariate = F)[[1]]
#  Calcualte DIC
double_cov_dic <- dic.samples(out_double_cov, n.iter = 500, thin = 1, type = "pD")


#################################################################################


# #  Save model objects to use later
save(single_fit, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit.RData")
save(single_fit_surf, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit_surf.RData")
save(single_fit_dive, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit_dive.RData")
save(single_fit_near, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit_near.RData")
save(single_fit_far, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit_far.RData")
save(single_fit_side, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit_side.RData")
save(single_fit_front, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit_front.RData")
save(double_fit, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/double_fit.RData")
save(double_cov_fit, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/double_cov_fit.RData")
################################################################################



