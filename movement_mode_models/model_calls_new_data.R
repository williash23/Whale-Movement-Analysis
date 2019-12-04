# Sara Williams
# 3/9/2016, last update: 10/28/2016; 1/17/2017
# Models movement states - script to run models.
#   Formualted after code from Morales et al. 2004.
#   Run in JAGS.
################################################################################

#  Load packages
library(jagsUI)
library(mcmcplots)
library(coda)
#  Load "glm" module for JAGS
load.module("glm")
load.module("dic")
################################################################################

#   MCMC settings
ni <- 25000
nt <- 1
nb <- 10000
nc <- 3
################################################################################


#  Run "single" model - DREM
#   Discrete Random Effects Model (different mean for each parameter according to indicator for ship-whale context)
#   Bundle data
jags.dat <- list(npts = npts_1, theta = theta, l = l, reg_ind = dist_ind)
 
#   Inits function
jags.inits <- function(){list(v = runif(3, 0.01,  5), 
	lambda = runif(3, 0.01, 5), 
	rho = runif(3, 0.01, 1), 
	mu = runif(3, -3.14159265359, 3.14159265359))
	}

#   Parameters to monitor
jags.parms <- c("mu", "rho", "v","lambda")

out_DREM <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "C:/Users/saraw/OneDrive/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/single_DREM.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 4)
	
S <- ggmcmc::ggs(out_DREM$samples) %>% dplyr::filter(Parameter != "deviance")

mu_plot <-  ggmcmc::ggs_density(S, family = "mu", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
rho_plot <-  ggmcmc::ggs_density(S, family = "rho", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
v_plot <-  ggmcmc::ggs_density(S, family = "v", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
lambda_plot <-  ggmcmc::ggs_density(S, family = "lambda", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
gridExtra::grid.arrange(mu_plot, rho_plot, v_plot, lambda_plot)
#mcmcplots::mcmcplot(out_DREM$samples)


#  Initialize model and go through adaptation 
out_single <- jagsUI(data = jags.dat,
	file = "C:/Users/saraw/OneDrive/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/single_DREM.txt", 
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




#  Run "single" model
#   Bundle data
jags.dat <- list(npts = npts_1, theta = theta) #l = l, 
 
#   Inits function
inits <- function(){list(#v = runif(1, 0.01,  5), 
	#lambda = runif(1, 0.01, 5), 
	rho = runif(1, 0.01, 1), 
	mu = runif(1, -3.14159265359, 3.14159265359))
	}

#   Parameters to monitor
params <- c("mu", "rho") #"v","lambda", 

#  Initialize model and go through adaptation 
out_single <- jags.model(data = jags.dat,
	file = "C:/Users/saraw/Documents/UM_current_work/GitHub/Whale-Movement-Analysis/models_ver_9_17/single_turn_only.txt", 
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











# NOTE: due to differences in formulation of Weibull distribution between JAGS and R:
###### scale = (1/lambda)^(1/v) ######
niter <- 11000
nsamp <- 8000

rwcauchy <- function(n, mu = 0, rho = 0) {
  u = runif(n)
  V = cos(2 * pi * u)
  c = 2 * rho/(1 + rho^2)
  t <- (sign(runif(n) - 0.5) * acos((V + c)/(1 + c * V)) + mu)%%(2 * pi)
  return(t)
}

#geom_vline(aes(xintercept=mean(weight)),
#color="blue", linetype="dashed", size=1)
################################################################################



#  All data in single movement mode model
#   Select rows from posterior distribution for each parameter estimate.
#   model_fit_object[rows, columns, sheets]
#   Insert model output object as "x"

x <- single_fit
keep_1 <- sample(1:niter, nsamp, replace = F)
keep_2 <- sample(1:niter, nsamp, replace = F)
keep_3 <- sample(1:niter, nsamp, replace = F)

chain_1 <- x[[1]]
sims_1 <- chain_1[keep_1, c(1, 2)]

chain_2 <- x[[2]]
sims_2 <- chain_2[keep_2, c(1, 2)]

chain_3 <- x[[3]]
sims_3 <- chain_3[keep_3, c(1, 2)]

sims <- rbind(sims_1, sims_2, sims_3)

#steps <- numeric(length = nrow(sims))
turns <- numeric(length = nrow(sims))

for(i in 1:nrow(sims)){
	#steps[i] <- rweibull(1, sims[i,3], (1/sims[i,4])^(1/sims[i,3]))
	turns[i] <- rwcauchy(1, sims[i,1], sims[i,2])
	#rwrappedcauchy(1, sims[i,1],  sims[i,2])
}

post_sims <- as.data.frame(cbind(steps, turns))
post_sims$turns[post_sims$turns>pi]=post_sims$turns[post_sims$turns>pi]-2*pi

post_sims_plot <- post_sims %>% 
	mutate(iter_num = 1:nrow(post_sims))


turns_all <- ggplot(post_sims_plot, aes(turns)) + 
	geom_density(size = 1.25) +
	geom_vline(aes(xintercept=mean(turns)), color="black", linetype="dashed", size=1) +
	xlab("\n Turn angle (rad)") +
	ylab("Frequency \n") +
	xlim(c(-0.1, 0.1)) +
	ylim(c(0, 1000)) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
	axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2))) 
	#theme(legend.position="none")
turns_all
















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




