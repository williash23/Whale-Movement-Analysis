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

#  Set WD
jags_wd <- "C:/Users/saraw/OneDrive/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/"
setwd(jags_wd)

#   MCMC settings
ni <- 60000
nt <- 10
nb <- 1000
nc <- 3
na <- 5000


## CONTINUOUS COVARIATE ON STATE PROBABILITY

# Continuous covariate influences prob. of being in state 1 or 2
jags.dat <- list(npts = n_pts, theta = theta, l = l, 
	cov = as.numeric(scale(mod_dat$ship_whale_dist, center = TRUE)))
jags.inits <- function(){list(v = runif(2, 0.01,  5), 
	lambda = c(runif(1, 0.01, 5), NA),
	rho = runif(2, 0.01, 1), 
	mu = runif(2, -3.14159265359, 3.14159265359),
	beta0 = runif(1, 0.01, 1), 
	beta1 = runif(1, 0.01, 1))
	}
jags.parms <- c("mu", "rho", "v","lambda", "beta0", "beta1")
out_cov_dist <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "cov_on_state.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.adapt = na,
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 10)
out_cov_dist
#mcmcplots::mcmcplot(out_cov_dist$samples)

# Continuous covariate influences prob. of being in state 1 or 2
jags.dat <- list(npts = n_pts, theta = theta, l = l, 
	cov = as.numeric(scale(abs(mod_dat$ship_whale_bear), center = TRUE)))
jags.inits <- function(){list(v = runif(2, 0.01,  5), 
	lambda = c(runif(1, 0.01, 5), NA),
	rho = runif(2, 0.01, 1), 
	mu = runif(2, -3.14159265359, 3.14159265359),
	beta0 = runif(1, 0.01, 1), 
	beta1 = runif(1, 0.01, 1))
	}
jags.parms <- c("mu", "rho", "v","lambda", "beta0", "beta1")
out_cov_bear <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "cov_on_state.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.adapt = na,
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 10)
out_cov_bear
#mcmcplots::mcmcplot(out_cov_bear$samples)



## ONE STATE ONLY
jags.dat <- list(npts = n_pts, theta = theta, l = l, 
	cov = as.numeric(scale(mod_dat$ship_whale_dist, center = TRUE)))
jags.inits <- function(){list(v = runif(1, 0.01,  5), 
	lambda = runif(1, 0.01, 5),
	rho = runif(1, 0.01, 1), 
	mu = runif(1, -3.14159265359, 3.14159265359))
	}
jags.parms <- c("mu", "rho", "v","lambda", "beta0", "beta1")
out_single <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "single.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.adapt = na,
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 10)
out_single
#mcmcplots::mcmcplot(out_single$samples)



### DISCRETE RANDOM EFFECT - 2 IDX
#   Different mean for each parameter according to indicator for ship-whale context

# Ship-whale position
jags.dat <- list(npts = n_pts, theta = theta, l = l, idx = bear_ind)
jags.inits <- function(){list(v = runif(2, 0.01,  5), 
	lambda = runif(2, 0.01, 5), 
	rho = runif(2, 0.01, 1), 
	mu = runif(2, -3.14159265259, 3.14159265259))
	}
jags.parms <- c("mu", "rho", "v","lambda")
out_bear <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "intercept_per_index_2.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.adapt = na,
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 10)
out_bear
#mcmcplots::mcmcplot(out_bear$samples)


# Ship-whale distance
jags.dat <-  list(npts = n_pts, theta = theta, l = l, idx = dist_ind)
jags.inits <- function(){list(v = runif(2, 0.01,  5), 
	lambda = runif(2, 0.01, 5), 
	rho = runif(2, 0.01, 1), 
	mu = runif(2, -3.14159265359, 3.14159265359))
	}
jags.parms <- c("mu", "rho", "v","lambda")
out_dist <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "intercept_per_index_2.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.adapt = na,
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 10)
out_dist
#mcmcplots::mcmcplot(out_dist$samples)


# Whale behavior
jags.dat <- list(npts = n_pts, theta = theta, l = l, idx = beh_ind)
jags.inits <- function(){list(v = runif(2, 0.01,  5), 
	lambda = runif(2, 0.01, 5), 
	rho = runif(2, 0.01, 1), 
	mu = runif(2, -3.14159265359, 3.14159265359))
	}
jags.parms <- c("mu", "rho", "v","lambda")
out_beh <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "intercept_per_index_2.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.adapt = na,
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 10)
out_beh
#mcmcplots::mcmcplot(out_beh$samples)








S <- ggmcmc::ggs(out_cov_dist$samples) %>% dplyr::filter(Parameter != "deviance")
mu_plot <-  ggmcmc::ggs_density(S, family = "mu", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(-2.2, 2.2) +
	ylab("Density") +
	xlab(NULL)
rho_plot <-  ggmcmc::ggs_density(S, family = "rho", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0, 07.5) +
	ylab(NULL) +	
	xlab(NULL)
v_plot <-  ggmcmc::ggs_density(S, family = "v", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0.6, 1.8) +
	ylab("Density") +
	xlab("Value")
lambda_plot <-  ggmcmc::ggs_density(S, family = "lambda", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) + 
	xlim(0, 15) +
	ylab(NULL) +
	xlab("Value")
beta_plot <-  ggmcmc::ggs_density(S, family = "beta", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) + 
	xlim(-5,7) +
	ylab(NULL) +
	xlab("Value")
gridExtra::grid.arrange(mu_plot, rho_plot, v_plot, lambda_plot, beta_plot)



S <- ggmcmc::ggs(out_dist$samples) %>% dplyr::filter(Parameter != "deviance")
mu_plot <-  ggmcmc::ggs_density(S, family = "mu", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(-3.3, 3.3) +
	ylab("Density") +
	xlab(NULL)
rho_plot <-  ggmcmc::ggs_density(S, family = "rho", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0, 07.5) +
	ylab(NULL) +	
	xlab(NULL)
v_plot <-  ggmcmc::ggs_density(S, family = "v", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0.6, 1.8) +
	ylab("Density") +
	xlab("Value")
lambda_plot <-  ggmcmc::ggs_density(S, family = "lambda", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) + 
	xlim(0, 10) +
	ylab(NULL) +
	xlab("Value")
gridExtra::grid.arrange(mu_plot, rho_plot, v_plot, lambda_plot)


S <- ggmcmc::ggs(out_beh$samples) %>% dplyr::filter(Parameter != "deviance")
mu_plot <-  ggmcmc::ggs_density(S, family = "mu", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(-3.3, 3.3) +
	ylab("Density") +
	xlab(NULL)
rho_plot <-  ggmcmc::ggs_density(S, family = "rho", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0, 07.5) +
	ylab(NULL) +	
	xlab(NULL)
v_plot <-  ggmcmc::ggs_density(S, family = "v", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0.6, 1.8) +
	ylab("Density") +
	xlab("Value")
lambda_plot <-  ggmcmc::ggs_density(S, family = "lambda", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) + 
	xlim(0, 15) +
	ylab(NULL) +
	xlab("Value")
gridExtra::grid.arrange(mu_plot, rho_plot, v_plot, lambda_plot)


S <- ggmcmc::ggs(out_single$samples) %>% dplyr::filter(Parameter != "deviance")
mu_plot <-  ggmcmc::ggs_density(S, family = "mu", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(-3.3, 3.3) +
	ylab("Density") +
	xlab(NULL)
rho_plot <-  ggmcmc::ggs_density(S, family = "rho", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0, 07.5) +
	ylab(NULL) +	
	xlab(NULL)
v_plot <-  ggmcmc::ggs_density(S, family = "v", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0.6, 1.8) +
	ylab("Density") +
	xlab("Value")
lambda_plot <-  ggmcmc::ggs_density(S, family = "lambda", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) + 
	xlim(0, 10) +
	ylab(NULL) +
	xlab("Value")
gridExtra::grid.arrange(mu_plot, rho_plot, v_plot, lambda_plot)



### DISCRETE RANDOM EFFECT - 3 IDX
#   Different mean for each parameter according to indicator for ship-whale context

# Ship-whale position
jags.dat <- list(npts = n_pts, theta = theta, l = l, idx = bear_ind)
jags.inits <- function(){list(v = runif(3, 0.01,  5), 
	lambda = runif(3, 0.01, 5), 
	rho = runif(3, 0.01, 1), 
	mu = runif(3, -3.14159265359, 3.14159265359))
	}
jags.parms <- c("mu", "rho", "v","lambda")
out_bear <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "intercept_per_index_3.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.adapt = na,
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 10)
out_bear
#mcmcplots::mcmcplot(out_bear$samples)


S <- ggmcmc::ggs(out_bear$samples) %>% dplyr::filter(Parameter != "deviance")

mu_plot <-  ggmcmc::ggs_density(S, family = "mu", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(-3.3, 3.3) +
	ylab("Density") +
	xlab(NULL)
rho_plot <-  ggmcmc::ggs_density(S, family = "rho", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0, 07.5) +
	ylab(NULL) +	
	xlab(NULL)
v_plot <-  ggmcmc::ggs_density(S, family = "v", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0.6, 1.8) +
	ylab("Density") +
	xlab("Value")
lambda_plot <-  ggmcmc::ggs_density(S, family = "lambda", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) + 
	xlim(0, 10) +
	ylab(NULL) +
	xlab("Value")
	
gridExtra::grid.arrange(mu_plot, rho_plot, v_plot, lambda_plot)

	

# Ship-whale distance
jags.dat <-  list(npts = n_pts, theta = theta, l = l, idx = dist_ind)
jags.inits <- function(){list(v = runif(3, 0.01,  5), 
	lambda = runif(3, 0.01, 5), 
	rho = runif(3, 0.01, 1), 
	mu = runif(3, -3.14159265359, 3.14159265359))
	}
jags.parms <- c("mu", "rho", "v","lambda")
out_dist <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "intercept_per_index_3.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.adapt = na,
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 10)
out_dist
#mcmcplots::mcmcplot(out_dist$samples)


S <- ggmcmc::ggs(out_dist$samples) %>% dplyr::filter(Parameter != "deviance")

mu_plot <-  ggmcmc::ggs_density(S, family = "mu", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(-3.3, 3.3) +
	ylab("Density") +
	xlab(NULL)
rho_plot <-  ggmcmc::ggs_density(S, family = "rho", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0, 07.5) +
	ylab(NULL) +	
	xlab(NULL)
v_plot <-  ggmcmc::ggs_density(S, family = "v", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0.6, 1.8) +
	ylab("Density") +
	xlab("Value")
lambda_plot <-  ggmcmc::ggs_density(S, family = "lambda", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) + 
	xlim(0, 10) +
	ylab(NULL) +
	xlab("Value")
	
gridExtra::grid.arrange(mu_plot, rho_plot, v_plot, lambda_plot)



# Whale behavior
jags.dat <- list(npts = n_pts, theta = theta, l = l, idx = beh_ind)
jags.inits <- function(){list(v = runif(3, 0.01,  5), 
	lambda = runif(3, 0.01, 5), 
	rho = runif(3, 0.01, 1), 
	mu = runif(3, -3.14159265359, 3.14159265359))
	}
jags.parms <- c("mu", "rho", "v","lambda")
out_beh <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "intercept_per_index_3.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.adapt = na,
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 10)
out_beh
mcmcplots::mcmcplot(out_beh$samples)





# 4 contexts
jags.dat <- list(npts = n_pts_ind, theta = theta_ind, l = l_ind, idx = ind_sm)
jags.inits <- function(){list(v = runif(4, 0.01,  5), 
	lambda = runif(4, 0.01, 5), 
	rho = runif(4, 0.01, 1), 
	mu = runif(4, -3.14159265359, 3.14159265359))
	}
jags.parms <- c("mu", "rho", "v","lambda")
out_ind <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "C:/Users/saraw/OneDrive/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/index_four_intercepts.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.adapt = na,
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 10)
out_ind
mcmcplots::mcmcplot(out_ind$samples)






S <- ggmcmc::ggs(out_ind$samples) %>% dplyr::filter(Parameter != "deviance")

mu_plot <-  ggmcmc::ggs_density(S, family = "mu", hpd = TRUE, greek = TRUE) + 
	#theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(-3.3, 3.3) +
	ylab("Density") +
	xlab(NULL)
rho_plot <-  ggmcmc::ggs_density(S, family = "rho", hpd = TRUE, greek = TRUE) + 
	#theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0, 07.5) +
	ylab(NULL) +	
	xlab(NULL)
v_plot <-  ggmcmc::ggs_density(S, family = "v", hpd = TRUE, greek = TRUE) + 
	#theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0.6, 1.8) +
	ylab(NULL) +
	xlab(NULL)
lambda_plot <-  ggmcmc::ggs_density(S, family = "lambda", hpd = TRUE, greek = TRUE) + 
	#theme_bw() +
	guides(colour = FALSE, fill = FALSE) + 
	xlim(0, 15) +
	ylab(NULL) +
	xlab(NULL)
	
gridExtra::grid.arrange(mu_plot, rho_plot, v_plot, lambda_plot, ncol = 4)





