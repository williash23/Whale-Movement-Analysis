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
ni <- 45000
nt <- 10
nb <- 1000
nc <- 3
na <- 5000



#  Run "single" model
#   Different mean for each parameter according to indicator for ship-whale context

# Ship-whale position
jags.dat <- list(npts = n_pts_bear, theta = theta_bear, l = l_bear, idx = bear_ind_sm)
jags.inits <- function(){list(v = runif(2, 0.01,  5), 
	lambda = runif(2, 0.01, 5), 
	rho = runif(2, 0.01, 1), 
	mu = runif(2, -3.14159265359, 3.14159265359))
	}
jags.parms <- c("mu", "rho", "v","lambda")
out_bear <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "C:/Users/saraw/OneDrive/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/index_multiple_intercepts.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.adapt = na,
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 10)
out_bear
mcmcplots::mcmcplot(out_bear$samples)


# Ship-whale distance
jags.dat <- list(npts = n_pts_dist, theta = theta_dist, l = l_dist, idx = dist_ind_sm)
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
mcmcplots::mcmcplot(out_dist$samples)


# Whale behavior
jags.dat <- list(npts = n_pts_beh, theta = theta_beh, l = l_beh, idx = beh_ind)
jags.inits <- function(){list(v = runif(2, 0.01,  5), 
	lambda = runif(2, 0.01, 5), 
	rho = runif(2, 0.01, 1), 
	mu = runif(2, -3.14159265359, 3.14159265359))
	}
jags.parms <- c("mu", "rho", "v","lambda")
out_beh <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "C:/Users/saraw/OneDrive/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/index_multiple_intercepts.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.adapt = na,
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 10)
out_beh
mcmcplots::mcmcplot(out_beh$samples)




S <- ggmcmc::ggs(out_bear$samples) %>% dplyr::filter(Parameter != "deviance")

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
	ylab("Density") +
	xlab("Value")
lambda_plot <-  ggmcmc::ggs_density(S, family = "lambda", hpd = TRUE, greek = TRUE) + 
	#theme_bw() +
	guides(colour = FALSE, fill = FALSE) + 
	xlim(0, 10) +
	ylab(NULL) +
	xlab("Value")
	
gridExtra::grid.arrange(mu_plot, rho_plot, v_plot, lambda_plot)




# 3 indicies per distance and position

# Ship-whale distance
jags.dat <- list(npts = n_pts, theta = theta, l = l, idx = dist_ind)
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
mcmcplots::mcmcplot(out_dist$samples)




# Covariate
# Ship-whale distance
jags.dat <- list(npts = n_pts, l = l, 
	dist = as.numeric(scale(mod_dat$bin1000_ship_whale_dist, center = TRUE)))
jags.inits <- function(){list(v = runif(2, 0.01,  5), 
	lambda = c(runif(1, 0.01, 5), NA),
	beta0 = dnorm(0, 0.001),
	beta1 = dnorm(0, 0.001))
	}
jags.parms <- c( "v","lambda", "beta1", "beta0")
out_dist_cov<- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "cov_step.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.adapt = na,
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 10)
out_dist_cov
mcmcplots::mcmcplot(out_dist_cov$samples)





# 4 contexts
jags.dat <- list(npts = n_pts_comb, theta = theta_comb, l = l_comb, idx = comb_ind_sm)
jags.inits <- function(){list(v = runif(4, 0.01,  5), 
	lambda = runif(4, 0.01, 5), 
	rho = runif(4, 0.01, 1), 
	mu = runif(4, -3.14159265359, 3.14159265359))
	}
jags.parms <- c("mu", "rho", "v","lambda")
out_comb <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "intercept_per_index_4.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.adapt = na,
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 10)
out_comb
mcmcplots::mcmcplot(out_comb$samples)






S <- ggmcmc::ggs(out_comb$samples) %>% dplyr::filter(Parameter != "deviance")

mu_plot <-  ggmcmc::ggs_density(S, family = "mu", hpd = TRUE, greek = TRUE) + 
	#theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	#xlim(-3.3, 3.3) +
	ylab("Density") +
	xlab(NULL)
rho_plot <-  ggmcmc::ggs_density(S, family = "rho", hpd = TRUE, greek = TRUE) + 
	#theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	#xlim(-0.5, 0.75) +
	ylab(NULL) +	
	xlab(NULL)
v_plot <-  ggmcmc::ggs_density(S, family = "v", hpd = TRUE, greek = TRUE) + 
	#theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0, 2) +
	ylab(NULL) +
	xlab(NULL)
lambda_plot <-  ggmcmc::ggs_density(S, family = "lambda", hpd = TRUE, greek = TRUE) + 
	#theme_bw() +
	guides(colour = FALSE, fill = FALSE) + 
	xlim(0, 30) +
	ylab(NULL) +
	xlab(NULL)
	
gridExtra::grid.arrange(mu_plot, rho_plot, v_plot, lambda_plot, ncol = 4)





