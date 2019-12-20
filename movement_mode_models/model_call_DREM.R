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
nt <- 2
nb <- 15000
nc <- 3



#  Run "single" model - DREM
#   Discrete Random Effects Model (different mean for each parameter according to indicator for 
#   ship-whale context)
#   Bundle data
jags.dat <- list(npts = npts150, theta = theta150, l = l150, idx = dist_ind150)
 
#   Inits function
jags.inits <- function(){list(v = runif(3, 0.01,  5), 
	lambda = runif(3, 0.01, 5), 
	rho = runif(3, 0.01, 1), 
	mu = runif(3, -3.14159265359, 3.14159265359))
	}

#   Parameters to monitor
jags.parms <- c("mu", "rho", "v","lambda")

out <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "C:/Users/saraw/OneDrive/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/index_multiple_intercepts.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 6)

out_rand_int_dist <- out
#out_rand_int_bear <- out

S <- ggmcmc::ggs(out_rand_int_bear$samples) %>% dplyr::filter(Parameter != "deviance")

mu_plot <-  ggmcmc::ggs_density(S, family = "mu", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	#xlim(-3, 3) +
	ylab("Density") +
	xlab(NULL)
rho_plot <-  ggmcmc::ggs_density(S, family = "rho", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	#xlim(0, 0.5) +
	ylab(NULL) +	
	xlab(NULL)
v_plot <-  ggmcmc::ggs_density(S, family = "v", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	#xlim(0.6, 1.2) +
	ylab("Density") +
	xlab("Value")
lambda_plot <-  ggmcmc::ggs_density(S, family = "lambda", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) + 
	#xlim(1.5, 4) +
	ylab(NULL) +
	xlab("Value")
	
gridExtra::grid.arrange(mu_plot, rho_plot, v_plot, lambda_plot)

ggmcmc::ci(S)
 
 

S <- ggmcmc::ggs(out_rand_int_dist$samples) %>% dplyr::filter(Parameter != "deviance")

mu_plot <-  ggmcmc::ggs_density(S, family = "mu", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(-3, 3) +
	ylab("Density") +
	xlab(NULL)
rho_plot <-  ggmcmc::ggs_density(S, family = "rho", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0, 0.6) +
	ylab(NULL) +	
	xlab(NULL)
v_plot <-  ggmcmc::ggs_density(S, family = "v", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0.8, 1.8) +
	ylab("Density") +
	xlab("Value")
lambda_plot <-  ggmcmc::ggs_density(S, family = "lambda", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) + 
	xlim(1, 10) +
	ylab(NULL) +
	xlab("Value")
gridExtra::grid.arrange(mu_plot, rho_plot, v_plot, lambda_plot)

ggmcmc::ci(S)












#  Run "double" model - DREM of bearing on movement parameters
#   Discrete Random Effects Model (random effect on each parameter according to indicator for ship-whale 
#    context)
#   Bundle data
# jags.dat <- list(npts = npts, nind = nind, nocc = nocc,
	# theta = theta_double, l = l_double, 
	# bear_ind = bear_ind)
 
 #   Bundle data
jags.dat <- list(npts = npts5, nind = nind5, nocc = nocc5,
	theta = theta_double5, l = l_double5, 
	bear_ind = bear_ind5)
 
 # #   Bundle data
# jags.dat <- list(npts = npts10, nind = nind10, nocc = nocc10,
	# theta = theta_double10, l = l_double10, 
	# bear_ind = bear_ind10)
 
 # #   Bundle data
# jags.dat <- list(npts = npts20, nind = nind20, nocc = nocc20,
	# theta = theta_double20, l = l_double20, 
	# bear_ind = bear_ind20)
 
 # #   Bundle data
# jags.dat <- list(npts = npts60, nind = nind60, nocc = nocc60,
	# theta = theta_double60, l = l_double60, 
	# bear_ind = bear_ind60)
 
  #   Bundle data
jags.dat <- list(npts = npts120, nind = nind120, nocc = nocc120,
	theta = theta_double120, l = l_double120, 
	bear_ind = bear_ind120)
	
 
#   Inits function
jags.inits <- function(){list(alpha_v = runif(1, 0.01,  2), 
	alpha_lambda = runif(1, 0.01, 30), 
	alpha_rho = runif(1, 0.01, 1), 
	alpha_mu = runif(1, -3.14159265359, 3.14159265359))
	}

#   Parameters to monitor
jags.parms <- c("alpha_mu", "alpha_rho", "alpha_v","alpha_lambda", 
	"dre_bear_mu", "dre_bear_rho", "dre_bear_v", "dre_bear_lambda")

out <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "./DREM_bear_cov.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 6)
	
out_DREM_bear120 <- out
save(out_DREM_bear120, file = "out_DREM_bear120.Rdata")



S <- ggmcmc::ggs(out$samples) %>% dplyr::filter(Parameter != "deviance")

mu_plot <-  ggmcmc::ggs_density(S, family = "alpha_mu", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
mu_re_plot <-  ggmcmc::ggs_density(S, family = "dre_bear_mu", hpd = TRUE, greek = FALSE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
gridExtra::grid.arrange(mu_plot, mu_re_plot, ncol = 2)

rho_plot <-  ggmcmc::ggs_density(S, family = "alpha_rho", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
rho_re_plot <-  ggmcmc::ggs_density(S, family = "dre_bear_rho", hpd = TRUE, greek = FALSE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
gridExtra::grid.arrange(rho_plot, rho_re_plot, ncol = 2)

v_plot <-  ggmcmc::ggs_density(S, family = "alpha_v", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
v_re_plot <-  ggmcmc::ggs_density(S, family = "dre_bear_v", hpd = TRUE, greek = FALSE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
gridExtra::grid.arrange(v_plot, v_re_plot, ncol = 2)
	
lambda_plot <-  ggmcmc::ggs_density(S, family = "alpha_lambda", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
lambda_re_plot <-  ggmcmc::ggs_density(S, family = "dre_bear_lambda", hpd = TRUE, greek = FALSE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
gridExtra::grid.arrange(lambda_plot, lambda_re_plot, ncol = 2)
	
	



#  Run "double" model - DREM of distance on movement parameters
#   Discrete Random Effects Model (random effect on each parameter according to indicator for ship-whale 
#    context)
#   Bundle data
jags.dat <- list(npts = npts, nind = nind, nocc = nocc,
	theta = theta_double, l = l_double, 
	dist_ind = dist_ind)
 
  #   Bundle data
jags.dat <- list(npts = npts5, nind = nind5, nocc = nocc5,
	theta = theta_double5, l = l_double5, 
	dist_ind = dist_ind5)
 
 #   Bundle data
jags.dat <- list(npts = npts10, nind = nind10, nocc = nocc10,
	theta = theta_double10, l = l_double10, 
	dist_ind = dist_ind10)
 
 #   Bundle data
jags.dat <- list(npts = npts20, nind = nind20, nocc = nocc20,
	theta = theta_double20, l = l_double20, 
	dist_ind = dist_ind20)
 
 #   Bundle data
jags.dat <- list(npts = npts60, nind = nind60, nocc = nocc60,
	theta = theta_double60, l = l_double60, 
	dist_ind = dist_ind60)
 
  #   Bundle data
jags.dat <- list(npts = npts120, nind = nind120, nocc = nocc120,
	theta = theta_double120, l = l_double120, 
	dist_ind = dist_ind120)
 
 
#   Inits function
jags.inits <- function(){list(alpha_v = runif(1, 0.01,  2), 
	alpha_lambda = runif(1, 0.01, 30), 
	alpha_rho = runif(1, 0.01, 1), 
	alpha_mu = runif(1, -3.14159265359, 3.14159265359))
	}

#   Parameters to monitor
jags.parms <- c("alpha_mu", "alpha_rho", "alpha_v","alpha_lambda", 
	"dre_dist_mu", "dre_dist_rho", "dre_dist_v", "dre_dist_lambda")


out <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "./DREM_dist_cov.txt", 
	n.chains = nc, 
	n.iter = ni, 
	n.burnin = nb,
	n.thin = nt, 
	parallel = TRUE,
	n.cores = 4)
	

S <- ggmcmc::ggs(out$samples) %>% dplyr::filter(Parameter != "deviance")

mu_plot <-  ggmcmc::ggs_density(S, family = "alpha_mu", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
mu_re_plot <-  ggmcmc::ggs_density(S, family = "dre_bear_mu", hpd = TRUE, greek = FALSE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
gridExtra::grid.arrange(mu_plot, mu_re_plot, ncol = 2)

rho_plot <-  ggmcmc::ggs_density(S, family = "alpha_rho", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
rho_re_plot <-  ggmcmc::ggs_density(S, family = "dre_bear_rho", hpd = TRUE, greek = FALSE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
gridExtra::grid.arrange(rho_plot, rho_re_plot, ncol = 2)

v_plot <-  ggmcmc::ggs_density(S, family = "alpha_v", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
v_re_plot <-  ggmcmc::ggs_density(S, family = "dre_bear_v", hpd = TRUE, greek = FALSE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
gridExtra::grid.arrange(v_plot, v_re_plot, ncol = 2)
	
lambda_plot <-  ggmcmc::ggs_density(S, family = "alpha_lambda", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
lambda_re_plot <-  ggmcmc::ggs_density(S, family = "dre_bear_lambda", hpd = TRUE, greek = FALSE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE)
gridExtra::grid.arrange(lambda_plot, lambda_re_plot, ncol = 2)
	
	
out_DREM_dist <- out
save(out_DREM_dist, file = "out_DREM_dist.Rdata")







