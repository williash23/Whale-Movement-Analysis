# Sara Williams
# 3/9/2016, last update: 10/28/2016; 1/17/2017
# Models movement states - script to run models.
#   Formualted after code from Morales et al. 2004.
#   Run in JAGS.
################################################################################



jags_wd <- "C:/Users/saraw/OneDrive/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/"
setwd(jags_wd)


#  Load packages
library(jagsUI)
library(mcmcplots)
library(coda)




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
jags.dat <- list(npts = npts_1, theta = theta, l = l, idx = bear_ind)
 
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
	n.cores = 4)

#out_rand_int_dist <- out
out_rand_int_bear <- out

S <- ggmcmc::ggs(out_rand_int_bear$samples) %>% dplyr::filter(Parameter != "deviance")

mu_plot <-  ggmcmc::ggs_density(S, family = "mu", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(-3, 3) +
	ylab("Density") +
	xlab(NULL)
rho_plot <-  ggmcmc::ggs_density(S, family = "rho", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0, 0.5) +
	ylab(NULL) +	
	xlab(NULL)
v_plot <-  ggmcmc::ggs_density(S, family = "v", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) +
	xlim(0.6, 1.2) +
	ylab("Density") +
	xlab("Value")
lambda_plot <-  ggmcmc::ggs_density(S, family = "lambda", hpd = TRUE, greek = TRUE) + 
	theme_bw() +
	guides(colour = FALSE, fill = FALSE) + 
	xlim(1.5, 4) +
	ylab(NULL) +
	xlab("Value")
	
gridExtra::grid.arrange(mu_plot, rho_plot, v_plot, lambda_plot)

ggmcmc::ci(S)
 # Parameter      low      Low  median  High  high
   # <fct>        <dbl>    <dbl>   <dbl> <dbl> <dbl>
 # 1 lambda[1]  2.13     2.17     2.34   2.51  2.55 
 # 2 lambda[2]  2.40     2.45     2.76   3.08  3.15 
 # 3 lambda[3]  2.34     2.40     2.74   3.12  3.19 
 # 4 mu[1]     -0.261   -0.206    0.0752 0.370 0.433
 # 5 mu[2]     -1.06    -0.769   -0.0278 0.652 0.921
 # 6 mu[3]     -2.89    -2.67    -0.577  2.55  2.85 
 # 7 rho[1]     0.170    0.187    0.274  0.356 0.371
 # 8 rho[2]     0.0250   0.0453   0.186  0.315 0.337
 # 9 rho[3]     0.00312  0.00632  0.0679 0.189 0.215
# 10 v[1]       0.943    0.953    1.01   1.06  1.07 
# 11 v[2]       0.848    0.862    0.934  1.01  1.02 
# 12 v[3]       0.776    0.790    0.863  0.937 0.953


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









#  Run "double" model - DREM
#   Discrete Random Effects Model (random effect on each parameter according to indicator for ship-whale context)
#   Bundle data
jags.dat <- list(npts = npts_1, nind = nind_1, nocc = nocc_1,
	theta = theta_double, l = l_double, 
	bear_ind = bear_ind)
 
#   Inits function
jags.inits <- function(){list(alpha_v = runif(1, 0.01,  2), 
	alpha_lambda = runif(1, 0.01, 30), 
	alpha_rho = runif(1, 0.01, 1), 
	alpha_mu = runif(1, -3.14159265359, 3.14159265359))
	}

#   Parameters to monitor
jags.parms <- c("alpha_mu", "alpha_rho", "alpha_v","alpha_lambda", 
	"dre_bear_mu", "dre_bear_rho", "dre_bear_v", "dre_bear_lambda")
	#"dre_dist_mu", "dre_dist_rho", "dre_dist_v", "dre_dist_lambda")


out <- jagsUI(jags.dat,
	jags.inits,
	jags.parms,
	model.file = "./DREM_cov.txt", 
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
	
	
out_DREM_bear <- out_DREM
save(out_DREM_dist, file = "out_DREM_bear.Rdata")



