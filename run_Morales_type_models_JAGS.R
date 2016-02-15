# Sara Williams
# 12/8/2015; updated 2/1/2016
# Model variety of CRWs after code from Morales et al. 2004.
#   Run in JAGS.
################################################################################

#  Load packages
library(R2jags)
library(dplyr)
library(adehabitatLT)		
library(mcmcplots)

#  Load and prep data
dat <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/Whales_0615_locations_clean.csv")

#  Generate step lengths and turning angles using ADEpackage
locs2_a <- arrange(dat, same_whale_ID, ob_order_time)
locs3_a<- dplyr::select(locs2_a, same_whale_ID, X_whale_UTM, Y_whale_UTM)
n_uni_a <- length(unique(locs3_a$same_whale_ID))
locs_tmp_a <- locs3_a %>% 
						distinct(same_whale_ID) %>%
						mutate(Name = seq(1,n_uni_a))
locs_tmp2_a <- full_join(locs3_a, locs_tmp_a, by = "same_whale_ID")				
locs_a <- locs_tmp2_a %>%
			  dplyr::select(Name, X_whale_UTM.x, Y_whale_UTM.x) %>%
			  dplyr::rename(X = X_whale_UTM.x, Y = Y_whale_UTM.x)
#   Create ltraj object				
whale_traj <- as.ltraj(xy = locs_a[,c("X","Y")], id = locs_a$Name, typeII = FALSE)			  
#   Convert into dataframe
traj_dat <- ld(whale_traj)			  
#   Connect traj_dat to dataframe with same_whale_ID
traj_dat_full <- cbind(traj_dat, dat)

#  Select only needed variables and rename
tmp <- traj_dat_full %>%
				 dplyr::select(id, x, y, dist, rel.angle, whale_dist_shore_m, ship_whale_dist) 
names(tmp) <- c("ID", "X", "Y", "steps", "turns", "shore_dist", "ship_dist")
tmp2 <- na.omit(tmp)
tmp3 <- filter(tmp2, steps < 5000)

#  Generate ordered ID variable
n_uni <- length(unique(tmp3$ID))
tmp4 <- tmp3 %>%
				distinct(ID) %>%
				mutate(ID_2 = seq(1, n_uni))
tmp5 <- full_join(tmp4, tmp3, by = "ID")			
obs <- tmp5 %>%
			 dplyr::select(ID, ID_2, X.x, Y.x, steps.x, turns.x, shore_dist.x, ship_dist.x) 
names(obs) <- c("ID", "ID_2", "X", "Y", "steps", "turns", "shore_dist", "ship_dist")

#  Data for models				
#   Indexing
npts <- nrow(obs)
raw_ID <- obs$ID_2
ID <- as.numeric(raw_ID)
#   Steps and turns
raw_l <- (obs$steps/1000)
l <-as.numeric(raw_l)
raw_theta <-  obs$turns
theta <- as.numeric(raw_theta)
#   Covariates
shore_raw <- as.numeric(obs$shore_dist)
ship_raw <- as.numeric(obs$ship_dist)
shore <-as.numeric(scale(shore_raw))
ship <- as.numeric(scale(ship_raw))
################################################################################

#   MCMC settings
nc <- 3					
ni <- 30000				
nb <- 10000				
nt <- 2			
################################################################################

#  Run "single" model

#   Bundle data
jags.data <- list("npts", "l", "theta")

#   Inits function
inits <- function(){list(v0=runif(1, 0.01, 10), 
											 lambda0=runif(1, 0.01, 10), 
											 rho0 = runif(1, 0.01, 1), 
											 mu0 = runif(1, -pi, pi))}

#   Parameters to estimate
params <- c("v0","lambda0", "mu0", "scale")

#   Unleash Gibbs sampler
out_single <- jags(data = jags.data, 
									 inits = inits, 
									 parameters.to.save = params, 
									 model.file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models/single.txt",
									 working.directory = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models",
									 n.thin = nt, 
									 n.chains = nc, 
									 n.burnin = nb, 
									 n.iter = ni)

print(out_single, dig = 3)
mcmcplot(out_single)

# sim_reps_b0_sing <- out_single$BUGS$sims.list$b0
# sim_reps_a0_sing <- out_single$BUGS$sims.list$a0
# sim_reps_mu0_sing <- out_single$BUGS$sims.list$mu0
# sim_reps_mean.wC_sing<- out_single$BUGS$sims.list$mean.wC
		
# save(sim_reps_b0_sing, file="sim_reps_b0_sing.RData")
# save(sim_reps_a0_sing, file="sim_reps_a0_sing.RData")
# save(sim_reps_mu0_sing, file="sim_reps_mu0_sing.RData")
# save(sim_reps_mean.wC_sing, file="sim_reps_mean.wC_sing.RData")
####################################################################################################

#  Run "double" model

#   Bundle data
jags.dat <- list("npts", "l", "theta", "ID")

#   Inits function
inits <- function(){list(v=runif(2, 0.01, 5), 
											 lambda=c(NA,runif(1, 0.01, 5)),
											 eps=runif(1, 0.01, 5),											 
											 rho=runif(2, 0.01, 1), 
											 mu=runif(2, -pi, pi),
											 beta0=runif(1, -5, 5))
											 }

#   Parameters to estimate
params <- c("v","lambda", "mu", "rho", "scale", "eps", "mean.q", "beta0")

#   Unleash Gibbs sampler
out_double <- jags(data = jags.dat, 
									    inits = inits, 
										parameters.to.save = params, 
										model.file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models/double.txt",
										working.directory = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models",
										n.thin = nt, 
										n.chains = nc, 
										n.burnin = nb, 
										n.iter = ni)

print(out_double, dig = 3)
mcmcplot(out_double)

# sim_reps_v1_doub <- out_double$BUGS$sims.list$v[1]
# sim_reps_v2_doub_sw <- out_double_sw$BUGS$sims.list$v[2]
# sim_reps_scale1_doub <- out_double$BUGS$sims.list$scale[1]
# sim_reps_scale2_doub <- out_double$BUGS$sims.list$scale[2]
# sim_reps_mu1_doub <- out_double$BUGS$sims.list$mu[1]
# sim_reps_mu2_doub <- out_double$BUGS$sims.list$mu[2]
# sim_reps_mean.idx_doub <- out_double$BUGS$sims.list$mean.idx
		
# save(sim_reps_v1_doub, file="sim_reps_v1_doub.RData")
# save(sim_reps_v2_doub, file="sim_reps_v2_doub.RData")
# save(sim_reps_scale1_doub, file="sim_reps_scale1_doub.RData")
# save(sim_reps_scale1_doub, file="sim_reps_scale1_doub.RData")
# save(sim_reps_mu1_doub, file="sim_reps_mu1_doub.RData")
# save(sim_reps_mu2_doub, file="sim_reps_mu2_doub.RData")
# save(sim_reps_mean.idx_doub, file="sim_reps_mean.idx_doub.RData")
####################################################################################################

#  Run "double switch" model

#   Bundle data
jags.dat <- list("npts", "l", "theta", "ID")

#   Inits function
inits <- function(){list(v=runif(2, 0.01,  5), 
											 lambda=c(NA, runif(1, 0.01, 5)), 
											 eps=runif(1, 0.01, 5),
											 rho = runif(2, 0.01, 1), 
											 mu = runif(2, -pi, pi),
											 phi = c(runif(1, 0.01, 1), NA),
											 beta0=runif(2, -5, 5),
											 idx = c(NA, rep(1, length.out = (npts-1))))
											 }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "scale", "beta0", "mean.q")

#   Unleash Gibbs sampler
out_double_sw <- jags(data = jags.dat, 
											inits = inits, 
											parameters.to.save = params, 
											model.file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models/double_switch.txt",
											working.directory = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models",
											n.thin = nt, 
											n.chains = nc, 
											n.burnin = nb, 
											n.iter = ni)

print(out_double_sw, dig = 3)
mcmcplot(out_double_sw)

# sim_reps_v1_doub_sw <- out_double_sw$BUGS$sims.list$v[1]
# sim_reps_v2_doub_sw <- out_double_sw$BUGS$sims.list$v[2]
# sim_reps_scale1_doub_sw <- out_double_sw$BUGS$sims.list$scale[1]
# sim_reps_scale2_doub_sw <- out_double_sw$BUGS$sims.list$scale[2]
# sim_reps_mu1_doub_sw <- out_double_sw$BUGS$sims.list$mu[1]
# sim_reps_mu2_doub_sw <- out_double_sw$BUGS$sims.list$mu[2]
# sim_reps_mean.idx_doub_sw <- out_double_sw$BUGS$sims.list$mean.idx
# sim_reps_qu1_doub_sw <- out_double_sw$BUGS$sims.list$qu[1]
# sim_reps_qu2_doub_sw <- out_double_sw$BUGS$sims.list$qu[2]
		
# save(sim_reps_v1_doub_sw, file="sim_reps_v1_doub_sw.RData")
# save(sim_reps_v2_doub_sw, file="sim_reps_v2_doub_sw.RData")
# save(sim_reps_scale1_doub_sw, file="sim_reps_scale1_doub_sw.RData")
# save(sim_reps_scale2_doub_sw, file="sim_reps_scale2_doub_sw.RData")
# save(sim_reps_mu1_doub_sw, file="sim_reps_mu1_doub_sw.RData")
# save(sim_reps_mu2_doub_sw, file="sim_reps_mu2_doub_sw.RData")
# save(sim_reps_mean.idx_doub_sw, file="sim_reps_mean.idx_doub_sw.RData")
# save(sim_reps_qu1_doub_sw, file="sim_reps_qu1_doub_sw.RData")
# save(sim_reps_qu2_doub_sw, file="sim_reps_qu2_doub_sw.RData")
####################################################################################################

#  Run "double  covariate" model

#   Bundle data
jags.dat <- list("npts", "l", "theta", "shore", "ID")

#   Inits function
inits <- function(){list(v=runif(2, 0.01, 5), 
											 lambda=c(NA,runif(1, 0.01, 5)),
											 eps=runif(1, 0.01, 5),											 
											 rho=runif(2, 0.01, 1), 
											 mu=runif(2, -pi, pi),
											 beta0=runif(1, -5, 5),
											 beta1=runif(1, -5, 5))
											 }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "scale", "beta0", "beta1", "mean.q")

#   Unleash Gibbs sampler
out_doub_cov <- jags(data = jags.dat, 
											inits = inits, 
											parameters.to.save = params, 
											model.file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models/double_cov.txt",
											working.directory = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models",
											n.thin = nt, 
											n.chains = nc, 
											n.burnin = nb, 
											n.iter = ni)

print(out_doub_cov, dig = 3)
mcmcplot(out_doub_cov)

# sim_reps_mean.nu_h_doub_cov <- out_doub_cov$BUGS$sims.list$mean.nu_h
# sim_reps_v1_doub_cov <- out_double_cov$BUGS$sims.list$v[1]
# sim_reps_v2_doub_cov <- out_double_cov$BUGS$sims.list$v[2]
# sim_reps_scale1_doub_cov <- out_double_cov$BUGS$sims.list$scale[1]
# sim_reps_scale2_doub_cov <- out_double_cov$BUGS$sims.list$scale[2]
# sim_reps_mu1_doub_cov <- out_double_cov$BUGS$sims.list$mu[1]
# sim_reps_mu2_doub_cov <- out_double_cov$BUGS$sims.list$mu[2]
# sim_reps_beta0_doub_cov <- out_double_cov$BUGS$sims.list$beta0
# sim_reps_beta1_doub_cov <- out_double_cov$BUGS$sims.list$beta1
# sim_reps_beta2_doub_cov <- out_double_cov$BUGS$sims.list$beta2

# save(sim_reps_v1_doub_cov, file="sim_reps_v1_doub_cov.RData")	
# save(sim_reps_v2_doub_cov, file="sim_reps_v2_doub_cov.RData")		
# save(sim_reps_scale1_doub_cov, file="sim_reps_scale1_doub_cov.RData")
# save(sim_reps_scale2_doub_cov, file="sim_reps_scale2_doub_cov.RData")
# save(sim_reps_mu1_doub_cov, file="sim_reps_mu1_doub_cov.RData")
# save(sim_reps_mu2_doub_cov, file="sim_reps_mu2_doub_cov.RData")
# save(sim_reps_beta0_doub_cov, file="sim_reps_beta0_doub_cov.RData")
# save(sim_reps_beta1_doub_cov, file="sim_reps_beta1_doub_cov.RData")
# save(sim_reps_beta2_doub_cov, file="sim_reps_beta2_doub_cov.RData")
# save(sim_reps_alpha_doub_cov, file="sim_reps_alpha_doub_cov.RData")
####################################################################################################

#  Run "double switch with covariate" model

#   Bundle data
jags.dat <- list("npts", "l", "theta", "ship", "ID")

#   Inits function
inits <- function(){list(v=runif(2, 0.01,  5), 
											 lambda=c(NA, runif(1, 0.01, 5)), 
											 eps=runif(1, 0.01, 5),
											 rho = runif(2, 0.01, 1), 
											 mu = runif(2, -pi, pi),
											 phi = c(runif(1, 0.01, 1), NA),
											 idx = c(NA, rep(1, length.out = (npts-1))),
											 beta0=runif(2, -5, 5),
											 beta1=runif(2, -5, 5))
											 }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "scale", "beta0", "beta1", "ship_eff", "mean.q")

#   Unleash Gibbs sampler
out_double_sw_cov <- jags(data = jags.dat, 
											inits = inits, 
											parameters.to.save = params, 
											model.file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models/double_switch_cov.txt",
											working.directory = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models",
											n.thin = nt, 
											n.chains = nc, 
											n.burnin = nb, 
											n.iter = ni)

print(out_double_sw_cov, dig = 3)
mcmcplot(out_double_sw_cov)

# sim_reps_v1_doub_sw <- out_double_sw$BUGS$sims.list$v[1]
# sim_reps_v2_doub_sw <- out_double_sw$BUGS$sims.list$v[2]
# sim_reps_scale1_doub_sw <- out_double_sw$BUGS$sims.list$scale[1]
# sim_reps_scale2_doub_sw <- out_double_sw$BUGS$sims.list$scale[2]
# sim_reps_mu1_doub_sw <- out_double_sw$BUGS$sims.list$mu[1]
# sim_reps_mu2_doub_sw <- out_double_sw$BUGS$sims.list$mu[2]
# sim_reps_mean.idx_doub_sw <- out_double_sw$BUGS$sims.list$mean.idx
# sim_reps_qu1_doub_sw <- out_double_sw$BUGS$sims.list$qu[1]
# sim_reps_qu2_doub_sw <- out_double_sw$BUGS$sims.list$qu[2]
		
# save(sim_reps_v1_doub_sw, file="sim_reps_v1_doub_sw.RData")
# save(sim_reps_v2_doub_sw, file="sim_reps_v2_doub_sw.RData")
# save(sim_reps_scale1_doub_sw, file="sim_reps_scale1_doub_sw.RData")
# save(sim_reps_scale2_doub_sw, file="sim_reps_scale2_doub_sw.RData")
# save(sim_reps_mu1_doub_sw, file="sim_reps_mu1_doub_sw.RData")
# save(sim_reps_mu2_doub_sw, file="sim_reps_mu2_doub_sw.RData")
# save(sim_reps_mean.idx_doub_sw, file="sim_reps_mean.idx_doub_sw.RData")
# save(sim_reps_qu1_doub_sw, file="sim_reps_qu1_doub_sw.RData")
# save(sim_reps_qu2_doub_sw, file="sim_reps_qu2_doub_sw.RData")
####################################################################################################

## model "Switch with covariates" used in:
## Extracting More out of Relocation Data: Building Movement Models as Mixtures of  Random Walks
## Juan Manuel Morales, Daniel T. Haydon, Jacqui Frair, Kent E. Holsinger and John M. Fryxell 
## contact juan.morales@uconn.edu

model{

   ## priors

   b[1] ~ dgamma(0.01,0.01) 	## shape parameter for fast movement
   b[2] ~ dgamma(0.01,0.01) 	## shape parameter for slow movement

   a[2] ~ dgamma(0.01, 0.01)	## scale parameter for slow movement
   eps ~ dnorm(0.0, 0.01)I(0.0,)	## a nonnegative variate
   a[1]  <- a[2] + eps		## scale parameter for fast movement

   rho[1] ~ dunif(0,1)		## mean cosine of turns for fast movement
   rho[2] ~ dunif(0,1)		## mean cosine of turns for slow movement

   mu[1] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for fast movement 
   mu[2] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for slow movement

   for(i in 1:10){
      for(j in 1:2){
         m[j,i] ~ dnorm(0,0.1)	# coefficients to relate distance to habitat i to switching rate
      }
   }	

   for (i in 1:10){
      m[1,i] <- 0
   }

   beta[1] ~ dnorm(0,0.1)       # intercept
   beta[2] ~ dnorm(0,0.1)

   phi[1] ~ dunif(0,1)
   phi[2] ~ 1-phi[1]
   idx[1] ~ dcat(phi[])		## asign state for first observation 
  
   Pi <- 3.14159265359		## define pi


   for (t in 2:npts) {

      logit.q[t] <- exp(beta[idx[t-1]] + m[idx[t-1],1]*water[t]/10000 + m[idx[t-1],2]*swamp[t]/10000 + m[idx[t-1],3]*otw[t]/10000 + m[idx[t-1],4]*openfor[t]/10000 + m[idx[t-1],5]*ntw[t]/10000 + m[idx[t-1],6]*mixfor[t]/10000 + m[idx[t-1],7]*dev[t]/10000 + m[idx[t-1],8]*ddf[t]/10000 + m[idx[t-1],9]*conif[t]/10000 + m[idx[t-1],10]*alvar[t]/10000)

      q[t] <- logit.q[t]/(1 + logit.q[t])

      nu[t,1] <- q[idx[t-1]]
      nu[t,2] <- 1-q[idx[t-1]]
      idx[t] ~ dcat(nu[t,])   ##  idx is the latent variable and the parameter index

      ## likelihood for steps
      l[t] ~ dweib(b[idx[t]], a[idx[t]])	# Weibull distriution for step length

      ## likelihood for turns.  
      ## use the “ones” trick (see WinBUGS manual) to sample from the Wrapped Cauchy distribution

      ones[t] <- 1
      ones[t] ~ dbern(wC[t])
      ## below is the pdf for Wrapped Cauchy distribution, divided by 500 (arbitrary) to ensure that wC[t] will be less than one
      wC[t] <- ( 1/(2*Pi)*(1-rho[idx[t]]*rho[idx[t]])/(1+rho[idx[t]]*rho[idx[t]]-2*rho[idx[t]]*cos(theta[t]-mu[idx[t]t])) )/500

  }
}
####################################################################################################
#  Run "triple switch" model

#   Bundle data
jags.dat <- list("npts", "l", "theta")

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "q", "qq", "scale")

#   Inits function
inits <- function(){list(v=runif(3, 0.01,  5), 
											 lambda=c(NA, NA, runif(1, 0.01, 5)), 
											 eps=runif(2, 0.01, 5),
											 rho = runif(3, 0.01, 1), 
											 mu = runif(3, -pi, pi),
											 qu = runif(3, 0.01, 1),
											 qq=runif(3, 0.01, 1),
											 #phi = c(runif(2, 0.01, 1), NA),
											 idx = c(NA, rep(1, length.out = (npts-1))))
											 }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "qu", "qq", "scale", "mean.idx")

#   Unleash Gibbs sampler
out_triple_sw <- jags(data = jags.dat, 
											inits = inits, 
											parameters.to.save = params, 
											model.file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models/triple_switch.txt",
											working.directory = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models",
											n.thin = nt, 
											n.chains = nc, 
											n.burnin = nb, 
											n.iter = ni)

print(out_triple_sw, dig = 3)
mcmcplot(out_triple_sw)

## # phi[1] <- 0.2
	#phi[2] <- 0..5
	#phi[3] <- 0.3
####################################################################################################


## model "Triple switch" used in:
## Extracting More out of Relocation Data: Building Movement Models as Mixtures of  Random Walks
## Juan Manuel Morales, Daniel T. Haydon, Jacqui Frair, Kent E. Holsinger and John M. Fryxell 
## contact juan.morales@uconn.edu

model{

   ## priors
   ## shape parameters for step length
   b[1] ~ dgamma(0.01,0.01)
   b[2] ~ dgamma(0.01,0.01)
   b[3] ~ dgamma(0.01,0.01) 

   eps1 ~  dnorm(0.0, 0.01)I(0.0,)
   eps2 ~ dnorm(0.0, 0.01)I(0.0,)
   ## scale parameters for step length
   a[3] ~ dgamma(0.01, 0.01) 
   a[2] <- a[3] + eps1 
   a[1] <- a[2] + eps2

   ## mean cosine for turns
   rho[1] ~ dunif(0,1)
   rho[2] ~ dunif(0,1)
   rho[3] ~ dunif(0,1)

   ## mean direction for turns
   mu[1] ~ dunif(-3.14159265359, 3.14159265359) 
   mu[2] ~ dunif(-3.14159265359, 3.14159265359)
   mu[3] ~ dunif(-3.14159265359, 3.14159265359)


   ##  priors for the probability of switching from anything to 1
   qq[1] ~ dunif(0,1)
   qq[2] ~ dunif(0,1)
   qq[3] ~ dunif(0,1)
   q[1] ~ dunif(0,1)
   q[2] ~ dunif(0,1)
   q[3] ~ dunif(0,1)
   ## asign state for first observation 
   idx[1] ~ dcat(phi[])		
  
   Pi <- 3.14159265359		## define pi


   for (t in 2:npts) {

      nu[t,1] <- q[idx[t-1]]
      nu[t,2] <- (1 -q [idx[t-1]] ) * qq[idx[t-1]] 
      nu[t,3] <- (1 -q [idx[t-1]] ) * (1-qq[idx[t-1]] )

      idx[t] ~ dcat(nu[t,])   ##  idx is the latent variable and the parameter index

      ## likelihood for steps
      l[t] ~ dweib(b[idx[t]], a[idx[t]])	# Weibull distriution for step length

      ## likelihood for turns.  
      ## use the “ones” trick (see WinBUGS manual) to sample from the Wrapped Cauchy distribution

      ones[t] <- 1
      ones[t] ~ dbern(wC[t])
      ## below is the pdf for Wrapped Cauchy distribution, divided by 500 (arbitrary) to ensure that wC[t] will be less than one
      wC[t] <- ( 1/(2*Pi)*(1-rho[idx[t]]*rho[idx[t]])/(1+rho[idx[t]]*rho[idx[t]]-2*rho[idx[t]]*cos(theta[t]-mu[idx[t]t])) )/500

  }
}
####################################################################################################

#  Hierarchical switch model
model{

# Hyperpriors. These are priors on the parameters of the prior distributions for individual
# level parameters.

# Amu.h and Atau.h are priors for the parameters of the Normal distribution used as prior
# for scale parameter in step length during fast movement. eps.h and epstau.h are priors 
# for the censored Normal used for the difference between scale parameter of step length
# in fast and slow movement. 

Amu.h ~  dnorm(0.0,0.01)I(0,)
Atau.h ~ dnorm(0.0,0.01)I(0,)
eps.h ~ dnorm(0,0.01)I(0,)
epstau.h ~ dnorm(0,0.01)I(0,)

# priors on the parameters of the Normal distributions used as priors on the shape 
# parameters of the Weibul distributions used to model step length
Bmu.h[1] ~ dnorm(0.0,0.01)I(0,)
Bmu.h[2] ~ dnorm(0.0,0.01)I(0,)
Btau.h[1] ~ dnorm(0.0,0.01)I(0,)
Btau.h[2] ~ dnorm(0.0,0.01)I(0,)

# priors on the parameters of the Normal distributions used as priors on the mean
# direction of turning angles
mumean.h[1] ~ dunif(0,6.28318530717959)
mumean.h[2] ~ dunif(0,6.28318530717959)
mutau.h[1] ~ dnorm(0,0.01)I(0.0,)
mutau.h[2] ~ dnorm(0,0.01)I(0.0,)

# priors on the parameters of the Beta distribution used as priors on the mean 
# cosine of turning angles
arho.h[1] ~ dnorm(0,0.001)I(0,)
arho.h[2] ~ dnorm(0,0.001)I(0,)
brho.h[1] ~ dnorm(0,0.001)I(0,)
brho.h[2] ~ dnorm(0,0.001)I(0,)

# priors on the parameters of the Beta distribution used as priors on switching rates
qa.h[1] ~ dnorm(0,0.01)I(0,)
qa.h[2] ~ dnorm(0,0.01)I(0,)
qb.h[1] ~ dnorm(0,0.01)I(0,)
qb.h[2] ~ dnorm(0,0.01)I(0,)

Pi <- 3.14159265359   # define pi
# assign initial state
for(i in 1:npaths){
idx[1,i]~dcat(phi[])
}

# iterate over movement paths
	for(k in 1:npaths){
		# individual level priors
		A[k,2] ~ dnorm(Amu.h, Asigma.h)
		eps[k] ~ dnorm(eps.h, epstau.h)I(0,)
		A[k,1] <- A[k,2] + eps[k]
		B[k,1] ~ dnorm(Bmu.h[1], Btau.h[1]) I(0,)
		B[k,2] ~ dnrom(Bmu.h[2], Btau.h[2]) I(0,)
		RHO[k,1] ~ dbeta(arho.h[1], brho.h[1])
		RHO[k,2] ~ dbeta(arho.h[2], brho.h[2])
		MU[k,1] ~ dnorm(mumean.h[1], mutau.h[1])
		MU[k,2] ~ dnorm(mumean.h[2], mutau.h[2])
		Q[k,1] ~ dbeta(qa.h[1],qb.h[1])
		Q[k,2] ~ dbeta(qa.h[2],qb.h[2])

# iterate over observations
		for (t in 2:npts) {
			l[t,k] ~ dweib(B[k,idx[t,k]], A[k,idx[t,k]]) 

# use 'ones'trick to simulate the wrapped Cauchy distribution
			ones[t,k] <- 1
			ones[t,k] ~ dbern( wc[t,k] )
			wc[t,k] <- (1/(2*Pi)*(1-rho[t,k]*rho[t,k])/(1+rho[t,k]*rho[t,k]-2*rho[t,k]*cos(theta[t,k]-mu.t[t,k])))/ 300

			theta[t,k] ~ dunif(0,6.28318530717959)
			rho[t,k]<-RHO[k,idx[t,k]]
			mu.t[t,k]<- MU[k,idx[t,k]]
			# the probability of being in movement type 1
			p[t,k,1] <- Q[k,idx[t-1,k]]
			p[t,k,2] <- 1-Q[k,idx[t-1,k]]
			idx[t,k] ~ dcat(p[t,k,])
		}
	}
}




####################################################################################################

#  JAGS parameterization of weibull distribution.
#   dweib(v = shape, lambda)
#  R parameterization of weibull distribution.
#   rweibull(n, shape, scale)
#  Transform lambda to R scale parameter.
#    scale[1] <- (1/lambda[1])^ (1/v[1])
#  Find mean step length.
#   mean <- scale*gamma((1+1/shape))

#  Plots.
n <- 1000000

#  Step lengths.

#  Simulate step lengths.
#   Uses parameters estimated from Morales et al. code.
sim_steps1 <- as.data.frame(rweibull(n, 1.063, 600))
sim_steps2 <- as.data.frame(rweibull(n, 1.210,  300))
names(sim_steps1)[1] <- "dist"
names(sim_steps2)[1] <- "dist"
#   Create dataframe of simulated movement 1 and movement 2 steps.
sim_steps <- cbind(sim_steps1, sim_steps2)
names(sim_steps)[1] <- "1"
names(sim_steps)[2] <- "2"
#  Look at density plots of both on same figure.
df_steps_sim <- melt(sim_steps)
compare_steps_sim <- ggplot(df_steps_sim) + 
											   geom_density(aes(x = value, colour = variable, fill = variable), alpha = .25) +
											   xlim(0, 5000) +
											   xlab("Step length (m)") +
											   theme_bw()
 compare_steps_sim


rhoD <- est.rho(turns.dive)
muD<- circ.mean(turns.dive)

sim_turns1 <- rwrpcauchy(n, muD, rhoD)


data{
	#  Generate vector of 1's for 1's trick for wrapped Cauchy distribution
	for(i in 2:npts){
		ones[i] <- 1
		}	
	}
