# Sara Williams
# 10/16/2015; updated 12/7/2015
# Model variety of CRWs after code from Morales et al. 2004
#   Run in WinBUGS
################################################################################
### model "Double with covariates" used in:
## Extracting More out of Relocation Data: Building Movement Models as Mixtures of  Random Walks
## Juan Manuel Morales, Daniel T. Haydon, Jacqui Frair, Kent E. Holsinger and John M. Fryxell 
## contact juan.morales@uconn.edu

sink("double_cov.txt")
cat("
model{
   ## priors

  b[1] ~ dgamma(0.01,0.01) 	## shape parameter for slow movement
  b[2] ~ dgamma(0.01,0.01) 	## shape parameter for fast movement

  a[2] ~ dgamma(0.01, 0.01)	## scale parameter for fast movement
  eps ~ dnorm(0.0, 0.01)I(0.0,)	## a nonnegative variate
  a[1]  <- a[2] + eps		## scale parameter for slow movement

  rho[1] ~ dunif(0,1)		## mean cosine of turns for slow movement
  rho[2] ~ dunif(0,1)		## mean cosine of turns for slow movement

  mu[1] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for slow movement 
  mu[2] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for fast movement

  ## priors for habitat types
  for (i in 1:10) {
  mu.phi[i] ~ dnorm(0.0, 0.01)
  }
    
  Pi <- 3.14159265359		## define pi
		for (t in 1:npts) {

		## movement state is related to current habitat type
		mu.type[t] <- mu.phi[typ[t]] ## typ is a variable that indicates habitat type at current location
		logit.nu[t] ~ dnorm(mu.type[t], 0.01)
		nu_h[t] <- exp(logit.nu[t])/(1 + exp(logit.nu[t]))  ## probability of being in movement type 1
		nu[t,1] <- nu_h[t]
		nu[t,2] <- 1 - nu_h[t]
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
   mean.wC <- mean(wC[])
   mean.nu_h <- mean(nu_h[])
   
 }

",fill=TRUE)
sink()

#   Bundle data
win.data <- list("npts", "l", "theta")

#   Inits function
inits <- function(){list(b=runif(1, 0.01, 10), 
											 a=runif(1, 0.01, 10), 
											 rho = runif(1, 0.01, 1), 
											 mu = runif(1, -pi, pi), 
											 eps=rnorm(1),
											 mu.ship=rnorm(1))}

#   Parameters to estimate
params <- c("b","a", "mu", "mean.wC", "mean.nu_h")

#   MCMC settings
nc <- 3					# Number of chains
ni <- 30000				# Number of draws from posterior for each chain
nb <- 15000				# Number of draws to discard as burn-in
nt <- 2					# Thinning rate

#   Unleash Gibbs sampler
out <- bugs(data = win.data, inits = inits, parameters = params, model = "double_cov.txt",
n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, debug = TRUE)

print(out, dig = 3)
####################################################################################################

## model "Double switch" used in:
## Extracting More out of Relocation Data: Building Movement Models as Mixtures of  Random Walks
## Juan Manuel Morales, Daniel T. Haydon, Jacqui Frair, Kent E. Holsinger and John M. Fryxell 
## contact juan.morales@uconn.edu

sink("double_sw.txt")
cat("
model{

   ## priors

   b[1] ~ dgamma(0.01,0.01) 	## shape parameter for slow movement
   b[2] ~ dgamma(0.01,0.01) 	## shape parameter for fast movement

   a[2] ~ dgamma(0.01, 0.01)	## scale parameter for fast movement
   eps ~ dnorm(0.0, 0.01)I(0.0,)	## a nonnegative variate
   a[1]  <- a[2] + eps		## scale parameter for slow movement

   rho[1] ~ dunif(0,1)		## mean cosine of turns for slow movement
   rho[2] ~ dunif(0,1)		## mean cosine of turns for slow movement

   mu[1] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for slow movement 
   mu[2] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for fast movement

   q[1] ~ dunif(0,1)		## probability of being in state 1 at t given that individual was in state 1 at time t-1
   q[2] ~ dunif(0,1)		## probability of being in state 1 at t given that individual was in state 2 at time t-1

   phi[1] ~ dunif(0,1)
   phi[2] <- 1-phi[1]
   idx[1] ~ dcat(phi[])		## asign state for first observation 
  
   Pi <- 3.14159265359		## define pi


   for (t in 2:npts) {

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
      wC[t] <- ( 1/(2*Pi)*(1-rho[idx[t]]*rho[idx[t]])/(1+rho[idx[t]]*rho[idx[t]]-2*rho[idx[t]]*cos(theta[t]-mu[idx[t]])) )/500

	}
}
",fill=TRUE)
sink()

#   Bundle data
win.data <- list("npts", "l", "theta")

#   Inits function
inits <- function(){list(b=runif(1, 0.01, 10), 
											 a=runif(1, 0.01, 10), 
											 rho = runif(1, 0.01, 1), 
											 mu = runif(1, -pi, pi), 
											 eps=rnorm(1),
											 q=runif(1, 0.01, 1))}

#   Parameters to estimate
params <- c("b","a", "mu", "rho", "q")

#   MCMC settings
nc <- 3					# Number of chains
ni <- 30000				# Number of draws from posterior for each chain
nb <- 15000				# Number of draws to discard as burn-in
nt <- 2					# Thinning rate

#   Unleash Gibbs sampler
out_doub_sw <- bugs(data = win.data, inits = inits, parameters = params, model = "double_sw.txt",
n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, debug = TRUE)

print(out_doub_sw, dig = 3)
mcmcplot(out_doub_sw)

sim_reps_b_doub_sw <- out_doub_sw$BUGS$sims.list$b
sim_reps_a_doub_sw <- out_doub_sw$BUGS$sims.list$a
sim_reps_mu_doub_sw <- out_doub_sw$BUGS$sims.list$mu
sim_reps_mean.wC_doub_sw<- out_doub_sw$BUGS$sims.list$mean.wC
sim_reps_mean.nu1_doub_sw<- out_doub_sw$BUGS$sims.list$mean.nu1
sim_reps_mean.nu2_doub_sw<- out_doub_sw$BUGS$sims.list$mean.nu2
		
save(sim_reps_b_doub_sw, file="sim_reps_b0_doub_sw.RData")
save(sim_reps_a_doub_sw, file="sim_reps_a0_doub_sw.RData")
save(sim_reps_mu_doub_sw, file="sim_reps_mu0_doub_sw.RData")
save(sim_reps_mean.wC_doub_sw, file="sim_reps_mean.wC_doub_sw.RData")
save(sim_reps_mean.nu1_doub_sw, file="sim_reps_mean.nu1_doub_sw.RData")
save(sim_reps_mean.nu2_doub_sw, file="sim_reps_mean.nu2_doub_sw.RData")
####################################################################################################

## model "Switch with covariates" used in:
## Extracting More out of Relocation Data: Building Movement Models as Mixtures of  Random Walks
## Juan Manuel Morales, Daniel T. Haydon, Jacqui Frair, Kent E. Holsinger and John M. Fryxell 
## contact juan.morales@uconn.edu

model{

   ## priors

   b[1] ~ dgamma(0.01,0.01) 	## shape parameter for slow movement
   b[2] ~ dgamma(0.01,0.01) 	## shape parameter for fast movement

   a[2] ~ dgamma(0.01, 0.01)	## scale parameter for fast movement
   eps ~ dnorm(0.0, 0.01)I(0.0,)	## a nonnegative variate
   a[1]  <- a[2] + eps		## scale parameter for slow movement

   rho[1] ~ dunif(0,1)		## mean cosine of turns for slow movement
   rho[2] ~ dunif(0,1)		## mean cosine of turns for slow movement

   mu[1] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for slow movement 
   mu[2] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for fast movement

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

## model "Double switch" used in:
## Extracting More out of Relocation Data: Building Movement Models as Mixtures of  Random Walks
## Juan Manuel Morales, Daniel T. Haydon, Jacqui Frair, Kent E. Holsinger and John M. Fryxell 
## contact juan.morales@uconn.edu

model{

   ## priors

   b[1] ~ dgamma(0.01,0.01) 	## shape parameter for slow movement
   b[2] ~ dgamma(0.01,0.01)I(1.1,) 	## shape parameter for fast movement forced to be greater than 1.1

   a[2] ~ dgamma(0.01, 0.01)	## scale parameter for fast movement
   eps ~ dnorm(0.0, 0.01)I(0.0,)	## a nonnegative variate
   a[1]  <- a[2] + eps		## scale parameter for slow movement

   rho[1] ~ dunif(0,1)		## mean cosine of turns for slow movement
   rho[2] ~ dunif(0,1)		## mean cosine of turns for slow movement

   mu[1] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for slow movement 
   mu[2] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for fast movement

   q[1] ~ dunif(0,1)		## probability of being in state 1 at t given that individual was in state 1 at time t-1
   q[2] ~ dunif(0,1)		## probability of being in state 1 at t given that individual was in state 2 at time t-1

   phi[1] ~ dunif(0,1)
   phi[2] ~ 1-phi[1]
   idx[1] ~ dcat(phi[])		## asign state for first observation 
  
   Pi <- 3.14159265359		## define pi


   for (t in 2:npts) {

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