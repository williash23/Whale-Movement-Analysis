# Sara Williams
# 10/16/2015
#  Model single CRW - as Morales et al. 2004
################################################################################

## model "Single" used in:
## Extracting More out of Relocation Data: Building Movement Models as Mixtures of  Random Walks
## Juan Manuel Morales, Daniel T. Haydon, Jacqui Frair, Kent E. Holsinger and John M. Fryxell 
## contact juan.morales@uconn.edu

model{

   ## priors

   b ~ dgamma(0.01,0.01) 	## shape parameter for step length distribution
   a ~ dgamma(0.01, 0.01)	## scale parameter for step length distribution

   rho ~ dunif(0,1)		## mean cosine of turns
   mu ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns
  
   Pi <- 3.14159265359		## define pi


   for (t in 1:npts) {

      ## likelihood for steps
      l[t] ~ dweib(b[t], a[t])	# Weibull distriution for step length

      ## likelihood for turns.  
      ## use the “ones” trick (see WinBUGS manual) to sample from the Wrapped Cauchy distribution

      ones[t] <- 1
      ones[t] ~ dbern(wC[t])
      ## below is the pdf for Wrapped Cauchy distribution, divided by 500 (arbitrary) to ensure that wC[t] will be less than one
      wC[t] <- ( 1/(2*Pi)*(1-rho[t]*rho[t])/(1+rho[t]*rho[t]-2*rho[t]*cos(theta[t]-mu[t])) )/500

  }
}


## model "Double" used in:
## Extracting More out of Relocation Data: Building Movement Models as Mixtures of  Random Walks
## Juan Manuel Morales, Daniel T. Haydon, Jacqui Frair, Kent E. Holsinger and John M. Fryxell 
## contact juan.morales@uconn.edu
npts <- length(steps)
l <- steps

sink("double.txt")
cat("
model{
   ## priors

   b[1] ~ dgamma(0.01,0.01) 	## shape parameter for slow movement
   b[2] ~ dgamma(0.01,0.01) 	## shape parameter for fast movement

   a[2] ~ dgamma(0.01, 0.01)	## scale parameter for fast movement
   eps ~ dnorm(0.0, 0.01)I(0.0,)	## a nonnegative variate
   a[1]  <- a[2] + eps		## scale parameter for slow movement

   rho[1] ~ dunif(0,1)		## mean cosine of turns for slow movement
   rho[2] ~ dunif(0,1)		## mean cosine of turns for fast movement

   mu[1] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for slow movement 
   mu[2] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for fast movement
  
   Pi <- 3.14159265359		## define pi

   for (t in 1:npts) {

      nu[t,1] ~ dunif(0,1)    ## probability of being in movement state 1 at time t
      nu[t,2] <- 1 - nu[t,1]
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

",fill=TRUE)
sink()


#   Bundle data
win.data <- list("npts", "l")

#   Inits function
inits <- function(){ list(b=rnorm(1), a=rnorm(1), rho = rlnorm(1), mu = rlnorm(1), nu = rlnorm(1))}

#   Parameters to estimate
params <- c("b","a", "rho", "mu")

#   MCMC settings
nc <- 3					# Number of chains
ni <- 2000				# Number of draws from posterior for each chain
nb <- 500				# Number of draws to discard as burn-in
nt <- 1					# Thinning rate

#   Unleash Gibbs sampler
out <- bugs(data = win.data, inits = inits, parameters = params, model = "double.txt",
n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, debug = TRUE)

print(out, dig = 3)