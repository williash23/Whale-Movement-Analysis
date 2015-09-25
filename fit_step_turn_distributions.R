#  Sara Williams
#  9/3/2015
#  Invesitgate step length and turning angle distributions, generate simulated data from
#   observered data distribution parameters, and practice CRW.
################################################################################

require(vcd)
require(MASS)
require(fitdistrplus)
require(CircStats)
require(qualityTools)

#  Load data
#   Use data.out and turning angles and step lengths from  
#   "Movement_parameters_analysis_script"

#  Use package MASS fitdistr() function or package fitdistrplus fitdist() function.
fit_steps_exp<- fitdist(all.steps, "exp", method = "mle", start = list(rate=0.005)) 
#fit_steps_gam <- fitdist(all.steps, "gamma", method = "mle") 
fit_steps_wei <- fitdist(all.steps, "weibull",  method = "mle") 
fit_steps_lnorm<- fitdist(all.steps, "lnorm", method = "mle") 

summary(fit_steps_exp)
#summary(fit_steps_gam)
summary(fit_steps_wei)
summary(fit_steps_lnorm)
#  AIC comparison indicates weibull distribution is best fit for observed data.

plotdist(all.steps)
qqcomp(list(fit_steps_exp, fit_steps_wei, fit_steps_lnorm), 
				legendtext=c("Exponential", "Weibull","LogNormal"),
				main="Step Lengths", xlab="Theo.",
				ylab="step lengths (m)", xlim = c(0,3500))

par(mfrow = c(2,2))

#   Examine turning angles using package CircStats
#  MLE of wrapped cauchy distribution parameters for observed data
fit_turns1 <- wrpcauchy.ml(all.turns, mu0, rho0)
	#   Get parameter estimates for wrapped Cauchy distribution using est.rho() and circ.mean()
	rho0 <- est.rho(all.turns)
	mu0<- circ.mean(all.turns)
#  MLE of Von Mises distribution parameters for observed data
fit_turns2 <- vm.ml(all.turns)	
	#   Get parameter estimates for Von Mises distribution using est.kappa() 
	kappa0 <- est.kappa(all.turns)
	
fit_turns1
fit_turns2

circ.disp(all.steps)
watson(all.steps, alpha=0.1, dist='uniform')

watson(all.steps, alpha=0.1, dist='vm')
pp.plot (all.steps, ref.line = TRUE)

v0.test(all.steps, mu0=0, degree=FALSE)

par(mfrow = c(2,2))
plot.edf(rwrpnorm(300, mu0, rho0))
plot.edf(rwrpcauchy(300, mu0, rho0))
plot.edf(dvm(300,mu0, kappa0))
plot.edf(all.turns)

circ.plot(all.steps)
################################################################################
	
#  Fit to distributions!	
#   Assuming step lengths are log-normal distribution (based on lowest AIC value
#   when fitting distributions above) and turning angles are a wrapped Cauchy distribution...
#   First look at histograms and density plots of both.
library(ggplot2)
steps_dens <- ggplot(data.out, aes(x=dist)) + 
									geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
									binwidth=100,
									colour="black", fill="grey") +
									geom_density(alpha=.2, fill="#FF6666") +
									theme_bw()
steps_dens

turns_dens <- ggplot(data.out, aes(x=turn.angle)) + 
									geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
									binwidth=0.25,
									colour="black", fill="grey") +
									geom_density(alpha=.2, fill="#FF6666") +
									xlim(-4, 4) +
									theme_bw()
turns_dens

#   Simulate random points from these distributions using parameters from observed data.
n <- 1000000

#  Step lengths.
sim_steps <- rlnorm(n, meanlog = 5.600105, sdlog =  1.171822)
dens_func_steps <- density(sim_steps)
dens_func_steps$y
sim_steps <- as.data.frame(sim_steps)
names(sim_steps)[1] <- "dist"
sim_steps <- filter(sim_steps, dist <= 5000)

sim_steps_dens <- ggplot(sim_steps, aes(x=dist)) + 
									geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
									binwidth=100,
									colour="black", fill="grey") +
									geom_density(alpha=.2, fill="#FF6666") +
									theme_bw()
sim_steps_dens

#   Turning angles.
#   Wrapped cauchy random number generation
sim_turns <- rwrpcauchy(n, mu0, rho0)
sim_turns[sim_turns>pi]=sim_turns[sim_turns>pi]-2*pi
dens_func_turns <- density(sim_turns)
dens_func_turns$y
sim_turns <- as.data.frame(sim_turns)
names(sim_turns)[1] <- "turn.angle"

sim_turns_dens <- ggplot(sim_turns, aes(x=turn.angle)) + 
									geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
									binwidth=0.25,
									colour="black", fill="grey") +
									geom_density(alpha=.2, fill="#FF6666") +
									xlim(-4, 4) +
									theme_bw()
sim_turns_dens

library(gridExtra)
grid.arrange(turns_dens, sim_turns_dens, steps_dens, sim_steps_dens, ncol=2, nrow=2)




################################################################################
#  Generate CRW, from:
#   http://wiki.cbr.washington.edu/qerm/index.php/R/Correlated_Random_Walk
library(circular)

# length of walk
  N<-10

# make weibull distributed steps
  steps <- rlnorm(N, meanlog = 5.600105, sdlog =  1.171822)
# check out their distribution
  hist(steps)

# make clustered turning angles
  theta <- rwrappedcauchy(N, .11, rho0) # transit
    theta <- rwrappedcauchy(N, .26, rho0) # station
# check out their distribution
  rose.diag(theta,bins=50)

# cumulative angle (absolute orientation)
  Phi <- cumsum(theta)

# step length components
  dX <- steps*cos(Phi)
  dY <- steps*sin(Phi)

# actual X-Y values
  X<-cumsum(dX)
  Y<-cumsum(dY)

# plot that puppy
  plot(X,Y,type="l")
	
	
for (i in 2:5){
	x[i] <- x[i-1] + steps*sin(turns)
	y[i] <- y[i-1] + steps*cos(turns)
	}
	
################################################################################

# #   Turning angles simulation option.
# #   Function from: https://github.com/benaug/move.HMM/blob/master/R/rwrpcauchy.R
# #   Wrapped cauchy random number generation
# #   This function generates random number from the wrapped cauchy distribution.
# #   It is modified from the function in the CircStats package so that it can
# #    evaluate multiple parameter combinations in the same call.
# #   @param n The number of random numbers to generate
# #   @param mu A value for the mu parameter
# #   @param rho A value for the concentration parameter in the interval [0,1]
# #   @return A vector of random numbers drawn from a wrapped cauchy distribution
# #   @export

# rwrpcauchy=function(n, mu = 0, rho = exp(-1)) 
# {
  # if(length(mu)!=length(rho))stop("Dimension of mu and rho must be equal")
  # if(length(mu)==1){
    # mu=rep(mu,n)
  # }
  # if(length(rho)==1){
    # rho=rep(rho,n)
  # }
    
  # result=rep(NA,n)
  # case1=which(rho==0)
  # case2=which(rho==1)
  # case3=1:n
  # if(length(c(case1,case2))>0){
    # case3=case3[-c(case1,case2)]
  # }
  # if(length(case1)>0){
    # result[case1]=runif(length(case1), 0, 2 * pi)
  # }
  # if(length(case2)>0){
    # result[case2]= rep(mu, length(case2))
  # }
  # if(length(case3)>0){
    # scale <- -log(rho[case3])
    # result[case3] <- rcauchy(length(case3), mu, scale)%%(2 * pi)
  # }
  # #shift support to (-pi,pi)
  # result[result>pi]=result[result>pi]-2*pi
  
  # result  
# }