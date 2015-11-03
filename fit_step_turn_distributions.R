#  Sara Williams
#  9/3/2015
#  Invesitgate step length and turning angle distributions, generate simulated data from
#   observered data distribution parameters, and practice CRW.
################################################################################

library(vcd)
library(MASS)
library(fitdistrplus)
library(CircStats)
library(qualityTools)
library(reshape)
library(ggplot2)
library(plyr)
library(dplyr)
################################################################################

#   Uses output - turning angles and step lengths - from  
#   "Generate_trajectory_ADEpackage"

#  Use package MASS fitdistr() function.  -- WEIBULL?

#   Transit observations
pos.steps.transit <- steps.transit[steps.transit > 0 & steps.transit <= quantile(steps.transit, 0.99)] 
pos.steps.transit <- sort(pos.steps.transit)
fit_steps_transit_exp <- fitdistr(pos.steps.transit, "exponential") 
fit_steps_transit_gam <- fitdistr(pos.steps.transit, "gamma")
fit_steps_transit_wei <- fitdistr(pos.steps.transit, "weibull", lower=0.001) 
fit_steps_transit_lnorm<- fitdist(pos.steps.transit, "lognormal") 


fit_steps_transit_exp$estimate
ks.test(pos.steps.transit, "pexp", fit_steps_transit_exp$estimate) # p-value < 0.05 -> distribution refused
fit_steps_transit_gam$estimate
ks.test(pos.steps.transit, "pgamma", fit_steps_transit_exp$estimate) # p-value < 0.05 -> distribution refused
fit_steps_transit_wei$estimate
ks.test(pos.steps.transit, "pweibull", fit_steps_transit_exp$estimate) # p-value < 0.05 -> distribution refused
fit_steps_transit_lnorm$estimate
ks.test(pos.steps.transit, "plnorm", fit_steps_transit_exp$estimate) # p-value < 0.05 -> distribution refused
ks.test(pos.steps.transit, pos.steps.station)

#  Use plots for diagnostics
plotdist(pos.steps.transit)
plot(ecdf(pos.steps.transit))
par(mfrow = c(2,2))
qqPlot(pos.steps.transit, "exponential", DB = TRUE, confbounds = FALSE)
qqPlot(pos.steps.transit, "log-normal", confbounds = FALSE)
qqPlot(pos.steps.transit, "weibull", confbounds = FALSE, lower=0.001)
qqPlot(pos.steps.transit, "gamma", confbounds = FALSE)

hist(pos.steps.transit, freq = FALSE, breaks = 100, xlim = c(0, quantile(pos.steps.transit, 0.99)))
curve(dweibull(x, fit_steps_transit_exp$estimate), col = "red", add = TRUE)
curve(dgamma(x, e), col = "red", add = TRUE)

#   Stationary observations
pos.steps.station <- steps.station[steps.station > 0 & steps.station <= quantile(steps.station, 0.99)]
pos.steps.station <- sort(pos.steps.station)
fit_steps_station_exp<- fitdistr(pos.steps.station, "exponential") 
fit_steps_station_gam <- fitdistr(pos.steps.station, "gamma") 
fit_steps_station_wei <- fitdistr(pos.steps.station, "weibull", lower=0.001) 
fit_steps_station_lnorm<- fitdistr(pos.steps.station, "lognormal") 

fit_steps_station_exp
fit_steps_station_gam
fit_steps_station_wei
fit_steps_station_lnorm

#  Use plots for diagnostics
plotdist(pos.steps.station)
par(mfrow = c(2,2))
qqPlot(pos.steps.station, "exponential",  DB = TRUE, confbounds = FALSE)
qqPlot(pos.steps.station, "log-normal",  confbounds = FALSE)
qqPlot(pos.steps.station, "weibull", confbounds = FALSE)
qqPlot(pos.steps.station, "gamma", confbounds = FALSE)

#   Dive observations
pos.steps.dive <- steps.dive[steps.dive > 0 & steps.dive <=  quantile(steps.dive, 0.99)]
pos.steps.dive <- sort(pos.steps.dive)
fit_steps_dive_exp<- fitdistr(pos.steps.dive, "exponential") 
fit_steps_dive_gam <- fitdistr(pos.steps.dive, "gamma") ## lowest AIC
fit_steps_dive_wei <- fitdistr(pos.steps.dive, "weibull",  lower=0.001) 
fit_steps_dive_lnorm<- fitdistr(pos.steps.dive, "lognormal") 

fit_steps_dive_exp
fit_steps_dive_gam
fit_steps_dive_wei
fit_steps_dive_lnorm

#  Use plots for diagnostics
plotdist(pos.steps.dive)
par(mfrow = c(2,2))
qqPlot(pos.steps.dive, "exponential", DB = TRUE, confbounds = FALSE, xlim = c(0, quantile(steps.dive, 0.99)))
qqPlot(pos.steps.dive, "log-normal", confbounds = FALSE)
qqPlot(pos.steps.dive, "weibull", confbounds = FALSE)
qqPlot(pos.steps.dive, "gamma", confbounds = FALSE)

#   Surface observations
pos.steps.surf <- steps.surf[steps.surf > 0 & quantile(steps.surf, 0.99)]
pos.steps.surf <- sort(pos.steps.surf)
fit_steps_surf_exp<- fitdistr(pos.steps.surf, "exponential") 
fit_steps_surf_gam <- fitdistr(pos.steps.surf, "gamma") 
fit_steps_surf_wei <- fitdistr(pos.steps.surf, "weibull", lower=0.001) 
fit_steps_surf_lnorm<- fitdistr(pos.steps.surf, "lognormal") 

AIC(fit_steps_surf_exp)
AIC(fit_steps_surf_gam)
AIC(fit_steps_surf_wei)
AIC(fit_steps_surf_lnorm)


#  Use plots for diagnostics
plotdist(pos.steps.surf)
par(mfrow = c(2,2))
qqPlot(pos.steps.surf, "exponential", xlim = c(0,8000), ylim =c(0,8000))
qqPlot(pos.steps.surf, "log-normal", xlim = c(0,8000), ylim =c(0,8000))
qqPlot(pos.steps.surf, "weibull", xlim = c(0,8000), ylim =c(0,8000))
qqPlot(pos.steps.surf, "gamma", xlim = c(0,8000), ylim =c(0,8000))
################################################################################

#   Examine turning angles using package CircStats

#   Transit observations
#   Get parameter estimates for wrapped Cauchy distribution using est.rho() and circ.mean()
rhoT <- est.rho(turns.transit)
muT<- circ.mean(turns.transit)
#  MLE of wrapped cauchy distribution parameters for observed data
fit_turns_transit_wrcauch <- wrpcauchy.ml(turns.transit, muT , rhoT)
# # # #   Get parameter estimates for Von Mises distribution using est.kappa() 
# # # kappaT <- est.kappa(turns.transit)
# # # #  MLE of Von Mises distribution parameters for observed data
# # # fit_turns_transit_vm <- vm.ml(turns.transit)	
	
fit_turns_transit_wrcauch
# # # fit_turns_transit_vm

# # # watson(turns.transit, alpha=0.05, dist='vm')

circ.disp(turns.transit)
pp.plot (turns.transit, ref.line = TRUE)

#   Stationary observations
#   Get parameter estimates for wrapped Cauchy distribution using est.rho() and circ.mean()
rhoS <- est.rho(turns.station)
muS<- circ.mean(turns.station)
#  MLE of wrapped cauchy distribution parameters for observed data
fit_turns_station_wrcauch <- wrpcauchy.ml(turns.station, muS, rhoS)
# # # #   Get parameter estimates for Von Mises distribution using est.kappa() 
# # # kappaS <- est.kappa(turns.station)
# # # #  MLE of Von Mises distribution parameters for observed data
# # # fit_turns_station_vm <- vm.ml(turns.station)	

fit_turns_station_wrcauch
# # # fit_turns_station_vm

# # # watson(turns.station, alpha=0.05, dist='vm')

circ.disp(turns.station)
pp.plot (turns.station, ref.line = TRUE)

#   Dive observations
#   Get parameter estimates for wrapped Cauchy distribution using est.rho() and circ.mean()
rhoD <- est.rho(turns.dive)
muD<- circ.mean(turns.dive)
#  MLE of wrapped cauchy distribution parameters for observed data
fit_turns_dive_wrcauch <- wrpcauchy.ml(turns.dive, muD, rhoD)
# # # #   Get parameter estimates for Von Mises distribution using est.kappa() 
# # # kappaD <- est.kappa(turns.dive)
# # # #  MLE of Von Mises distribution parameters for observed data
# # # fit_turns_dive_vm <- vm.ml(turns.dive)	

fit_turns_dive_wrcauch
# # # fit_turns_dive_vm

# # # watson(turns.dive, alpha=0.05, dist='vm')

circ.disp(turns.dive)
pp.plot (turns.dive, ref.line = TRUE)

#   Surface observations
#   Get parameter estimates for wrapped Cauchy distribution using est.rho() and circ.mean()
rhoSurf<- est.rho(turns.surf)
muSurf<- circ.mean(turns.surf)
#  MLE of wrapped cauchy distribution parameters for observed data
fit_turns_surf_wrcauch <- wrpcauchy.ml(turns.surf, muSurf, rhoSurf)
# # # #   Get parameter estimates for Von Mises distribution using est.kappa() 
# # # kappasurf <- est.kappa(turns.surf)
# # # #  MLE of Von Mises distribution parameters for observed data
# # # fit_turns_surf_vm <- vm.ml(turns.surf)	

fit_turns_surf_wrcauch
# # # fit_turns_surf_vm

# # # watson(turns.surf, alpha=0.05, dist='vm')

circ.disp(turns.surf)
pp.plot (turns.surf, ref.line = TRUE)

################################################################################

#  Other diagnostic distribution tests? 
watson(steps.dive, alpha=0.05, dist='uniform')
watson(steps, alpha=0.05, dist='vm')
# par(mfrow = c(2,1))
# #plot.edf(rwrpnorm(300, mu0, rho0))
# plot.edf(rwrpcauchy(300, mu0, rho0))
# #plot.edf(dvm(300,mu0, kappa0))
# plot.edf(turns)
################################################################################
	
# # #  Fit to distributions!	
# # #   Assuming step lengths are weibull distribution (based on lowest AIC value
# # #   when fitting distributions above) and turning angles are a wrapped Cauchy distribution...
# # #   First look at histograms and density plots of both.
# # library(ggplot2)
# # steps_dens <- ggplot(data.out, aes(x=dist)) + 
									# # geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
									# # binwidth=100,
									# # colour="black", fill="grey") +
									# # geom_density(alpha=.2, fill="#FF6666") +
									# # theme_bw()
# # steps_dens

# # turns_dens <- ggplot(data.out, aes(x=turn.angle)) + 
									# # geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
									# # binwidth=0.25,
									# # colour="black", fill="grey") +
									# # geom_density(alpha=.2, fill="#FF6666") +
									# # xlim(-4, 4) +
									# # theme_bw()
# # turns_dens
################################################################################

#####  SIMULATION  ##############################################################

#   Simulate random points from these distributions using parameters from observed data.
n <- 1000000

#  Step lengths.
#   Get shape and scale parameter estimates:
#   For transit: summary(fit_steps_transit_gam)
#   For stationary: summary(fit_steps_station_wei)
#  Transit:
sim_steps_transit <- rgamma(n, coef(fit_steps_transit_gam)[1], coef(fit_steps_transit_gam)[2])
dens_func_steps_t <- density(sim_steps_transit)
dens_func_steps_t$y
sim_steps_transit <- as.data.frame(sim_steps_transit)
names(sim_steps_transit)[1] <- "dist"
sim_steps_transit <- filter(sim_steps_transit, dist <= 5000)
sim_steps_transit <- dplyr::sample_n(sim_steps_transit, 50000, replace = TRUE)

sim_steps_dens_transit <- ggplot(sim_steps_transit, aes(x=dist)) + 
														geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
														color="black", fill="grey", 
														binwidth=250) +
														geom_density(alpha=.2, fill="#FF6666") +
														xlim(0, 4000) +
														ylim(0,0.003) +
														theme_bw()
sim_steps_dens_transit

#  Dive:
sim_steps_dive <- rgamma(n, coef(fit_steps_dive_gam)[1], coef(fit_steps_dive_gam)[2])
dens_func_steps_d <- density(sim_steps_dive)
dens_func_steps_d$y
sim_steps_dive <- as.data.frame(sim_steps_dive)
names(sim_steps_dive)[1] <- "dist"
sim_steps_dive <- filter(sim_steps_dive, dist <= 5000)
sim_steps_dive<- dplyr::sample_n(sim_steps_dive, 50000, replace = TRUE)

sim_steps_dens_dive <- ggplot(sim_steps_dive, aes(x=dist)) + 
														geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
														color="black", fill="grey", 
														binwidth=250) +
														geom_density(alpha=.2, fill="#FF6666") +
														xlim(0, 4000) +
														ylim(0,0.003) +
														theme_bw()
sim_steps_dens_dive

#   Stationary:
sim_steps_station <- rexp(n, rate = coef(fit_steps_station_exp)[1])
dens_func_steps_s <- density(sim_steps_station)
dens_func_steps_s$y
sim_steps_station <- as.data.frame(sim_steps_station)
names(sim_steps_station)[1] <- "dist"
sim_steps_station <- filter(sim_steps_station, dist <= 5000)
sim_steps_station <- dplyr::sample_n(sim_steps_station, 50000, replace = TRUE)

sim_steps_dens_station<- ggplot(sim_steps_station, aes(x=dist)) + 
														geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
														color="black", fill="grey", 
														binwidth=250) +
														geom_density(alpha=.2, fill="#FF6666") +
														xlim(0, 4000) +
														ylim(0,0.003) +
													theme_bw()
sim_steps_dens_station

#   Turning angles via Von Mises random generation

#   For transit: 
sim_turns_transit <- rwrpcauchy(n, muT, rhoT)
sim_turns_transit[sim_turns_transit>pi]=sim_turns_transit[sim_turns_transit>pi]-2*pi
dens_func_turns_transit <- density(sim_turns_transit)
dens_func_turns_transit$y
sim_turns_transit <- as.data.frame(sim_turns_transit)
names(sim_turns_transit)[1] <- "turn.angle"
sim_turns_transit <- dplyr::sample_n(sim_turns_transit, 50000, replace = TRUE)

sim_turns_dens_transit <- ggplot(sim_turns_transit, aes(x=turn.angle)) + 
														geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
														color="black", fill="grey", 
														binwidth=.1) +
														geom_density(alpha=.2, fill="#FF6666") +
														xlim(-3.5, 3.5) +
														ylim(0,1.0) +
														theme_bw()
sim_turns_dens_transit

#   For dive: 
sim_turns_dive <- rwrpcauchy(n, muD, rhoD)
sim_turns_dive[sim_turns_dive>pi]=sim_turns_dive[sim_turns_dive>pi]-2*pi
dens_func_turns_dive <- density(sim_turns_dive)
dens_func_turns_dive$y
sim_turns_dive <- as.data.frame(sim_turns_dive)
names(sim_turns_dive)[1] <- "turn.angle"
sim_turns_dive <- dplyr::sample_n(sim_turns_dive, 50000, replace = TRUE)

sim_turns_dens_dive <- ggplot(sim_turns_dive, aes(x=turn.angle)) + 
														geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
														color="black", fill="grey", 
														binwidth=.1) +
														geom_density(alpha=.2, fill="#FF6666") +
														xlim(-3.5, 3.5) +
														ylim(0,1.0) +
														theme_bw()
sim_turns_dens_dive

#   For stationary:
sim_turns_station <- rwrpcauchy(n, muS, rhoS)
sim_turns_station[sim_turns_station>pi]=sim_turns_station[sim_turns_station>pi]-2*pi
dens_func_turns_station <- density(sim_turns_station)
dens_func_turns_station$y
sim_turns_station <- as.data.frame(sim_turns_station)
names(sim_turns_station)[1] <- "turn.angle"
sim_turns_station <- dplyr::sample_n(sim_turns_station, 50000, replace = TRUE)

sim_turns_dens_station <- ggplot(sim_turns_station, aes(x=turn.angle)) + 
														geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
														color="black", fill="grey", 
														binwidth=.1) +
														geom_density(alpha=.2, fill="#FF6666") +
														xlim(-3.5, 3.5) +
														ylim(0,1.0) +
														theme_bw()
sim_turns_dens_station
################################################################################

#  Plots
#  Create dataframe of simulated transit steps and simulated stationary steps
sim_steps <- cbind(sim_steps_dive, sim_steps_station)
names(sim_steps)[1] <- "dist_dive"
names(sim_steps)[2] <- "dist_station"
#  Look at density plots of both on same figure
df_steps_sim <- melt(sim_steps)
compare_steps_sim <- ggplot(df_steps_sim) + 
											geom_density(aes(x = value, colour = variable)) +
											xlim(0, 4000) +
											theme_bw()
 compare_steps_sim
 
 #  Create dataframe of simulated transit steps and simulated stationary steps
sim_turns <- cbind(sim_turns_dive, sim_turns_station)
names(sim_turns)[1] <- "turns_dive"
names(sim_turns)[2] <- "turns_station"
#  Look at density plots of both on same figure
df_turns_sim <- melt(sim_turns)
compare_turns_sim <- ggplot(df_turns_sim) + 
											geom_density(aes(x = value, colour = variable)) +
											xlim(-4, 4) +
											theme_bw()
 compare_turns_sim
################################################################################

#  Generate CRW, from:
#   http://wiki.cbr.washington.edu/qerm/index.php/R/Correlated_Random_Walk
library(circular)

# length of walk
  N<-10

 par(mfrow=c(1,2))
 
# make weibull distributed steps
  steps_dive <-(rweibull(N, coef(fit_steps_dive_wei)[1], coef(fit_steps_dive_wei)[2]))/1000
  steps_station <- (rweibull(N, coef(fit_steps_station_wei)[1], coef(fit_steps_station_wei)[2]))/1000
  hist(steps_dive)
  hist(steps_station)

# make clustered turning angles
  theta_dive <- rvm(N, muD, kappaD)
  theta_station <- rvm(N, muS, kappaS) 
# check out their distribution
  rose.diag(theta_dive,bins=50)
  rose.diag(theta_station,bins=50)

# cumulative angle (absolute orientation)
  phi_dive <- cumsum(theta_dive)
  phi_station <- cumsum(theta_station)

# step length components
  dX_dive <- steps_dive*cos(phi_dive)
  dY_dive <- steps_dive*sin(phi_dive)
  dX_station <- steps_station*cos(phi_station)
  dY_station <- steps_station*sin(phi_station)

# actual X-Y values
  X_dive<-cumsum(dX_dive)
  Y_dive<-cumsum(dY_dive)
  X_station<-cumsum(dX_station)
  Y_station<-cumsum(dY_station)
  
# plot that puppy
   plot(X_dive,Y_dive,type="l", xlim=c(-5, 5), ylim = c(-5, 5))
   plot(X_station,Y_station,type="l",  xlim=c(-5,5),  ylim = c(-5, 5))	
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