#  Sara Williams
#  9/3/2015; updated 12/4/2015
#  Invesitgate step length and turning angle distributions, generate simulated data from
#   observered data distribution parameters, and practice CRW.
#   Uses output - turning angles and step lengths - from "Generate_trajectory_ADEpackage"
#   script.
################################################################################

#  Load packages
library(vcd)
library(MASS)
library(fitdistrplus)
library(CircStats)
library(qualityTools)
library(reshape)
library(ggplot2)
library(plyr)
library(dplyr)

#  Use package MASS fitdistr() function.  -- WEIBULL?

# #   Transit observations
# pos.steps.transit <- steps.transit[steps.transit > 0 & steps.transit <= 5000] 
# pos.steps.transit <- sort(pos.steps.transit)
# fit_steps_transit_exp <- fitdistr(pos.steps.transit, "exponential") 
# fit_steps_transit_gam <- fitdistr(pos.steps.transit, "gamma")
# fit_steps_transit_wei <- fitdistr(pos.steps.transit, "weibull", lower=0.001) 
# fit_steps_transit_lnorm<- fitdist(pos.steps.transit, "lognormal") 


# fit_steps_transit_exp$estimate
# ks.test(pos.steps.transit, "pexp", fit_steps_transit_exp$estimate) # p-value < 0.05 -> distribution refused
# fit_steps_transit_gam$estimate
# ks.test(pos.steps.transit, "pgamma", fit_steps_transit_exp$estimate) # p-value < 0.05 -> distribution refused
# fit_steps_transit_wei$estimate
# ks.test(pos.steps.transit, "pweibull", fit_steps_transit_exp$estimate) # p-value < 0.05 -> distribution refused
# fit_steps_transit_lnorm$estimate
# ks.test(pos.steps.transit, "plnorm", fit_steps_transit_exp$estimate) # p-value < 0.05 -> distribution refused
# ks.test(pos.steps.transit, pos.steps.station)

# #  Use plots for diagnostics
# plotdist(pos.steps.transit)
# #plot(ecdf(pos.steps.transit))
# par(mfrow = c(2,2))
# qqPlot(pos.steps.transit, "exponential", DB = TRUE, confbounds = FALSE, xlim=c(0,5000))
# qqPlot(pos.steps.transit, "log-normal", confbounds = FALSE, xlim=c(0,5000))
# qqPlot(pos.steps.transit, "weibull", confbounds = FALSE, xlim=c(0,5000))
# qqPlot(pos.steps.transit, "gamma", confbounds = FALSE, xlim=c(0,5000))

# hist(pos.steps.transit, freq = FALSE, breaks = 100, xlim = c(0, quantile(pos.steps.transit, 0.99)))
# curve(dweibull(x, fit_steps_transit_exp$estimate), col = "red", add = TRUE)
# curve(dgamma(x, e), col = "red", add = TRUE)

# #   Stationary observations --- Weibull?
# pos.steps.station <- steps.station[steps.station > 0 & steps.station <= 5000]
# pos.steps.station <- sort(pos.steps.station)
# fit_steps_station_exp<- fitdistr(pos.steps.station, "exponential") 
# fit_steps_station_gam <- fitdistr(pos.steps.station, "gamma") 
# fit_steps_station_wei <- fitdistr(pos.steps.station, "weibull", lower=0.001) 
# fit_steps_station_lnorm<- fitdistr(pos.steps.station, "lognormal") 

# fit_steps_station_exp
# fit_steps_station_gam
# fit_steps_station_wei
# fit_steps_station_lnorm

# #  Use plots for diagnostics
# plotdist(pos.steps.station)
# par(mfrow = c(2,2))
# qqPlot(pos.steps.station, "exponential",  DB = TRUE, confbounds = FALSE, xlim=c(0,5000))
# qqPlot(pos.steps.station, "log-normal",  confbounds = FALSE, xlim=c(0,5000))
# qqPlot(pos.steps.station, "weibull", confbounds = FALSE, xlim=c(0,5000))
# qqPlot(pos.steps.station, "gamma", confbounds = FALSE, xlim=c(0,5000))

#   Dive observations --- Exponential?
pos.steps.dive <- steps.dive[steps.dive > 0 & steps.dive <=  5000]
pos.steps.dive <- sort(pos.steps.dive)
fit_steps_dive_exp<- fitdistr(pos.steps.dive, "exponential") 
fit_steps_dive_gam <- fitdistr(pos.steps.dive, "gamma") ## lowest AIC
fit_steps_dive_wei <- fitdistr(pos.steps.dive, "weibull",  lower=0.001) 
fit_steps_dive_lnorm<- fitdistr(pos.steps.dive, "lognormal") 

fit_steps_dive_exp
fit_steps_dive_gam
fit_steps_dive_wei
fit_steps_dive_lnorm

AIC(fit_steps_dive_exp)
AIC(fit_steps_dive_gam)
AIC(fit_steps_dive_wei)
AIC(fit_steps_dive_lnorm)

#  Use plots for diagnostics 
#plotdist(pos.steps.dive)
par(mfrow = c(2,2))
qqPlot(pos.steps.dive, "exponential", DB = TRUE, confbounds = FALSE,  xlim=c(0,5000))
qqPlot(pos.steps.dive, "log-normal", confbounds = FALSE, xlim=c(0,5000))
qqPlot(pos.steps.dive, "weibull", confbounds = FALSE, xlim=c(0,5000))
qqPlot(pos.steps.dive, "gamma", confbounds = FALSE, xlim=c(0,5000))

#   Surface observations --- Weibull?
pos.steps.surf <- steps.surf[steps.surf > 0 & steps.surf <=  5000]
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
#plotdist(pos.steps.surf)
par(mfrow = c(2,2))
qqPlot(pos.steps.surf, "exponential", confbounds = FALSE, xlim = c(0,5000))
qqPlot(pos.steps.surf, "log-normal", confbounds = FALSE, xlim = c(0,5000))
qqPlot(pos.steps.surf, "weibull", confbounds = FALSE, xlim = c(0,5000))
qqPlot(pos.steps.surf, "gamma", confbounds = FALSE, xlim = c(0,5000))
################################################################################

#   Examine turning angles using package CircStats

# #   Transit observations
# #   Get parameter estimates for wrapped Cauchy distribution using est.rho() and circ.mean()
# rhoT <- est.rho(turns.transit)
# muT<- circ.mean(turns.transit)
# #  MLE of wrapped cauchy distribution parameters for observed data
# fit_turns_transit_wrcauch <- wrpcauchy.ml(turns.transit, muT , rhoT)
# # # # #   Get parameter estimates for Von Mises distribution using est.kappa() 
# # # # kappaT <- est.kappa(turns.transit)
# # # # #  MLE of Von Mises distribution parameters for observed data
# # # # fit_turns_transit_vm <- vm.ml(turns.transit)	
	
# fit_turns_transit_wrcauch
# # # # fit_turns_transit_vm

# # # # watson(turns.transit, alpha=0.05, dist='vm')

# circ.disp(turns.transit)
# pp.plot (turns.transit, ref.line = TRUE)

# #   Stationary observations
# #   Get parameter estimates for wrapped Cauchy distribution using est.rho() and circ.mean()
# rhoS <- est.rho(turns.station)
# muS<- circ.mean(turns.station)
# #  MLE of wrapped cauchy distribution parameters for observed data
# fit_turns_station_wrcauch <- wrpcauchy.ml(turns.station, muS, rhoS)
# # # # #   Get parameter estimates for Von Mises distribution using est.kappa() 
# # # # kappaS <- est.kappa(turns.station)
# # # # #  MLE of Von Mises distribution parameters for observed data
# # # # fit_turns_station_vm <- vm.ml(turns.station)	

# fit_turns_station_wrcauch
# # # # fit_turns_station_vm

# # # # watson(turns.station, alpha=0.05, dist='vm')

# circ.disp(turns.station)
# pp.plot (turns.station, ref.line = TRUE)

#   Dive observations
#   Get parameter estimates for wrapped Cauchy distribution using est.rho() and circ.mean()
rhoD <- est.rho(turns.dive)
muD<- circ.mean(turns.dive)
#  MLE of wrapped cauchy distribution parameters for observed data
fit_turns_dive_wrcauch <- wrpcauchy.ml(turns.dive, muD, rhoD)
fit_turns_dive_wrcauch

# # # #   Get parameter estimates for Von Mises distribution using est.kappa() 
# # # kappaD <- est.kappa(turns.dive)
# # # #  MLE of Von Mises distribution parameters for observed data
# # # fit_turns_dive_vm <- vm.ml(turns.dive)	
# # # fit_turns_dive_vm
# # # watson(turns.dive, alpha=0.05, dist='vm')
# # # circ.disp(turns.dive)
# # # pp.plot (turns.dive, ref.line = TRUE)

#   Surface observations
#   Get parameter estimates for wrapped Cauchy distribution using est.rho() and circ.mean()
rhoSurf<- est.rho(turns.surf)
muSurf<- circ.mean(turns.surf)
#  MLE of wrapped cauchy distribution parameters for observed data
fit_turns_surf_wrcauch <- wrpcauchy.ml(turns.surf, muSurf, rhoSurf)
fit_turns_surf_wrcauch

# # # #   Get parameter estimates for Von Mises distribution using est.kappa() 
# # # kappasurf <- est.kappa(turns.surf)
# # # #  MLE of Von Mises distribution parameters for observed data
# # # fit_turns_surf_vm <- vm.ml(turns.surf)	
# # # fit_turns_surf_vm
# # # watson(turns.surf, alpha=0.05, dist='vm')
# # # circ.disp(turns.surf)
# # # pp.plot (turns.surf, ref.line = TRUE)

################################################################################

#  Other diagnostic distribution tests? 
# watson(steps.dive, alpha=0.05, dist='uniform')
# watson(steps, alpha=0.05, dist='vm')
# par(mfrow = c(2,1))
# #plot.edf(rwrpnorm(300, mu0, rho0))
# plot.edf(rwrpcauchy(300, mu0, rho0))
# #plot.edf(dvm(300,mu0, kappa0))
# plot.edf(turns)
################################################################################
	
#  Fit to distributions!	
#   Assuming  dive steps are exponential distribution and surface steps are weibull distribution
#    (based on visual inspection of Q-Q plots) and turning angles are a wrapped Cauchy distribution.

# #   First look at histograms and density plots of both.
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
# sim_steps_transit <- rgamma(n, coef(fit_steps_transit_gam)[1], coef(fit_steps_transit_gam)[2])
# dens_func_steps_t <- density(sim_steps_transit)
# dens_func_steps_t$y
# sim_steps_transit <- as.data.frame(sim_steps_transit)
# names(sim_steps_transit)[1] <- "dist"
# sim_steps_transit <- filter(sim_steps_transit, dist <= 5000)
# sim_steps_transit <- dplyr::sample_n(sim_steps_transit, 50000, replace = TRUE)

# sim_steps_dens_transit <- ggplot(sim_steps_transit, aes(x=dist)) + 
														# geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
														# color="black", fill="grey", 
														# binwidth=250) +
														# geom_density(alpha=.2, fill="#FF6666") +
														# xlim(0, 4000) +
														# ylim(0,0.003) +
														# theme_bw()
# sim_steps_dens_transit

#  Dive:
sim_steps_dive <- rexp(n, coef(fit_steps_dive_exp)[1])
dens_func_steps_d <- density(sim_steps_dive)
dens_func_steps_d$y
sim_steps_dive <- as.data.frame(sim_steps_dive)
names(sim_steps_dive)[1] <- "dist"
sim_steps_dive <- filter(sim_steps_dive, dist <= 5000)
sim_steps_dive<- dplyr::sample_n(sim_steps_dive, 50000, replace = TRUE)

sim_steps_dens_dive <- ggplot(sim_steps_dive, aes(x=dist)) + 
														geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
														color="black", fill="grey", 
														binwidth=1) +
														#geom_density(alpha=.2, fill="#FF6666") +
														xlim(0, 5000) +
														ylim(0,0.003) +
														xlab("Step length (m)") +
														theme_bw()
sim_steps_dens_dive

#   Surface:
sim_steps_surf <- rweibull(n, coef(fit_steps_surf_wei)[1], coef(fit_steps_surf_wei)[2])
dens_func_steps_s <- density(sim_steps_surf)
dens_func_steps_s$y
sim_steps_surf <- as.data.frame(sim_steps_surf)
names(sim_steps_surf)[1] <- "dist"
sim_steps_surf <- filter(sim_steps_surf, dist <= 5000)
sim_steps_surf <- dplyr::sample_n(sim_steps_surf, 50000, replace = TRUE)

sim_steps_dens_surf <- ggplot(sim_steps_surf, aes(x=dist)) + 
												 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
												 color="black", fill="grey", 
												 binwidth=1) +
												 #geom_density(alpha=.2, fill="#FF6666") +
												 xlim(0, 5000) +
												 ylim(0,0.003) +
												 xlab("Step length (m)") +
												theme_bw()
sim_steps_dens_surf

#   Turning angles via Wrapped Cauchy random generation

# #   For transit: 
# sim_turns_transit <- rwrpcauchy(n, muT, rhoT)
# sim_turns_transit[sim_turns_transit>pi]=sim_turns_transit[sim_turns_transit>pi]-2*pi
# dens_func_turns_transit <- density(sim_turns_transit)
# dens_func_turns_transit$y
# sim_turns_transit <- as.data.frame(sim_turns_transit)
# names(sim_turns_transit)[1] <- "turn.angle"
# sim_turns_transit <- dplyr::sample_n(sim_turns_transit, 50000, replace = TRUE)

# sim_turns_dens_transit <- ggplot(sim_turns_transit, aes(x=turn.angle)) + 
														# geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
														# color="black", fill="grey", 
														# binwidth=.1) +
														# geom_density(alpha=.2, fill="#FF6666") +
														# xlim(-3.5, 3.5) +
														# ylim(0,1.0) +
														# theme_bw()
# sim_turns_dens_transit

#  Dive: 
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
												  binwidth=.001) +
												  #geom_density(alpha=.2, fill="#FF6666") +
												  xlim(-3.5, 3.5) +
												  ylim(0,1.0) +
												  xlab("Turn angle (radian)") +
												  theme_bw()
sim_turns_dens_dive

#  Surf:
sim_turns_surf <- rwrpcauchy(n, muSurf, rhoSurf)
sim_turns_surf[sim_turns_surf>pi]=sim_turns_surf[sim_turns_surf>pi]-2*pi
dens_func_turns_surf <- density(sim_turns_surf)
dens_func_turns_surf$y
sim_turns_surf <- as.data.frame(sim_turns_surf)
names(sim_turns_surf)[1] <- "turn.angle"
sim_turns_surf <- dplyr::sample_n(sim_turns_surf, 50000, replace = TRUE)

sim_turns_dens_surf <- ggplot(sim_turns_surf, aes(x=turn.angle)) + 
												 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
												 color="black", fill="grey", 
												 binwidth=.001) +
												 #geom_density(alpha=.2, fill="#FF6666") +
												 xlim(-3.5, 3.5) +
												 ylim(0,1.0) +
												 xlab("Turn angle (radian)") +
												 theme_bw()
sim_turns_dens_surf
################################################################################

#  Plots
#  Create dataframe of simulated dive steps and simulated surface steps
sim_steps <- cbind(sim_steps_dive, sim_steps_surf)
names(sim_steps)[1] <- "Dive"
names(sim_steps)[2] <- "Surf"
#  Look at density plots of both on same figure
df_steps_sim <- melt(sim_steps)
compare_steps_sim <- ggplot(df_steps_sim) + 
											   geom_density(aes(x = value, colour = variable, fill = variable), alpha = .25) +
											   xlim(0, 4000) +
											   xlab("Step length (m)") +
											   theme_bw()
 compare_steps_sim
  
#  Compute shared area under curves???
#   Generate kernel densities
ddive <- density(sim_steps$Dive, from=0, to=4000)
dsurf <- density(sim_steps$Surf, from=0, to=4000)
d <- data.frame(x=ddive$x, a=ddive$y, b=dsurf$y)
#   Calculate intersection densities
d$w <- pmin(d$a, d$b)
#   Integrate areas under curves
library(sfsmisc)
total <- integrate.xy(d$x, d$a) + integrate.xy(d$x, d$b)
intersection <- integrate.xy(d$x, d$w)
#   Compute overlap coefficient
overlap <- 2 * intersection / total
 
#  Create dataframe of simulated dive turns and simulated surface turns
sim_turns <- cbind(sim_turns_dive, sim_turns_surf)
names(sim_turns)[1] <- "Dive"
names(sim_turns)[2] <- "Surf"
#  Look at density plots of both on same figure
df_turns_sim <- melt(sim_turns)
compare_turns_sim <- ggplot(df_turns_sim) + 
											   geom_density(aes(x = value, colour = variable, fill = variable), alpha = .25) +
											   xlim(-4, 4) +
											   xlab("Turn angle (radian)") +
											   theme_bw()
 compare_turns_sim
 
#  Compute shared area under curves???
#   Generate kernel densities
ddive_s <- density(sim_turns$Dive, from=-4, to=4)
dsurf_s <- density(sim_turns$Surf, from=-4, to=4)
d_s <- data.frame(x=ddive_s$x, a=ddive_s$y, b=dsurf_s$y)
#   Calculate intersection densities
d_s$w <- pmin(d_s$a, d_s$b)
#   Integrate areas under curves
library(sfsmisc)
total_s <- integrate.xy(d_s$x, d_s$a) + integrate.xy(d_s$x, d_s$b)
intersection_s <- integrate.xy(d_s$x, d_s$w)
#   Compute overlap coefficient
overlap_s <- 2 * intersection_s / total_s
################################################################################

#  Generate 1 single CRW (code from link below)
#   http://wiki.cbr.washington.edu/qerm/index.php/R/Correlated_Random_Walk
library(circular)

# length of walk
N<-5

# make distributed steps
  steps_dive <-rexp(N, coef(fit_steps_dive_exp)[1])/1000
  steps_surf <- rweibull(N, coef(fit_steps_surf_wei)[1], coef(fit_steps_surf_wei)[2])/1000
# check out their distribution  
  par(mfrow=c(1,2))
  hist(steps_dive)
  hist(steps_surf)

# make clustered turning angles
  theta_dive <- rwrpcauchy(N, muD, rhoD)
  theta_surf <- rwrpcauchy(N, muSurf, rhoSurf)
# check out their distribution
  par(mfrow=c(1,2))
  # rose.diag(theta_dive,bins=50)
  # rose.diag(theta_surf,bins=50)

# cumulative angle (absolute orientation)
  phi_dive <- cumsum(theta_dive)
  phi_surf <- cumsum(theta_surf)

# step length components
  dX_dive <- steps_dive*cos(phi_dive)
  dY_dive <- steps_dive*sin(phi_dive)
  dX_surf <- steps_surf*cos(phi_surf)
  dY_surf <- steps_surf*sin(phi_surf)

# actual X-Y values
  X_dive<-cumsum(dX_dive)
  Y_dive<-cumsum(dY_dive)
  X_surf<-cumsum(dX_surf)
  Y_surf<-cumsum(dY_surf)
  
dive_df <- as.data.frame(cbind(X_dive, Y_dive)) 
dive_df <- mutate(dive_df, group = "Dive")    
names(dive_df)[1] <- "X"
names(dive_df)[2] <- "Y"   

surf_df <- as.data.frame(cbind(X_surf, Y_surf))     
surf_df <- mutate(surf_df, group = "Surf")
names(surf_df)[1] <- "X"
names(surf_df)[2] <- "Y"   
traj_df <- as.data.frame(bind_rows(dive_df, surf_df))

compare_sim_traj <- ggplot(traj_df) + 
										   geom_line(aes(x = X, y = Y, colour = group)) +
										   xlab("X (km)") +
										   ylab("Y (km)") +
										   theme_bw()
compare_sim_traj  
################################################################################

#  Generate many simulations of N steps of movement in diving or surface mode and compare.
#  Number of simulations
sims <- 1000
#   Length of walk
N <- 4
#   Dive
dive_mat_X <- matrix(nrow=N, ncol = sims)
dive_mat_Y <- matrix(nrow=N, ncol = sims)
for(j in 1: sims){
	 # make distributed steps
	 steps_dive <-rexp(N, coef(fit_steps_dive_exp)[1])/1000
	 # make clustered turning angles
	 theta_dive <- rwrpcauchy(N, muD, rhoD)
	 # cumulative angle (absolute orientation)
	 phi_dive <- cumsum(theta_dive)
	 # step length components
	 dX_dive <- steps_dive*cos(phi_dive)
	 dY_dive <- steps_dive*sin(phi_dive)
	 # actual X-Y values
	 X_dive <- as.matrix(cumsum(dX_dive))
	 Y_dive <- as.matrix(cumsum(dY_dive))
	 dive_mat_X[,j] <- X_dive
	 dive_mat_Y[,j] <- Y_dive
}
#   Surface
surf_mat_X <- matrix(nrow=N, ncol = sims)
surf_mat_Y <- matrix(nrow=N, ncol = sims)
for(j in 1: sims){
	 # make distributed steps
	 steps_surf <- rweibull(N, coef(fit_steps_surf_wei)[1], coef(fit_steps_surf_wei)[2])/1000
	 # make clustered turning angles
	 theta_surf <- rwrpcauchy(N, muSurf, rhoSurf)
	 # cumulative angle (absolute orientation)
	 phi_surf <- cumsum(theta_surf)
	 # step length components
	 dX_surf <- steps_surf*cos(phi_surf)
	 dY_surf <- steps_surf*sin(phi_surf)
	 # actual X-Y values
	 X_surf <- as.matrix(cumsum(dX_surf))
	 Y_surf <- as.matrix(cumsum(dY_surf))
	 surf_mat_X[,j] <- X_surf
	 surf_mat_Y[,j] <- Y_surf
}

dive_df_X_tmp <- as.data.frame(dive_mat_X)
dive_df_Y_tmp <- as.data.frame(dive_mat_Y)
dive_df_X <- melt(dive_df_X_tmp)
dive_df_Y <- melt(dive_df_Y_tmp)
dive_df_XY_tmp <- cbind(dive_df_X, dive_df_Y)
names(dive_df_XY_tmp) <- c("walk_num", "X", "walk_num_rep", "Y")
dive_df_XY <- dive_df_XY_tmp %>%
							 dplyr::select(walk_num, X, Y) %>%
							 mutate(mov_mod = "Dive")

surf_df_X_tmp <- as.data.frame(surf_mat_X)
surf_df_Y_tmp <- as.data.frame(surf_mat_Y)
surf_df_X <- melt(surf_df_X_tmp)
surf_df_Y <- melt(surf_df_Y_tmp)
surf_df_XY_tmp <- cbind(surf_df_X,surf_df_Y)
names(surf_df_XY_tmp) <- c("walk_num", "X", "walk_num_rep", "Y")
surf_df_XY <- surf_df_XY_tmp %>%
							 dplyr::select(walk_num, X, Y) %>%
							 mutate(mov_mod = "Surf")

sim_crw_df <- as.data.frame(bind_rows(dive_df_XY , surf_df_XY))
#testlengths$replane <- paste(testlengths$replicate, testlengths$lane, sep="_")

compare_sim_crw <- ggplot(sim_crw_df) + 
										   geom_line(aes(x = X, y = Y, 
																		   group = walk_num, 
																		   colour = mov_mod), alpha = 0.3) +
										   xlab("X (km)") +
										   ylab("Y (km)") +
										   theme_bw()
compare_sim_crw 
################################################################################

