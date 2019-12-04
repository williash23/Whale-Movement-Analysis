# Sara Williams
# 2/22/2016; updated 10/10/2019
# Generate figures from parameter estimates from movement state models
#   and run in JAGS. For data seperated into surface interval and deep dive by time between
#   observations.
################################################################################

#  Load packages
#library(CircStats)
#library(reshape)
library(ggplot2)
library(dplyr)
library(rjags)



################################################################################
# Load output from MCMC posterior distribtuion iterations (JAGS object)
load("C:/Users/saraw/OneDrive/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/out_rand_int_dist.Rdata")
load("C:/Users/saraw/OneDrive/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/out_rand_int_bear.Rdata")
################################################################################

# NOTE: due to differences in formulation of Weibull distribution between JAGS and R:
###### scale = (1/lambda)^(1/v) ######

rwcauchy <- function(n, mu = 0, rho = 0) {
  u = runif(n)
  V = cos(2 * pi * u)
  c = 2 * rho/(1 + rho^2)
  t <- (sign(runif(n) - 0.5) * acos((V + c)/(1 + c * V)) + mu)%%(2 * pi)
  return(t)
}
################################################################################



# Parameter values - Distance
mu_dist1 <- c(out_rand_int_dist$samples[[1]][,1], out_rand_int_dist$samples[[2]][,1], 
	out_rand_int_dist$samples[[3]][,1])
mu_dist2 <- c(out_rand_int_dist$samples[[1]][,2], out_rand_int_dist$samples[[2]][,2], 
	out_rand_int_dist$samples[[3]][,2])
mu_dist3 <- c(out_rand_int_dist$samples[[1]][,3], out_rand_int_dist$samples[[2]][,3], 
	out_rand_int_dist$samples[[3]][,3])

rho_dist1 <- c(out_rand_int_dist$samples[[1]][,4], out_rand_int_dist$samples[[2]][,4], 
	out_rand_int_dist$samples[[3]][,4])
rho_dist2 <- c(out_rand_int_dist$samples[[1]][,5], out_rand_int_dist$samples[[2]][,5], 
	out_rand_int_dist$samples[[3]][,5])
rho_dist3 <- c(out_rand_int_dist$samples[[1]][,6], out_rand_int_dist$samples[[2]][,6], 
	out_rand_int_dist$samples[[3]][,6])

v_dist1 <- c(out_rand_int_dist$samples[[1]][,7], out_rand_int_dist$samples[[2]][,7], 
	out_rand_int_dist$samples[[3]][,7])
v_dist2 <- c(out_rand_int_dist$samples[[1]][,8], out_rand_int_dist$samples[[2]][,8], 
	out_rand_int_dist$samples[[3]][,8])
v_dist3 <- c(out_rand_int_dist$samples[[1]][,9], out_rand_int_dist$samples[[2]][,9], 
	out_rand_int_dist$samples[[3]][,9])

lambda_dist1 <- c(out_rand_int_dist$samples[[1]][,10], out_rand_int_dist$samples[[2]][,10], 
	out_rand_int_dist$samples[[3]][,10])
lambda_dist2 <- c(out_rand_int_dist$samples[[1]][,11], out_rand_int_dist$samples[[2]][,11], 
	out_rand_int_dist$samples[[3]][,11])
lambda_dist3 <- c(out_rand_int_dist$samples[[1]][,12], out_rand_int_dist$samples[[2]][,12], 
	out_rand_int_dist$samples[[3]][,12])
	

steps_dist1 <- numeric(length = length(mu_dist1))
turns_dist1 <- numeric(length = length(mu_dist1))
steps_dist2 <- numeric(length = length(mu_dist1))
turns_dist2 <- numeric(length = length(mu_dist1))
steps_dist3 <- numeric(length = length(mu_dist1))
turns_dist3 <- numeric(length = length(mu_dist1))

for(i in 1:length(steps_dist1)){
	steps_dist1[i] <- rweibull(1, v_dist1[i], (1/lambda_dist1[i])^(1/v_dist1[i]))
	turns_dist1[i] <- rwcauchy(1, mu_dist1[i], rho_dist1[i])

	steps_dist2[i] <- rweibull(1, v_dist2[i], (1/lambda_dist2[i])^(1/v_dist2[i]))
	turns_dist2[i] <- rwcauchy(1, mu_dist2[i], rho_dist2[i])
	
	steps_dist3[i] <- rweibull(1, v_dist3[i], (1/lambda_dist3[i])^(1/v_dist3[i]))
	turns_dist3[i] <- rwcauchy(1, mu_dist3[i], rho_dist3[i])
	}


sims_dist1 <- as.data.frame(cbind(steps_dist1, turns_dist1))
sims_dist1$turns_dist1[sims_dist1$turns_dist1 > pi] = sims_dist1$turns_dist1[sims_dist1$turns_dist1 > pi] - 2*pi

sims_dist2 <- as.data.frame(cbind(steps_dist2, turns_dist2))
sims_dist2$turns_dist2[sims_dist2$turns_dist2 > pi] = sims_dist2$turns_dist2[sims_dist2$ turns_dist2 > pi] - 2*pi

sims_dist3 <- as.data.frame(cbind(steps_dist3, turns_dist3))
sims_dist3$turns_dist3[sims_dist3$turns_dist3 > pi] = sims_dist3$turns_dist3[sims_dist3$turns_dist3 > pi] - 2*pi

keep <- sample(1:45000, 1000, replace = F)
sims_dist1 <- sims_dist1[keep,] %>%
	mutate(index = 1) %>%
	mutate(iter = 1:nrow(.))
colnames(sims_dist1) <- c("steps", "turns", "index", "iter")
sims_dist2 <- sims_dist2[keep,] %>%
	mutate(index = 2) %>%
	mutate(iter= 1:nrow(.))
colnames(sims_dist2) <- c("steps", "turns", "index", "iter")
sims_dist3 <- sims_dist3[keep,] %>%
	mutate(index = 3) %>%
	mutate(iter= 1:nrow(.))
colnames(sims_dist3) <- c("steps", "turns", "index", "iter")

sims <- sims_dist1 %>%
	rbind(sims_dist2) %>%
	rbind(sims_dist3) %>%
	mutate(turns_deg = NISTunits::NISTradianTOdeg(turns))


#  Density plots
steps_p <- ggplot(sims, aes(steps, group = factor(index))) + 
	geom_density(aes(linetype = factor(index)), size = 0.75) +
	scale_linetype_manual(name = "Distance Index", values = c(1, 6, 3)) +
	xlab("\n Step length (km)") +
	ylab("Density \n") +
	xlim(c(0, 3)) +
	ylim(c(0, 2)) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		legend.position = c(0.88, 0.88)) 
steps_p

turns_p <- ggplot(sims, aes(turns_deg, group = factor(index))) + 
	geom_density(aes(linetype = factor(index)), size = 0.75) +
	scale_linetype_manual(name = "Distance Index", values = c(1, 6, 3)) +
	xlab("\n Turn angle (deg)") +
	ylab("Density \n") +
	xlim(c(-180, 180)) +
	ylim(c(0, 0.005)) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		legend.position = c(0.88, 0.88)) 
turns_p




# Parameter values - Position
mu_bear1 <- c(out_rand_int_bear$samples[[1]][,1], out_rand_int_bear$samples[[2]][,1], 
	out_rand_int_bear$samples[[3]][,1])
mu_bear2 <- c(out_rand_int_bear$samples[[1]][,2], out_rand_int_bear$samples[[2]][,2], 
	out_rand_int_bear$samples[[3]][,2])
mu_bear3 <- c(out_rand_int_bear$samples[[1]][,3], out_rand_int_bear$samples[[2]][,3], 
	out_rand_int_bear$samples[[3]][,3])

rho_bear1 <- c(out_rand_int_bear$samples[[1]][,4], out_rand_int_bear$samples[[2]][,4], 
	out_rand_int_bear$samples[[3]][,4])
rho_bear2 <- c(out_rand_int_bear$samples[[1]][,5], out_rand_int_bear$samples[[2]][,5], 
	out_rand_int_bear$samples[[3]][,5])
rho_bear3 <- c(out_rand_int_bear$samples[[1]][,6], out_rand_int_bear$samples[[2]][,6], 
	out_rand_int_bear$samples[[3]][,6])

v_bear1 <- c(out_rand_int_bear$samples[[1]][,7], out_rand_int_bear$samples[[2]][,7], 
	out_rand_int_bear$samples[[3]][,7])
v_bear2 <- c(out_rand_int_bear$samples[[1]][,8], out_rand_int_bear$samples[[2]][,8], 
	out_rand_int_bear$samples[[3]][,8])
v_bear3 <- c(out_rand_int_bear$samples[[1]][,9], out_rand_int_bear$samples[[2]][,9], 
	out_rand_int_bear$samples[[3]][,9])

lambda_bear1 <- c(out_rand_int_bear$samples[[1]][,10], out_rand_int_bear$samples[[2]][,10], 
	out_rand_int_bear$samples[[3]][,10])
lambda_bear2 <- c(out_rand_int_bear$samples[[1]][,11], out_rand_int_bear$samples[[2]][,11], 
	out_rand_int_bear$samples[[3]][,11])
lambda_bear3 <- c(out_rand_int_bear$samples[[1]][,12], out_rand_int_bear$samples[[2]][,12], 
	out_rand_int_bear$samples[[3]][,12])
	

steps_bear1 <- numeric(length = length(mu_bear1))
turns_bear1 <- numeric(length = length(mu_bear1))
steps_bear2 <- numeric(length = length(mu_bear1))
turns_bear2 <- numeric(length = length(mu_bear1))
steps_bear3 <- numeric(length = length(mu_bear1))
turns_bear3 <- numeric(length = length(mu_bear1))

for(i in 1:length(steps_bear1)){
	steps_bear1[i] <- rweibull(1, v_bear1[i], (1/lambda_bear1[i])^(1/v_bear1[i]))
	turns_bear1[i] <- rwcauchy(1, mu_bear1[i], rho_bear1[i])

	steps_bear2[i] <- rweibull(1, v_bear2[i], (1/lambda_bear2[i])^(1/v_bear2[i]))
	turns_bear2[i] <- rwcauchy(1, mu_bear2[i], rho_bear2[i])
	
	steps_bear3[i] <- rweibull(1, v_bear3[i], (1/lambda_bear3[i])^(1/v_bear3[i]))
	turns_bear3[i] <- rwcauchy(1, mu_bear3[i], rho_bear3[i])
	}


sims_bear1 <- as.data.frame(cbind(steps_bear1, turns_bear1))
sims_bear1$turns_bear1[sims_bear1$turns_bear1 > pi] = sims_bear1$turns_bear1[sims_bear1$turns_bear1 > pi] - 2*pi

sims_bear2 <- as.data.frame(cbind(steps_bear2, turns_bear2))
sims_bear2$turns_bear2[sims_bear2$turns_bear2 > pi] = sims_bear2$turns_bear2[sims_bear2$ turns_bear2 > pi] - 2*pi

sims_bear3 <- as.data.frame(cbind(steps_bear3, turns_bear3))
sims_bear3$turns_bear3[sims_bear3$turns_bear3 > pi] = sims_bear3$turns_bear3[sims_bear3$turns_bear3 > pi] - 2*pi

keep <- sample(1:45000, 1000, replace = F)
sims_bear1 <- sims_bear1[keep,] %>%
	mutate(index = 1) %>%
	mutate(iter = 1:nrow(.))
colnames(sims_bear1) <- c("steps", "turns", "index", "iter")
sims_bear2 <- sims_bear2[keep,] %>%
	mutate(index = 2) %>%
	mutate(iter= 1:nrow(.))
colnames(sims_bear2) <- c("steps", "turns", "index", "iter")
sims_bear3 <- sims_bear3[keep,] %>%
	mutate(index = 3) %>%
	mutate(iter= 1:nrow(.))
colnames(sims_bear3) <- c("steps", "turns", "index", "iter")

sims_bear <- sims_bear1 %>%
	rbind(sims_bear2) %>%
	rbind(sims_bear3) %>%
	mutate(turns_deg = NISTunits::NISTradianTOdeg(turns))


#  Density plots
steps_p <- ggplot(sims_bear, aes(steps, group = factor(index))) + 
	geom_density(aes(linetype = factor(index)), size = 0.75) +
	scale_linetype_manual(name = "Position Index", values = c(1, 6, 3)) +
	xlab("\n Step length (km)") +
	ylab("Density \n") +
	xlim(c(0, 3)) +
	ylim(c(0, 2.5)) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		legend.position = c(0.88, 0.88)) 
steps_p

turns_p <- ggplot(sims_bear, aes(turns_deg, group = factor(index))) + 
	geom_density(aes(linetype = factor(index)), size = 0.75) +
	scale_linetype_manual(name = "Distance Index", values = c(1, 6, 3)) +
	xlab("\n Turn angle (deg)") +
	ylab("Density \n") +
	xlim(c(-180, 180)) +
	ylim(c(0, 0.005)) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		legend.position = c(0.88, 0.88)) 
turns_p
