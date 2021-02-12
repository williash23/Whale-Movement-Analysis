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
library(patchwork)


################################################################################
# Load output from MCMC posterior distribtuion iterations (JAGS object)
load("C:/Users/saraw/OneDrive/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/out_dist.Rdata")
load("C:/Users/saraw/OneDrive/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/models/out_bear.Rdata")
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
mu_dist1 <- c(out_dist$samples[[1]][,1], out_dist$samples[[2]][,1], 
	out_dist$samples[[3]][,1])
mu_dist2 <- c(out_dist$samples[[1]][,2], out_dist$samples[[2]][,2], 
	out_dist$samples[[3]][,2])
mu_dist3 <- c(out_dist$samples[[1]][,3], out_dist$samples[[2]][,3], 
	out_dist$samples[[3]][,3])

rho_dist1 <- c(out_dist$samples[[1]][,4], out_dist$samples[[2]][,4], 
	out_dist$samples[[3]][,4])
rho_dist2 <- c(out_dist$samples[[1]][,5], out_dist$samples[[2]][,5], 
	out_dist$samples[[3]][,5])
rho_dist3 <- c(out_dist$samples[[1]][,6], out_dist$samples[[2]][,6], 
	out_dist$samples[[3]][,6])

v_dist1 <- c(out_dist$samples[[1]][,7], out_dist$samples[[2]][,7], 
	out_dist$samples[[3]][,7])
v_dist2 <- c(out_dist$samples[[1]][,8], out_dist$samples[[2]][,8], 
	out_dist$samples[[3]][,8])
v_dist3 <- c(out_dist$samples[[1]][,9], out_dist$samples[[2]][,9], 
	out_dist$samples[[3]][,9])

lambda_dist1 <- c(out_dist$samples[[1]][,10], out_dist$samples[[2]][,10], 
	out_dist$samples[[3]][,10])
lambda_dist2 <- c(out_dist$samples[[1]][,11], out_dist$samples[[2]][,11], 
	out_dist$samples[[3]][,11])
lambda_dist3 <- c(out_dist$samples[[1]][,12], out_dist$samples[[2]][,12], 
	out_dist$samples[[3]][,12])
	

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

keep <- sample(1:29700, 10000, replace = F)
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
	scale_linetype_manual(name = "Distance Index", values = c(1, 6, 3),
		labels = c("Near", "Intermediate", "Far")) +
	xlab("\n Step length (km)") +
	ylab("Density \n") +
	xlim(c(0, 3)) +
	ylim(c(0, 6)) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		legend.position = c(0.85, 0.88)) 
steps_p

turns_p <- ggplot(sims, aes(turns_deg, group = factor(index))) + 
	geom_density(aes(linetype = factor(index)), size = 0.75) +
	scale_linetype_manual(name = "Distance Index", values = c(1, 6, 3),
		labels = c("Near", "Intermediate", "Far")) +
	xlab("\n Turn angle (deg)") +
	ylab("Density \n") +
	xlim(c(-180, 180)) +
	ylim(c(0, 0.005)) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		legend.position = c(0.88, 0.88)) 
turns_p

steps_p / turns_p +  plot_layout(guides = "collect") & theme(legend.position = 'top')




# Parameter values - Position
mu_bear1 <- c(out_bear$samples[[1]][,1], out_bear$samples[[2]][,1], 
	out_bear$samples[[3]][,1])
mu_bear2 <- c(out_bear$samples[[1]][,2], out_bear$samples[[2]][,2], 
	out_bear$samples[[3]][,2])
mu_bear3 <- c(out_bear$samples[[1]][,3], out_bear$samples[[2]][,3], 
	out_bear$samples[[3]][,3])

rho_bear1 <- c(out_bear$samples[[1]][,4], out_bear$samples[[2]][,4], 
	out_bear$samples[[3]][,4])
rho_bear2 <- c(out_bear$samples[[1]][,5], out_bear$samples[[2]][,5], 
	out_bear$samples[[3]][,5])
rho_bear3 <- c(out_bear$samples[[1]][,6], out_bear$samples[[2]][,6], 
	out_bear$samples[[3]][,6])

v_bear1 <- c(out_bear$samples[[1]][,7], out_bear$samples[[2]][,7], 
	out_bear$samples[[3]][,7])
v_bear2 <- c(out_bear$samples[[1]][,8], out_bear$samples[[2]][,8], 
	out_bear$samples[[3]][,8])
v_bear3 <- c(out_bear$samples[[1]][,9], out_bear$samples[[2]][,9], 
	out_bear$samples[[3]][,9])

lambda_bear1 <- c(out_bear$samples[[1]][,10], out_bear$samples[[2]][,10], 
	out_bear$samples[[3]][,10])
lambda_bear2 <- c(out_bear$samples[[1]][,11], out_bear$samples[[2]][,11], 
	out_bear$samples[[3]][,11])
lambda_bear3 <- c(out_bear$samples[[1]][,12], out_bear$samples[[2]][,12], 
	out_bear$samples[[3]][,12])
	

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

keep <- sample(1:29700, 10000, replace = F)
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
	scale_linetype_manual(name = "Position Index", values = c(1, 6, 3),
		labels = c("Forward", "Mid", "Abeam")) +
	xlab("\n Step length (km)") +
	ylab("Density \n") +
	xlim(c(0, 3)) +
	ylim(c(0, 3)) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		legend.position = c(0.88, 0.88)) 
steps_p

turns_p <- ggplot(sims_bear, aes(turns_deg, group = factor(index))) + 
	geom_density(aes(linetype = factor(index)), size = 0.75) +
	scale_linetype_manual(name = "Position Index", values = c(1, 6, 3),
		labels = c("Forward", "Mid", "Abeam")) +
	xlab("\n Turn angle (deg)") +
	ylab("Density \n") +
	xlim(c(-180, 180)) +
	ylim(c(0, 0.005)) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		legend.position = c(0.88, 0.88)) 
turns_p



steps_p / turns_p +  plot_layout(guides = "collect") & theme(legend.position = 'top')





# Parameter values - Combined (4 indices)
mu_comb1 <- c(out_comb$samples[[1]][,1], out_comb$samples[[2]][,1], 
	out_comb$samples[[3]][,1])
mu_comb2 <- c(out_comb$samples[[1]][,2], out_comb$samples[[2]][,2], 
	out_comb$samples[[3]][,2])
mu_comb3 <- c(out_comb$samples[[1]][,3], out_comb$samples[[2]][,3], 
	out_comb$samples[[3]][,3])
mu_comb4 <- c(out_comb$samples[[1]][,4], out_comb$samples[[2]][,4], 
	out_comb$samples[[3]][,4])
	
rho_comb1 <- c(out_comb$samples[[1]][,5], out_comb$samples[[2]][,5], 
	out_comb$samples[[3]][,5])
rho_comb2 <- c(out_comb$samples[[1]][,6], out_comb$samples[[2]][,6], 
	out_comb$samples[[3]][,6])
rho_comb3 <- c(out_comb$samples[[1]][,7], out_comb$samples[[2]][,7], 
	out_comb$samples[[3]][,7])
rho_comb4 <- c(out_comb$samples[[1]][,8], out_comb$samples[[2]][,8], 
	out_comb$samples[[3]][,8])

v_comb1 <- c(out_comb$samples[[1]][,9], out_comb$samples[[2]][,9], 
	out_comb$samples[[3]][,9])
v_comb2 <- c(out_comb$samples[[1]][,10], out_comb$samples[[2]][,10], 
	out_comb$samples[[3]][,10])
v_comb3 <- c(out_comb$samples[[1]][,11], out_comb$samples[[2]][,11], 
	out_comb$samples[[3]][,11])
v_comb4 <- c(out_comb$samples[[1]][,12], out_comb$samples[[2]][,12], 
	out_comb$samples[[3]][,12])

lambda_comb1 <- c(out_comb$samples[[1]][,13], out_comb$samples[[2]][,13], 
	out_comb$samples[[3]][,13])
lambda_comb2 <- c(out_comb$samples[[1]][,14], out_comb$samples[[2]][,14], 
	out_comb$samples[[3]][,14])
lambda_comb3 <- c(out_comb$samples[[1]][,15], out_comb$samples[[2]][,15], 
	out_comb$samples[[3]][,15])
lambda_comb4 <- c(out_comb$samples[[1]][,16], out_comb$samples[[2]][,16], 
	out_comb$samples[[3]][,16])
	

steps_comb1 <- numeric(length = length(mu_comb1))
turns_comb1 <- numeric(length = length(mu_comb1))
steps_comb2 <- numeric(length = length(mu_comb1))
turns_comb2 <- numeric(length = length(mu_comb1))
steps_comb3 <- numeric(length = length(mu_comb1))
turns_comb3 <- numeric(length = length(mu_comb1))
steps_comb4 <- numeric(length = length(mu_comb1))
turns_comb4 <- numeric(length = length(mu_comb1))


for(i in 1:length(steps_comb1)){
	steps_comb1[i] <- rweibull(1, v_comb1[i], (1/lambda_comb1[i])^(1/v_comb1[i]))
	turns_comb1[i] <- rwcauchy(1, mu_comb1[i], rho_comb1[i])

	steps_comb2[i] <- rweibull(1, v_comb2[i], (1/lambda_comb2[i])^(1/v_comb2[i]))
	turns_comb2[i] <- rwcauchy(1, mu_comb2[i], rho_comb2[i])
	
	steps_comb3[i] <- rweibull(1, v_comb3[i], (1/lambda_comb3[i])^(1/v_comb3[i]))
	turns_comb3[i] <- rwcauchy(1, mu_comb3[i], rho_comb3[i])
	
	steps_comb4[i] <- rweibull(1, v_comb4[i], (1/lambda_comb4[i])^(1/v_comb4[i]))
	turns_comb4[i] <- rwcauchy(1, mu_comb4[i], rho_comb4[i])
	}


sims_comb1 <- as.data.frame(cbind(steps_comb1, turns_comb1))
sims_comb1$turns_comb1[sims_comb1$turns_comb1 > pi] = sims_comb1$turns_comb1[sims_comb1$turns_comb1 > pi] - 2*pi

sims_comb2 <- as.data.frame(cbind(steps_comb2, turns_comb2))
sims_comb2$turns_comb2[sims_comb2$turns_comb2 > pi] = sims_comb2$turns_comb2[sims_comb2$ turns_comb2 > pi] - 2*pi

sims_comb3 <- as.data.frame(cbind(steps_comb3, turns_comb3))
sims_comb3$turns_comb3[sims_comb3$turns_comb3 > pi] = sims_comb3$turns_comb3[sims_comb3$turns_comb3 > pi] - 2*pi

sims_comb4 <- as.data.frame(cbind(steps_comb4, turns_comb4))
sims_comb4$turns_comb4[sims_comb4$turns_comb4 > pi] = sims_comb4$turns_comb4[sims_comb4$turns_comb4 > pi] - 2*pi


keep <- sample(1:13200, 10000, replace = F)
sims_comb1 <- sims_comb1[keep,] %>%
	mutate(index = "Near/Forward") %>%
	mutate(iter = 1:nrow(.))
colnames(sims_comb1) <- c("steps", "turns", "index", "iter")
sims_comb2 <- sims_comb2[keep,] %>%
	mutate(index = "Near/Abeam") %>%
	mutate(iter= 1:nrow(.))
colnames(sims_comb2) <- c("steps", "turns", "index", "iter")
sims_comb3 <- sims_comb3[keep,] %>%
	mutate(index = "Far/Forward") %>%
	mutate(iter= 1:nrow(.))
colnames(sims_comb3) <- c("steps", "turns", "index", "iter")
sims_comb4 <- sims_comb4[keep,] %>%
	mutate(index = "Far/Abeam") %>%
	mutate(iter= 1:nrow(.))
colnames(sims_comb4) <- c("steps", "turns", "index", "iter")


sims_comb <- sims_comb1 %>%
	rbind(sims_comb2) %>%
	rbind(sims_comb3) %>%
	rbind(sims_comb4) %>%
	mutate(turns_deg = NISTunits::NISTradianTOdeg(turns)) %>%
	mutate(abs_turns_deg = abs(turns_deg))


#  Density plots
#steps_p <- ggplot(sims_comb, aes(steps, group = factor(index))) + 
steps_p <- ggplot(sims_comb, aes(steps)) + 
	#geom_density(aes(linetype = factor(index)), size = 0.75) +
	geom_density(size = 0.5) +
	#scale_linetype_manual(name = "Position Index", values = c(1, 6, 3),
		#labels = c("Forward", "Mid", "Abeam")) +
	xlab("\n Surfacing length (km)") +
	ylab("Density \n") +
	scale_x_continuous(expand = c(0, 0)) + 
	scale_y_continuous(expand = c(0, 0), limits = c(0, 12)) +
	#xlim(c(0, 3)) +
	#ylim(c(0, 3)) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		legend.position = c(0.88, 0.88)) +
	facet_wrap(~index)
steps_p

#turns_p <- ggplot(sims_comb, aes(turns_deg), group = factor(index))) + 
turns_p <- ggplot(sims_comb, aes(turns_deg)) + 
	#geom_density(aes(linetype = factor(index)), size = 0.75) +
	#geom_histogram() +
	geom_density(size = 0.5) +
	#scale_linetype_manual(name = "Position Index", values = c(1, 6, 3),
		#labels = c("Forward", "Mid", "Abeam")) +
	xlab("\n Absolute surfacing angle (deg)") +
	ylab("Density \n") +
	#xlim(c(-180, 180)) +
	#ylim(c(0, 0.005)) +
	scale_x_continuous(expand = c(0, 0), limits = c(-185, 185)) + 
	scale_y_continuous(expand = c(0, 0), limits = c(0, 30)) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		legend.position = c(0.88, 0.88)) +
	facet_wrap(~index)
turns_p



steps_p / turns_p +  plot_layout(guides = "collect") & theme(legend.position = 'top')









############################################################################

# Simulation 
ship_spd <- 9.8 #13 kts
# Number of whale surfacings
n_surface <- 2
# Ship movement
num_ship_pts <- 2
# Mean XY UTM location of all whales
init_x <- 438456
init_y <- 6486336 

sim_seq <- seq(from = 3, to = 10000, by = 2)

x <- raster(ncol=1000, nrow=1000, xmn=435000, xmx=441000, ymn=6485000, ymx=6488000)
projection(x) <- "+proj=utm +zone=8 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"



sims_bear1 <- as.data.frame(cbind(steps_bear1, turns_bear1))
sims_bear1$turns_bear1[sims_bear1$turns_bear1 > pi] = sims_bear1$turns_bear1[sims_bear1$turns_bear1 > pi] - 2*pi

sims_bear3 <- as.data.frame(cbind(steps_bear3, turns_bear3))
sims_bear3$turns_bear3[sims_bear3$turns_bear3 > pi] = sims_bear3$turns_bear3[sims_bear3$turns_bear3 > pi] - 2*pi



# Simulate over bearing index 1 

beta <- NISTunits::NISTdegTOradian(10)
radial_dist <- 3000
vert_dist <- radial_dist * sin(beta)
perp_dist <- sqrt((radial_dist^2 - vert_dist^2))


# Ship starting location
init_ship_x <- init_x - perp_dist
init_ship_y <- init_y + vert_dist

colnames(ship_init_loc_df) <- c("x", "y")
ship_loc1 <- c(init_ship_x, init_ship_y)
ship_loc2 <- c(ship_loc1[1] + (perp_dist/2), init_ship_y)
ship_loc3 <- c(ship_loc2[1] + (perp_dist/2), init_ship_y)

ship_loc_df <- as.data.frame(rbind(ship_loc1, ship_loc2, ship_loc3))
colnames(ship_loc_df) <- c("x", "y")
coordinates(ship_loc_df) <- ~x+y
ship_loc_sf <- st_as_sf(ship_loc_df)  %>%
	st_set_crs("+proj=utm +zone=8 +datum=WGS84") 	
ship_loc_sf_bear1 <- ship_loc_sf	
	
strikes <-numeric()



# WHALE
rnd_iter1 <- runif(1, 1, 29700)
rnd_iter2 <- runif(1, 1, 29700)
step1 <- sims_bear1$step[rnd_iter1]*1000
step2 <- sims_bear1$step[rnd_iter2]*1000
turn1 <- sims_bear1$turn[rnd_iter1]
turn2 <- sims_bear1$turn[rnd_iter2]

whale_loc1 <- c(init_x, init_y)
whale_move1_perp_dist <- (step1 * sin(turn1))
whale_move1_vert_dist <- sqrt(step1^2 - whale_move1_perp_dist^2)
whale_loc2 <- c(whale_loc1[1] + whale_move1_perp_dist, whale_loc1[2] + whale_move1_vert_dist)
whale_move2_perp_dist <- (step2 * sin(turn2))
whale_move2_vert_dist <- sqrt(step2^2 - whale_move2_perp_dist^2)
whale_loc3 <- c(whale_loc2[1] + whale_move2_perp_dist, whale_loc2[2] + whale_move2_vert_dist)

whale_loc_df <- as.data.frame(rbind(whale_loc1, whale_loc2, whale_loc3))
colnames(whale_loc_df) <- c("x", "y")
coordinates(whale_loc_df) <- ~x+y
whale_loc_sf <- st_as_sf(whale_loc_df)  %>%
	st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
	mutate(obs_num = 1:n()) %>%
	mutate(iter_num = 1)
	
#strike1 <- st_intersects(st_buffer(ship_loc_sf[2,], 50), st_buffer(add_sf[2,], 10), sparse = FALSE)
#strike1 <- FALSE
strike2 <- st_intersects(st_buffer(ship_loc_sf[3,], 50), st_buffer(add_sf[3,], 10), sparse = FALSE)
#strikes_tmp <- c(strike1, strike2)
strikes_tmp <- strike2
n_strike <- length(which(strikes_tmp == TRUE))

strikes[1] <- n_strike


for(i in 1:10000){
	rnd_iter1 <- runif(1, 1, 29700)
	rnd_iter2 <- runif(1, 1, 29700)
	step1 <- sims_bear1$step[rnd_iter1]*1000
	step2 <- sims_bear1$step[rnd_iter2]*1000
	turn1 <- sims_bear1$turn[rnd_iter1]
	turn2 <- sims_bear1$turn[rnd_iter2]

	whale_loc1 <- c(init_x, init_y)
	whale_move1_perp_dist <- (step1 * sin(turn1))
	whale_move1_vert_dist <- sqrt(step1^2 - whale_move1_perp_dist^2)
	whale_loc2 <- c(whale_loc1[1] + whale_move1_perp_dist, whale_loc1[2] + whale_move1_vert_dist)
	whale_move2_perp_dist <- (step2 * sin(turn2))
	whale_move2_vert_dist <- sqrt(step2^2 - whale_move2_perp_dist^2)
	whale_loc3 <- c(whale_loc2[1] + whale_move2_perp_dist, whale_loc2[2] + whale_move2_vert_dist)

	whale_loc_df <- as.data.frame(rbind(whale_loc1, whale_loc2, whale_loc3))
	colnames(whale_loc_df) <- c("x", "y")
	coordinates(whale_loc_df) <- ~x+y
	add_sf <- st_as_sf(whale_loc_df)  %>%
		st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
		mutate(obs_num = 1:n()) %>%
		mutate(iter_num = i)
	
	#strike1 <- st_intersects(st_buffer(ship_loc_sf[2,], 50), st_buffer(add_sf[2,], 10), sparse = FALSE)
	#strike1 <- FALSE
	strike2 <- st_intersects(st_buffer(ship_loc_sf[3,], 50), st_buffer(add_sf[3,], 10), sparse = FALSE)
	#strikes_tmp <- c(strike1, strike2)
	strikes_tmp <- strike2
	n_strike <- length(which(strikes_tmp == TRUE))
	strikes[i] <- n_strike

	whale_loc_sf <- rbind(whale_loc_sf, add_sf)
	
	}

whale_loc_sf_bear1 <- whale_loc_sf 	
strk_bear1 <- sum(strikes)




# Simulate over bearing index 3
beta <- NISTunits::NISTdegTOradian(40)
radial_dist <- 810
vert_dist <- radial_dist * sin(beta)
perp_dist <- sqrt((radial_dist^2 - vert_dist^2))


# Ship starting location
init_ship_x <- init_x - perp_dist
init_ship_y <- init_y + vert_dist

colnames(ship_init_loc_df) <- c("x", "y")
ship_loc1 <- c(init_ship_x, init_ship_y)
ship_loc2 <- c(ship_loc1[1] + (perp_dist/2), init_ship_y)
ship_loc3 <- c(ship_loc2[1] + (perp_dist/2), init_ship_y)

ship_loc_df <- as.data.frame(rbind(ship_loc1, ship_loc2, ship_loc3))
colnames(ship_loc_df) <- c("x", "y")
coordinates(ship_loc_df) <- ~x+y
ship_loc_sf <- st_as_sf(ship_loc_df)  %>%
	st_set_crs("+proj=utm +zone=8 +datum=WGS84") 	
ship_loc_sf_bear3 <- ship_loc_sf	

strikes <- numeric()

# WHALE
rnd_iter1 <- runif(1, 1, 29700)
rnd_iter2 <- runif(1, 1, 29700)
step1 <- sims_bear3$step[rnd_iter1]*1000
step2 <- sims_bear3$step[rnd_iter2]*1000
turn1 <- sims_bear3$turn[rnd_iter1]
turn2 <- sims_bear3$turn[rnd_iter2]

whale_loc1 <- c(init_x, init_y)
whale_move1_perp_dist <- (step1 * sin(turn1))
whale_move1_vert_dist <- sqrt(step1^2 - whale_move1_perp_dist^2)
whale_loc2 <- c(whale_loc1[1] + whale_move1_perp_dist, whale_loc1[2] + whale_move1_vert_dist)
whale_move2_perp_dist <- (step2 * sin(turn2))
whale_move2_vert_dist <- sqrt(step2^2 - whale_move2_perp_dist^2)
whale_loc3 <- c(whale_loc2[1] + whale_move2_perp_dist, whale_loc2[2] + whale_move2_vert_dist)

whale_loc_df <- as.data.frame(rbind(whale_loc1, whale_loc2, whale_loc3))
colnames(whale_loc_df) <- c("x", "y")
coordinates(whale_loc_df) <- ~x+y
whale_loc_sf <- st_as_sf(whale_loc_df)  %>%
	st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
	mutate(obs_num = 1:n()) %>%
	mutate(iter_num = 1)
	
#strike1 <- st_intersects(st_buffer(ship_loc_sf[2,], 50), st_buffer(add_sf[2,], 10), sparse = FALSE)
#strike1 <- FALSE
strike2 <- st_intersects(st_buffer(ship_loc_sf[3,], 50), st_buffer(add_sf[3,], 10), sparse = FALSE)
#strikes_tmp <- c(strike1, strike2)
strikes_tmp <- strike2
n_strike <- length(which(strikes_tmp == TRUE))

strikes[1] <- n_strike

for(i in 1:10000){
	rnd_iter1 <- runif(1, 1, 29700)
	rnd_iter2 <- runif(1, 1, 29700)
	
	step1 <- sims_bear3$step[rnd_iter1]*1000
	step2 <- sims_bear3$step[rnd_iter2]*1000
	turn1 <- sims_bear3$turn[rnd_iter1]
	turn2 <- sims_bear3$turn[rnd_iter2]

	whale_loc1 <- c(init_x, init_y)
	whale_move1_perp_dist <- (step1 * sin(turn1))
	whale_move1_vert_dist <- sqrt(step1^2 - whale_move1_perp_dist^2)
	whale_loc2 <- c(whale_loc1[1] + whale_move1_perp_dist, whale_loc1[2] + whale_move1_vert_dist)
	whale_move2_perp_dist <- (step2 * sin(turn2))
	whale_move2_vert_dist <- sqrt(step2^2 - whale_move2_perp_dist^2)
	whale_loc3 <- c(whale_loc2[1] + whale_move2_perp_dist, whale_loc2[2] + whale_move2_vert_dist)

	whale_loc_df <- as.data.frame(rbind(whale_loc1, whale_loc2, whale_loc3))
	colnames(whale_loc_df) <- c("x", "y")
	coordinates(whale_loc_df) <- ~x+y
	add_sf <- st_as_sf(whale_loc_df)  %>%
		st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
		mutate(obs_num = 1:n()) %>%
		mutate(iter_num = i)
	
	#strike1 <- st_intersects(st_buffer(ship_loc_sf[2,], 50), st_buffer(add_sf[2,], 10), sparse = FALSE)
	#strike1 <- FALSE
	strike2 <- st_intersects(st_buffer(ship_loc_sf[3,], 50), st_buffer(add_sf[3,], 10), sparse = FALSE)
	#strikes_tmp <- c(strike1, strike2)
	strikes_tmp <- strike2
	n_strike <- length(which(strikes_tmp == TRUE))

	strikes[i] <- n_strike

	whale_loc_sf <- rbind(whale_loc_sf, add_sf)
	
	}

whale_loc_sf_bear3 <- whale_loc_sf 	%>%
	mutate(val = 1)
strk_bear3 <- sum(strikes)


bear1_r <- rasterize(whale_loc_sf_bear1, x, field = 1, fun = 'count')
bear1_r01 <- climateStability::rescale0to1(bear1_r)
bear1_sf <- st_as_sf(rasterToPolygons(bear1_r01)) 

bear3_r <- rasterize(whale_loc_sf_bear3, x, field = 1, fun = 'count')
bear3_r01 <- climateStability::rescale0to1(bear3_r)
bear3_sf <- st_as_sf(rasterToPolygons(bear3_r01)) 



p_bear1 <- ggplot() +
	geom_sf(data = bear1_sf, color = "transparent", fill = "darkred", aes(alpha = layer)) +
	geom_sf(data = ship_loc_sf_bear1, fill = "black", color = "black", shape = 16) +
	geom_sf(data = whale_loc_sf_bear1[1,], fill = "black", color = "black", shape = 17, size = 2) +
	scale_alpha(name = "Relative prob.\n of surfacing",
		range = c(0.05,1)) +
	coord_sf(datum = 32608) +
	#xlim(c(435500, 442000)) +
	#ylim(c(6486000, 6490000)) +
	theme_bw()
#p_bear1


p_bear3 <- ggplot() +
	geom_sf(data = bear3_sf, color = "transparent", fill = "darkred", aes(alpha = layer)) +
	geom_sf(data = ship_loc_sf_bear3, fill = "black", color = "black", shape = 16) +
	geom_sf(data = whale_loc_sf_bear3[1,], fill = "black", color = "black", shape = 17, size = 2) +
	scale_alpha(name = "Relative prob.\n of surfacing",
		range = c(0.05,1)) +
	coord_sf(datum = 32608) +
	#xlim(c(435500, 442000)) +
	#ylim(c(6486000, 6490000)) +
	theme_bw() + 
	theme(axis.text.y = element_blank())
#p_bear3


p_bear1 + p_bear3  + plot_layout(guides = "collect")




sims_dist1 <- as.data.frame(cbind(steps_dist1, turns_dist1))
sims_dist1$turns_dist1[sims_dist1$turns_dist1 > pi] = sims_dist1$turns_dist1[sims_dist1$turns_dist1 > pi] - 2*pi

sims_dist3 <- as.data.frame(cbind(steps_dist3, turns_dist3))
sims_dist3$turns_dist3[sims_dist3$turns_dist3 > pi] = sims_dist3$turns_dist3[sims_dist3$turns_dist3 > pi] - 2*pi


# Simulate over distance index 1 

beta <- NISTunits::NISTdegTOradian(10)
radial_dist <- 500
vert_dist <- radial_dist * sin(beta)
perp_dist <- sqrt((radial_dist^2 - vert_dist^2))


# Ship starting location
init_ship_x <- init_x - perp_dist
init_ship_y <- init_y + vert_dist

colnames(ship_init_loc_df) <- c("x", "y")
ship_loc1 <- c(init_ship_x, init_ship_y)
ship_loc2 <- c(ship_loc1[1] + (perp_dist/2), init_ship_y)
ship_loc3 <- c(ship_loc2[1] + (perp_dist/2), init_ship_y)

ship_loc_df <- as.data.frame(rbind(ship_loc1, ship_loc2, ship_loc3))
colnames(ship_loc_df) <- c("x", "y")
coordinates(ship_loc_df) <- ~x+y
ship_loc_sf <- st_as_sf(ship_loc_df)  %>%
	st_set_crs("+proj=utm +zone=8 +datum=WGS84") 	
ship_loc_sf_dist1 <- ship_loc_sf	
	
strikes <-numeric()



# WHALE
rnd_iter1 <- runif(1, 1, 29700)
rnd_iter2 <- runif(1, 1, 29700)
step1 <- sims_dist1$step[rnd_iter1]*1000
step2 <- sims_dist1$step[rnd_iter2]*1000
turn1 <- sims_dist1$turn[rnd_iter1]
turn2 <- sims_dist1$turn[rnd_iter2]

whale_loc1 <- c(init_x, init_y)
whale_move1_perp_dist <- (step1 * sin(turn1))
whale_move1_vert_dist <- sqrt(step1^2 - whale_move1_perp_dist^2)
whale_loc2 <- c(whale_loc1[1] + whale_move1_perp_dist, whale_loc1[2] + whale_move1_vert_dist)
whale_move2_perp_dist <- (step2 * sin(turn2))
whale_move2_vert_dist <- sqrt(step2^2 - whale_move2_perp_dist^2)
whale_loc3 <- c(whale_loc2[1] + whale_move2_perp_dist, whale_loc2[2] + whale_move2_vert_dist)

whale_loc_df <- as.data.frame(rbind(whale_loc1, whale_loc2, whale_loc3))
colnames(whale_loc_df) <- c("x", "y")
coordinates(whale_loc_df) <- ~x+y
whale_loc_sf <- st_as_sf(whale_loc_df)  %>%
	st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
	mutate(obs_num = 1:n()) %>%
	mutate(iter_num = 1)
	
#strike1 <- st_intersects(st_buffer(ship_loc_sf[2,], 50), st_buffer(add_sf[2,], 10), sparse = FALSE)
#strike1 <- FALSE
strike2 <- st_intersects(st_buffer(ship_loc_sf[3,], 50), st_buffer(add_sf[3,], 10), sparse = FALSE)
#strikes_tmp <- c(strike1, strike2)
strikes_tmp <- strike2
n_strike <- length(which(strikes_tmp == TRUE))

strikes[1] <- n_strike


for(i in 1:10000){
	rnd_iter1 <- runif(1, 1, 29700)
	rnd_iter2 <- runif(1, 1, 29700)
	step1 <- sims_dist1$step[rnd_iter1]*1000
	step2 <- sims_dist1$step[rnd_iter2]*1000
	turn1 <- sims_dist1$turn[rnd_iter1]
	turn2 <- sims_dist1$turn[rnd_iter2]

	whale_loc1 <- c(init_x, init_y)
	whale_move1_perp_dist <- (step1 * sin(turn1))
	whale_move1_vert_dist <- sqrt(step1^2 - whale_move1_perp_dist^2)
	whale_loc2 <- c(whale_loc1[1] + whale_move1_perp_dist, whale_loc1[2] + whale_move1_vert_dist)
	whale_move2_perp_dist <- (step2 * sin(turn2))
	whale_move2_vert_dist <- sqrt(step2^2 - whale_move2_perp_dist^2)
	whale_loc3 <- c(whale_loc2[1] + whale_move2_perp_dist, whale_loc2[2] + whale_move2_vert_dist)

	whale_loc_df <- as.data.frame(rbind(whale_loc1, whale_loc2, whale_loc3))
	colnames(whale_loc_df) <- c("x", "y")
	coordinates(whale_loc_df) <- ~x+y
	add_sf <- st_as_sf(whale_loc_df)  %>%
		st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
		mutate(obs_num = 1:n()) %>%
		mutate(iter_num = i)
	
	#strike1 <- st_intersects(st_buffer(ship_loc_sf[2,], 50), st_buffer(add_sf[2,], 10), sparse = FALSE)
	#strike1 <- FALSE
	strike2 <- st_intersects(st_buffer(ship_loc_sf[3,], 50), st_buffer(add_sf[3,], 10), sparse = FALSE)
	#strikes_tmp <- c(strike1, strike2)
	strikes_tmp <- strike2
	n_strike <- length(which(strikes_tmp == TRUE))
	strikes[i] <- n_strike

	whale_loc_sf <- rbind(whale_loc_sf, add_sf)
	
	}

whale_loc_sf_dist1 <- whale_loc_sf 	
strk_dist1 <- sum(strikes)




# Simulate over disting index 3
beta <- NISTunits::NISTdegTOradian(30)
radial_dist <- 4500
vert_dist <- radial_dist * sin(beta)
perp_dist <- sqrt((radial_dist^2 - vert_dist^2))


# Ship starting location
init_ship_x <- init_x - perp_dist
init_ship_y <- init_y + vert_dist

colnames(ship_init_loc_df) <- c("x", "y")
ship_loc1 <- c(init_ship_x, init_ship_y)
ship_loc2 <- c(ship_loc1[1] + (perp_dist/2), init_ship_y)
ship_loc3 <- c(ship_loc2[1] + (perp_dist/2), init_ship_y)

ship_loc_df <- as.data.frame(rbind(ship_loc1, ship_loc2, ship_loc3))
colnames(ship_loc_df) <- c("x", "y")
coordinates(ship_loc_df) <- ~x+y
ship_loc_sf <- st_as_sf(ship_loc_df)  %>%
	st_set_crs("+proj=utm +zone=8 +datum=WGS84") 	
ship_loc_sf_dist3 <- ship_loc_sf	

strikes <- numeric()

# WHALE
rnd_iter1 <- runif(1, 1, 29700)
rnd_iter2 <- runif(1, 1, 29700)
step1 <- sims_dist3$step[rnd_iter1]*1000
step2 <- sims_dist3$step[rnd_iter2]*1000
turn1 <- sims_dist3$turn[rnd_iter1]
turn2 <- sims_dist3$turn[rnd_iter2]

whale_loc1 <- c(init_x, init_y)
whale_move1_perp_dist <- (step1 * sin(turn1))
whale_move1_vert_dist <- sqrt(step1^2 - whale_move1_perp_dist^2)
whale_loc2 <- c(whale_loc1[1] + whale_move1_perp_dist, whale_loc1[2] + whale_move1_vert_dist)
whale_move2_perp_dist <- (step2 * sin(turn2))
whale_move2_vert_dist <- sqrt(step2^2 - whale_move2_perp_dist^2)
whale_loc3 <- c(whale_loc2[1] + whale_move2_perp_dist, whale_loc2[2] + whale_move2_vert_dist)

whale_loc_df <- as.data.frame(rbind(whale_loc1, whale_loc2, whale_loc3))
colnames(whale_loc_df) <- c("x", "y")
coordinates(whale_loc_df) <- ~x+y
whale_loc_sf <- st_as_sf(whale_loc_df)  %>%
	st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
	mutate(obs_num = 1:n()) %>%
	mutate(iter_num = 1)
	
#strike1 <- st_intersects(st_buffer(ship_loc_sf[2,], 50), st_buffer(add_sf[2,], 10), sparse = FALSE)
#strike1 <- FALSE
strike2 <- st_intersects(st_buffer(ship_loc_sf[3,], 50), st_buffer(add_sf[3,], 10), sparse = FALSE)
#strikes_tmp <- c(strike1, strike2)
strikes_tmp <- strike2
n_strike <- length(which(strikes_tmp == TRUE))

strikes[1] <- n_strike

for(i in 1:10000){
	rnd_iter1 <- runif(1, 1, 29700)
	rnd_iter2 <- runif(1, 1, 29700)
	
	step1 <- sims_dist3$step[rnd_iter1]*1000
	step2 <- sims_dist3$step[rnd_iter2]*1000
	turn1 <- sims_dist3$turn[rnd_iter1]
	turn2 <- sims_dist3$turn[rnd_iter2]

	whale_loc1 <- c(init_x, init_y)
	whale_move1_perp_dist <- (step1 * sin(turn1))
	whale_move1_vert_dist <- sqrt(step1^2 - whale_move1_perp_dist^2)
	whale_loc2 <- c(whale_loc1[1] + whale_move1_perp_dist, whale_loc1[2] + whale_move1_vert_dist)
	whale_move2_perp_dist <- (step2 * sin(turn2))
	whale_move2_vert_dist <- sqrt(step2^2 - whale_move2_perp_dist^2)
	whale_loc3 <- c(whale_loc2[1] + whale_move2_perp_dist, whale_loc2[2] + whale_move2_vert_dist)

	whale_loc_df <- as.data.frame(rbind(whale_loc1, whale_loc2, whale_loc3))
	colnames(whale_loc_df) <- c("x", "y")
	coordinates(whale_loc_df) <- ~x+y
	add_sf <- st_as_sf(whale_loc_df)  %>%
		st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
		mutate(obs_num = 1:n()) %>%
		mutate(iter_num = i)
	
	#strike1 <- st_intersects(st_buffer(ship_loc_sf[2,], 50), st_buffer(add_sf[2,], 10), sparse = FALSE)
	#strike1 <- FALSE
	strike2 <- st_intersects(st_buffer(ship_loc_sf[3,], 50), st_buffer(add_sf[3,], 10), sparse = FALSE)
	#strikes_tmp <- c(strike1, strike2)
	strikes_tmp <- strike2
	n_strike <- length(which(strikes_tmp == TRUE))

	strikes[i] <- n_strike

	whale_loc_sf <- rbind(whale_loc_sf, add_sf)
	
	}

whale_loc_sf_dist3 <- whale_loc_sf 	%>%
	mutate(val = 1)
strk_dist3 <- sum(strikes)


dist1_r <- rasterize(whale_loc_sf_dist1, x, field = 1, fun = 'count')
dist1_r01 <- climateStability::rescale0to1(dist1_r)
dist1_sf <- st_as_sf(rasterToPolygons(dist1_r01)) 

dist3_r <- rasterize(whale_loc_sf_dist3, x, field = 1, fun = 'count')
dist3_r01 <- climateStability::rescale0to1(dist3_r)
dist3_sf <- st_as_sf(rasterToPolygons(dist3_r01)) 



p_dist1 <- ggplot() +
	geom_sf(data = dist1_sf, color = "transparent", fill = "darkred", aes(alpha = layer)) +
	geom_sf(data = ship_loc_sf_dist1, fill = "black", color = "black", shape = 16) +
	geom_sf(data = whale_loc_sf_dist1[1,], fill = "black", color = "black", shape = 17, size = 2) +
	scale_alpha(name = "Relative prob.\n of surfacing",
		range = c(0.05,1)) +
	coord_sf(datum = 32608) +
	#xlim(c(435500, 442000)) +
	#ylim(c(6486000, 6490000)) +
	theme_bw()
#p_dist1


p_dist3 <- ggplot() +
	geom_sf(data = dist3_sf, color = "transparent", fill = "darkred", aes(alpha = layer)) +
	geom_sf(data = ship_loc_sf_dist3, fill = "black", color = "black", shape = 16) +
	geom_sf(data = whale_loc_sf_dist3[1,], fill = "black", color = "black", shape = 17, size = 2) +
	scale_alpha(name = "Relative prob.\n of surfacing",
		range = c(0.05,1)) +
	coord_sf(datum = 32608) +
	#xlim(c(435500, 442000)) +
	#ylim(c(6486000, 6490000)) +
	theme_bw() + 
	theme(axis.text.y = element_blank())
#p_dist3


p_dist1 + p_dist3  + plot_layout(guides = "collect")




























# Simulate over bearing index 1 
#mean(mod_dat$ship_whale_dist)
temp <- mod_dat %>% 
	dplyr::filter(bear_ind == 1) %>% 
	group_by(same_whale_ID) %>% 
	slice(1) %>% 
	as.data.frame()
mean(temp$ship_whale_dist)
mean(abs(temp$ship_whale_bear))

# Mean ship-whale bearing (absolute) for bearing index 1 = 10.3 deg
# Mean ship-whale distance for bearing index 1 = 2231 m
beta <- NISTunits::NISTdegTOradian(10.2)
radial_dist <- 2527
vert_dist <- radial_dist * sin(beta)
perp_dist <- sqrt((radial_dist^2 - vert_dist^2))


# Ship starting location
init_ship_x <- init_x - perp_dist
init_ship_y <- init_y + vert_dist

colnames(ship_init_loc_df) <- c("x", "y")
ship_loc1 <- c(init_ship_x, init_ship_y)
ship_loc2 <- c(ship_loc1[1] + (perp_dist/2), init_ship_y)
ship_loc3 <- c(ship_loc2[1] + (perp_dist/2), init_ship_y)

ship_loc_df <- as.data.frame(rbind(ship_loc1, ship_loc2, ship_loc3))
colnames(ship_loc_df) <- c("x", "y")
coordinates(ship_loc_df) <- ~x+y
ship_loc_sf <- st_as_sf(ship_loc_df)  %>%
	st_set_crs("+proj=utm +zone=8 +datum=WGS84") 	
ship_loc_sf_bear1 <- ship_loc_sf	
	
strikes <-numeric()


# WHALE
rnd_iter1 <- 1
rnd_iter2 <- 2
step1 <- sims_bear1$step[rnd_iter1]*1000
step2 <- sims_bear1$step[rnd_iter2]*1000
turn1 <- sims_bear1$turn[rnd_iter1]
turn2 <- sims_bear1$turn[rnd_iter2]

whale_loc1 <- c(init_x, init_y)
whale_move1_perp_dist <- (step1 * sin(turn1))
whale_move1_vert_dist <- sqrt(step1^2 - whale_move1_perp_dist^2)
whale_loc2 <- c(whale_loc1[1] + whale_move1_perp_dist, whale_loc1[2] + whale_move1_vert_dist)
whale_move2_perp_dist <- (step2 * sin(turn2))
whale_move2_vert_dist <- sqrt(step2^2 - whale_move2_perp_dist^2)
whale_loc3 <- c(whale_loc2[1] + whale_move2_perp_dist, whale_loc2[2] + whale_move2_vert_dist)

whale_loc_df <- as.data.frame(rbind(whale_loc1, whale_loc2, whale_loc3))
colnames(whale_loc_df) <- c("x", "y")
coordinates(whale_loc_df) <- ~x+y
whale_loc_sf <- st_as_sf(whale_loc_df)  %>%
	st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
	mutate(obs_num = 1:n()) %>%
	mutate(iter_num = 1)
	
#strike1 <- st_intersects(st_buffer(ship_loc_sf[2,], 50), st_buffer(add_sf[2,], 10), sparse = FALSE)
#strike1 <- FALSE
strike2 <- st_intersects(st_buffer(ship_loc_sf[3,], 50), st_buffer(add_sf[3,], 10), sparse = FALSE)
#strikes_tmp <- c(strike1, strike2)
strikes_tmp <- strike2
n_strike <- length(which(strikes_tmp == TRUE))

strikes[1] <- n_strike


for(i in 1:length(sim_seq)){
	rnd_iter1 <- sim_seq[i]
	rnd_iter2 <- rnd_iter1 + 1
	step1 <- sims_bear1$step[rnd_iter1]*1000
	step2 <- sims_bear1$step[rnd_iter2]*1000
	turn1 <- sims_bear1$turn[rnd_iter1]
	turn2 <- sims_bear1$turn[rnd_iter2]

	whale_loc1 <- c(init_x, init_y)
	whale_move1_perp_dist <- (step1 * sin(turn1))
	whale_move1_vert_dist <- sqrt(step1^2 - whale_move1_perp_dist^2)
	whale_loc2 <- c(whale_loc1[1] + whale_move1_perp_dist, whale_loc1[2] + whale_move1_vert_dist)
	whale_move2_perp_dist <- (step2 * sin(turn2))
	whale_move2_vert_dist <- sqrt(step2^2 - whale_move2_perp_dist^2)
	whale_loc3 <- c(whale_loc2[1] + whale_move2_perp_dist, whale_loc2[2] + whale_move2_vert_dist)

	whale_loc_df <- as.data.frame(rbind(whale_loc1, whale_loc2, whale_loc3))
	colnames(whale_loc_df) <- c("x", "y")
	coordinates(whale_loc_df) <- ~x+y
	add_sf <- st_as_sf(whale_loc_df)  %>%
		st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
		mutate(obs_num = 1:n()) %>%
		mutate(iter_num = i)
	
	#strike1 <- st_intersects(st_buffer(ship_loc_sf[2,], 50), st_buffer(add_sf[2,], 10), sparse = FALSE)
	#strike1 <- FALSE
	strike2 <- st_intersects(st_buffer(ship_loc_sf[3,], 50), st_buffer(add_sf[3,], 10), sparse = FALSE)
	#strikes_tmp <- c(strike1, strike2)
	strikes_tmp <- strike2
	n_strike <- length(which(strikes_tmp == TRUE))
	strikes[i] <- n_strike

	whale_loc_sf <- rbind(whale_loc_sf, add_sf)
	
	}

whale_loc_sf_bear1 <- whale_loc_sf 	
strk_bear1 <- sum(strikes)



# Simulate over bearing index 2 
mean(mod_dat$ship_whale_dist)
temp <- mod_dat %>% 
	dplyr::filter(bear_ind == 2) %>% 
	group_by(same_whale_ID) %>% 
	slice(1) %>% 
	as.data.frame()
mean(temp$ship_whale_dist)
mean(abs(temp$ship_whale_bear))

# Mean ship-whale bearing (absolute) for bearing index 1 = 10.3 deg
# Mean ship-whale distance for bearing index 1 = 2231 m
beta <- NISTunits::NISTdegTOradian(32.2)
radial_dist <- 1809
vert_dist <- radial_dist * sin(beta)
perp_dist <- sqrt((radial_dist^2 - vert_dist^2))


# Ship starting location
init_ship_x <- init_x - perp_dist
init_ship_y <- init_y + vert_dist

colnames(ship_init_loc_df) <- c("x", "y")
ship_loc1 <- c(init_ship_x, init_ship_y)
ship_loc2 <- c(ship_loc1[1] + (perp_dist/2), init_ship_y)
ship_loc3 <- c(ship_loc2[1] + (perp_dist/2), init_ship_y)

ship_loc_df <- as.data.frame(rbind(ship_loc1, ship_loc2, ship_loc3))
colnames(ship_loc_df) <- c("x", "y")
coordinates(ship_loc_df) <- ~x+y
ship_loc_sf <- st_as_sf(ship_loc_df)  %>%
	st_set_crs("+proj=utm +zone=8 +datum=WGS84") 	
ship_loc_sf_bear2 <- ship_loc_sf	
	
strikes <-numeric()


# WHALE
rnd_iter1 <- 1
rnd_iter2 <- 2
step1 <- sims_bear2$step[rnd_iter1]*1000
step2 <- sims_bear2$step[rnd_iter2]*1000
turn1 <- sims_bear2$turn[rnd_iter1]
turn2 <- sims_bear2$turn[rnd_iter2]

whale_loc1 <- c(init_x, init_y)
whale_move1_perp_dist <- (step1 * sin(turn1))
whale_move1_vert_dist <- sqrt(step1^2 - whale_move1_perp_dist^2)
whale_loc2 <- c(whale_loc1[1] + whale_move1_perp_dist, whale_loc1[2] + whale_move1_vert_dist)
whale_move2_perp_dist <- (step2 * sin(turn2))
whale_move2_vert_dist <- sqrt(step2^2 - whale_move2_perp_dist^2)
whale_loc3 <- c(whale_loc2[1] + whale_move2_perp_dist, whale_loc2[2] + whale_move2_vert_dist)

whale_loc_df <- as.data.frame(rbind(whale_loc1, whale_loc2, whale_loc3))
colnames(whale_loc_df) <- c("x", "y")
coordinates(whale_loc_df) <- ~x+y
whale_loc_sf <- st_as_sf(whale_loc_df)  %>%
	st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
	mutate(obs_num = 1:n()) %>%
	mutate(iter_num = 1)
	
#strike1 <- st_intersects(st_buffer(ship_loc_sf[2,], 50), st_buffer(add_sf[2,], 10), sparse = FALSE)
#strike1 <- FALSE
strike2 <- st_intersects(st_buffer(ship_loc_sf[3,], 50), st_buffer(add_sf[3,], 10), sparse = FALSE)
#strikes_tmp <- c(strike1, strike2)
strikes_tmp <- strike2
n_strike <- length(which(strikes_tmp == TRUE))

strikes[1] <- n_strike


for(i in 1:length(sim_seq)){
	rnd_iter1 <- sim_seq[i]
	rnd_iter2 <- rnd_iter1 + 1
	step1 <- sims_bear1$step[rnd_iter1]*1000
	step2 <- sims_bear1$step[rnd_iter2]*1000
	turn1 <- sims_bear1$turn[rnd_iter1]
	turn2 <- sims_bear1$turn[rnd_iter2]

	whale_loc1 <- c(init_x, init_y)
	whale_move1_perp_dist <- (step1 * sin(turn1))
	whale_move1_vert_dist <- sqrt(step1^2 - whale_move1_perp_dist^2)
	whale_loc2 <- c(whale_loc1[1] + whale_move1_perp_dist, whale_loc1[2] + whale_move1_vert_dist)
	whale_move2_perp_dist <- (step2 * sin(turn2))
	whale_move2_vert_dist <- sqrt(step2^2 - whale_move2_perp_dist^2)
	whale_loc3 <- c(whale_loc2[1] + whale_move2_perp_dist, whale_loc2[2] + whale_move2_vert_dist)

	whale_loc_df <- as.data.frame(rbind(whale_loc1, whale_loc2, whale_loc3))
	colnames(whale_loc_df) <- c("x", "y")
	coordinates(whale_loc_df) <- ~x+y
	add_sf <- st_as_sf(whale_loc_df)  %>%
		st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
		mutate(obs_num = 1:n()) %>%
		mutate(iter_num = i)
	
	#strike1 <- st_intersects(st_buffer(ship_loc_sf[2,], 50), st_buffer(add_sf[2,], 10), sparse = FALSE)
	#strike1 <- FALSE
	strike2 <- st_intersects(st_buffer(ship_loc_sf[3,], 50), st_buffer(add_sf[3,], 10), sparse = FALSE)
	#strikes_tmp <- c(strike1, strike2)
	strikes_tmp <- strike2
	n_strike <- length(which(strikes_tmp == TRUE))
	strikes[i] <- n_strike

	whale_loc_sf <- rbind(whale_loc_sf, add_sf)
	
	}

whale_loc_sf_bear2 <- whale_loc_sf 	
strk_bear2 <- sum(strikes)



# Simulate over bearing index 3
mean(mod_dat$ship_whale_dist)
temp <- mod_dat %>% 
	dplyr::filter(bear_ind == 3) %>% 
	group_by(same_whale_ID) %>% 
	slice(1) %>% 
	as.data.frame()
mean(temp$ship_whale_dist)
mean(abs(temp$ship_whale_bear))

beta <- NISTunits::NISTdegTOradian(57.6)
radial_dist <- 2022
vert_dist <- radial_dist * sin(beta)
perp_dist <- sqrt((radial_dist^2 - vert_dist^2))


# Ship starting location
init_ship_x <- init_x - perp_dist
init_ship_y <- init_y + vert_dist

colnames(ship_init_loc_df) <- c("x", "y")
ship_loc1 <- c(init_ship_x, init_ship_y)
ship_loc2 <- c(ship_loc1[1] + (perp_dist/2), init_ship_y)
ship_loc3 <- c(ship_loc2[1] + (perp_dist/2), init_ship_y)

ship_loc_df <- as.data.frame(rbind(ship_loc1, ship_loc2, ship_loc3))
colnames(ship_loc_df) <- c("x", "y")
coordinates(ship_loc_df) <- ~x+y
ship_loc_sf <- st_as_sf(ship_loc_df)  %>%
	st_set_crs("+proj=utm +zone=8 +datum=WGS84") 	
ship_loc_sf_bear3 <- ship_loc_sf	

strikes <-numeric()


# WHALE
rnd_iter1 <- 1
rnd_iter2 <- 2
step1 <- sims_bear3$step[rnd_iter1]*1000
step2 <- sims_bear3$step[rnd_iter2]*1000
turn1 <- sims_bear3$turn[rnd_iter1]
turn2 <- sims_bear3$turn[rnd_iter2]

whale_loc1 <- c(init_x, init_y)
whale_move1_perp_dist <- (step1 * sin(turn1))
whale_move1_vert_dist <- sqrt(step1^2 - whale_move1_perp_dist^2)
whale_loc2 <- c(whale_loc1[1] + whale_move1_perp_dist, whale_loc1[2] + whale_move1_vert_dist)
whale_move2_perp_dist <- (step2 * sin(turn2))
whale_move2_vert_dist <- sqrt(step2^2 - whale_move2_perp_dist^2)
whale_loc3 <- c(whale_loc2[1] + whale_move2_perp_dist, whale_loc2[2] + whale_move2_vert_dist)

whale_loc_df <- as.data.frame(rbind(whale_loc1, whale_loc2, whale_loc3))
colnames(whale_loc_df) <- c("x", "y")
coordinates(whale_loc_df) <- ~x+y
whale_loc_sf <- st_as_sf(whale_loc_df)  %>%
	st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
	mutate(obs_num = 1:n()) %>%
	mutate(iter_num = 1)
	
#strike1 <- st_intersects(st_buffer(ship_loc_sf[2,], 50), st_buffer(add_sf[2,], 10), sparse = FALSE)
#strike1 <- FALSE
strike2 <- st_intersects(st_buffer(ship_loc_sf[3,], 50), st_buffer(add_sf[3,], 10), sparse = FALSE)
#strikes_tmp <- c(strike1, strike2)
strikes_tmp <- strike2
n_strike <- length(which(strikes_tmp == TRUE))

strikes[1] <- n_strike

for(i in 1:length(sim_seq)){
	rnd_iter1 <- sim_seq[i]
	rnd_iter2 <- rnd_iter1 + 1
	
	step1 <- sims_bear3$step[rnd_iter1]*1000
	step2 <- sims_bear3$step[rnd_iter2]*1000
	turn1 <- sims_bear3$turn[rnd_iter1]
	turn2 <- sims_bear3$turn[rnd_iter2]

	whale_loc1 <- c(init_x, init_y)
	whale_move1_perp_dist <- (step1 * sin(turn1))
	whale_move1_vert_dist <- sqrt(step1^2 - whale_move1_perp_dist^2)
	whale_loc2 <- c(whale_loc1[1] + whale_move1_perp_dist, whale_loc1[2] + whale_move1_vert_dist)
	whale_move2_perp_dist <- (step2 * sin(turn2))
	whale_move2_vert_dist <- sqrt(step2^2 - whale_move2_perp_dist^2)
	whale_loc3 <- c(whale_loc2[1] + whale_move2_perp_dist, whale_loc2[2] + whale_move2_vert_dist)

	whale_loc_df <- as.data.frame(rbind(whale_loc1, whale_loc2, whale_loc3))
	colnames(whale_loc_df) <- c("x", "y")
	coordinates(whale_loc_df) <- ~x+y
	add_sf <- st_as_sf(whale_loc_df)  %>%
		st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
		mutate(obs_num = 1:n()) %>%
		mutate(iter_num = i)
	
	#strike1 <- st_intersects(st_buffer(ship_loc_sf[2,], 50), st_buffer(add_sf[2,], 10), sparse = FALSE)
	#strike1 <- FALSE
	strike2 <- st_intersects(st_buffer(ship_loc_sf[3,], 50), st_buffer(add_sf[3,], 10), sparse = FALSE)
	#strikes_tmp <- c(strike1, strike2)
	strikes_tmp <- strike2
	n_strike <- length(which(strikes_tmp == TRUE))

	strikes[i] <- n_strike

	whale_loc_sf <- rbind(whale_loc_sf, add_sf)
	
	}

whale_loc_sf_bear3 <- whale_loc_sf 	
strk_bear3 <- sum(strikes)




bear1UD <- kernelUD(as(whale_loc_sf_bear1, "Spatial"), grid = 500)
rbear1 <- raster(as(bear1UD, "SpatialPixelsDataFrame"))
#rbear1[rbear1 ==0] <- NA
rbear1_01 <- climateStability::rescale0to1(rbear1)
bear1UD_sf <- st_as_sf(rasterToPolygons(rbear1_01)) 



bear2UD <- kernelUD(as(whale_loc_sf_bear2, "Spatial"), grid = 500)
rbear2 <- raster(as(bear2UD, "SpatialPixelsDataFrame"))
rbear2[rbear2 ==0] <- NA
rbear2_01 <- climateStability::rescale0to1(rbear2)
bear2UD_sf <- st_as_sf(rasterToPolygons(rbear2_01)) 

bear3UD <- kernelUD(as(whale_loc_sf_bear3, "Spatial"), grid = 500)
rbear3 <- raster(as(bear3UD, "SpatialPixelsDataFrame"))
#rbear3[rbear3 ==0] <- NA
rbear3_01 <- climateStability::rescale0to1(rbear3)
bear3UD_sf <- st_as_sf(rasterToPolygons(rbear3_01)) 




p_bear1 <- ggplot() +
	geom_sf(data = bear1UD_sf, color = "transparent", fill = "darkred", aes(alpha = ud)) +
	geom_sf(data = ship_loc_sf_bear1, fill = "black", color = "black", shape = 16) +
	geom_sf(data = whale_loc_sf_bear1[1,], fill = "black", color = "black", shape = 17, size = 2) +
	scale_alpha(name = "Relative prob.\n of surfacing",
		range = c(0.05,1)) +
	coord_sf(datum = 32608) +
	#xlim(c(435500, 442000)) +
	#ylim(c(6486000, 6490000)) +
	theme_bw()
#p_bear1


p_bear2 <- ggplot() +
	geom_sf(data = bear2UD_sf, color = "transparent", fill = "darkred", aes(alpha = ud)) +
	geom_sf(data = ship_loc_sf_bear2, fill = "black", color = "black", shape = 16) +
	geom_sf(data = whale_loc_sf_bear2[1,], fill = "black", color = "black", shape = 17, size = 2) +
	scale_alpha(name = "Relative prob.\n of surfacing",
		range = c(0.05,1)) +
	coord_sf(datum = 32608) +
	#xlim(c(435500, 442000)) +
	#ylim(c(6486000, 6490000)) +
	theme_bw() + 
	theme(axis.text.y = element_blank())
#p_bear2


p_bear3 <- ggplot() +
	geom_sf(data = bear3UD_sf, color = "transparent", fill = "darkred", aes(alpha = ud)) +
	geom_sf(data = ship_loc_sf_bear3, fill = "black", color = "black", shape = 16) +
	geom_sf(data = whale_loc_sf_bear3[1,], fill = "black", color = "black", shape = 17, size = 2) +
	scale_alpha(name = "Relative prob.\n of surfacing",
		range = c(0.05,1)) +
	coord_sf(datum = 32608) +
	#xlim(c(435500, 442000)) +
	#ylim(c(6486000, 6490000)) +
	theme_bw() + 
	theme(axis.text.y = element_blank())
#p_bear3


p_bear1 + p_bear2 + p_bear3  + plot_layout(guides = "collect")


