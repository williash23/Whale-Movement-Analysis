# Sara Williams
# 12/8/2015; updated 2/1/2016, 3/9/2016
# Data prep for model running script.
################################################################################

#  Load packages
library(dplyr)
library(adehabitatLT)    
library(tidyr)
################################################################################

#  Load data
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
#   Connect traj_dat to dataframe 
traj_dat_full_tmp <- cbind(traj_dat, dat)

#  Add location where whale was spotted (as grid ID) to data set
gridID <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Sighting-Density/data/adjusted_density_by_first_sighting.csv")
gridcovs <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Sighting-Density/data/env_ship_covs_sighting_density_by_gridID.csv")
grid_dat <- left_join(gridID, gridcovs, by = "grid_ID")

#  Join covariates from grid cell to whale observations, by the grid that the whale was first sighted in.
traj_dat_full <-  left_join(traj_dat_full_tmp, grid_dat, by = "same_whale_ID")

traj_dat_full_transit <- filter(traj_dat_full, whale_behavior.x == "DF-Dive-fluke-up")
traj_dat_full_station <- filter(traj_dat_full, whale_behavior.x == "LF-Lunge-feed" | subarea_1 == "RE-Resting" | subarea_1 == "SA-Surface-active")
################################################################################

#  Sort only ID's that have 1+ observations. Used only "single", "double", 
#   "double covariate" models.

#  Select only needed variables and rename
tmp_transit <- traj_dat_full_transit %>%
         dplyr::select(id, x, y, dist, rel.angle, grid_ID, whale_dist_shore_m, ship_whale_dist, ship_speed_scaled,
                                trk_length_sum_km, sst_clim, chlor_clim, bath, bath_buff_500, cell_area) 
names(tmp_transit) <- c("ID", "X", "Y", "steps", "turns", "gridID", "shore_dist", "ship_dist", "ship_speed_scaled", 
                              "trk_length_sum_km", "sst_clim", "chlor_clim", "bath", "bath_buff_500", "cell_area")
tmp2_transit <- filter(tmp_transit,!is.na(steps))
tmp3_transit <- filter(tmp2_transit,!is.na(turns), !is.na(gridID))
tmp4_transit <- filter(tmp3_transit, steps < 5000)

obs_1_transit <- tmp4_transit %>%
               group_by(ID) %>%
               mutate(occ = 1:n()) %>%
               ungroup() %>%
               mutate(ID_new = as.numeric(as.factor(as.character(ID)))) %>%
               arrange(ID_new) %>%
               as.data.frame()

#   Indexing
npts_1_transit <- nrow(obs_1_transit)
ind_1_transit <- obs_1_transit$ID_new
nind_1_transit <- length(unique(obs_1_transit$ID_new))
nocc_1_transit <- obs_1_transit %>%
                 group_by(ID_new) %>%
                 summarise(nocc = n()) %>%
                 .$nocc

#  Data for SINGLE model only
l_sing_transit <- (obs_1_transit$steps)/1000
theta_sing_transit <- obs_1_transit$turns


#  Stationary
#  Select only needed variables and rename
tmp_station <- traj_dat_full_station %>%
         dplyr::select(id, x, y, dist, rel.angle, grid_ID, whale_dist_shore_m, ship_whale_dist, ship_speed_scaled,
                                trk_length_sum_km, sst_clim, chlor_clim, bath, bath_buff_500, cell_area) 
names(tmp_station) <- c("ID", "X", "Y", "steps", "turns", "gridID", "shore_dist", "ship_dist", "ship_speed_scaled", 
                              "trk_length_sum_km", "sst_clim", "chlor_clim", "bath", "bath_buff_500", "cell_area")
tmp2_station <- filter(tmp_station,!is.na(steps))
tmp3_station <- filter(tmp2_station,!is.na(turns), !is.na(gridID))
tmp4_station <- filter(tmp3_station, steps < 5000)

obs_1_station <- tmp4_station %>%
               group_by(ID) %>%
               mutate(occ = 1:n()) %>%
               ungroup() %>%
               mutate(ID_new = as.numeric(as.factor(as.character(ID)))) %>%
               arrange(ID_new) %>%
               as.data.frame()

#   Indexing
npts_1_station <- nrow(obs_1_station)
ind_1_station <- obs_1_station$ID_new
nind_1_station <- length(unique(obs_1_station$ID_new))
nocc_1_station <- obs_1_station %>%
                 group_by(ID_new) %>%
                 summarise(nocc = n()) %>%
                 .$nocc

#  Data for SINGLE model only
l_sing_station <- (obs_1_station$steps)/1000
theta_sing_station <- obs_1_station$turns
                  
#  Load packages
library(rjags)
library(mcmcplots)
#  Load "glm" module for JAGS
load.module("glm")
################################################################################

#   MCMC settings
nc <- 3
ni <- 40000
nb <- 10000
nt <- 2
na <- 5000

#  Run "single_transit" model
#   Bundle data
jags.dat <- list(npts = npts_1_transit, l = l_sing_transit, theta = theta_sing_transit, nind = nind_1_transit)
 
#   Inits function
inits <- function(){list(v0 = runif(1, 0.01,  5), 
                                         lambda0 = runif(1, 0.01, 5), 
                                         rho0 = runif(1, 0.01, 1), 
                                         mu0 = runif(1, -3.14159265359, 3.14159265359),
                                         tau_alpha = runif(1, 0, 50))
                                         }

#   Parameters to monitor
params <- c("mean.v","mean.lambda", "mean.mu", "mean.rho", "scale", "tau_alpha")

out_single_transit <- jags.model(data = jags.dat,
                                            file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_state/models/single.txt", 
                                            inits = inits, 
                                            n.chains = nc, 
                                            n.adapt = na)

update(out_single_transit, n.iter = nb)

single_fit_transit <- coda.samples(out_single_transit,
                                              variable.names= params, 
                                              n.iter = ni, 
                                              thin = nt)

mcmcplot(single_fit_transit)
summary(single_fit_transit)
################################################################################

#  Run "single_station" model
#   Bundle data
jags.dat <- list(npts = npts_1_station, l = l_sing_station, theta = theta_sing_station, nind = nind_1_station)
 
#   Inits function
inits <- function(){list(v0 = runif(1, 0.01,  5), 
                                         lambda0 = runif(1, 0.01, 5), 
                                         rho0 = runif(1, 0.01, 1), 
                                         mu0 = runif(1, -3.14159265359, 3.14159265359),
                                         tau_alpha = runif(1, 0, 50))
                                         }

#   Parameters to monitor
params <- c("mean.v","mean.lambda", "mean.mu", "mean.rho", "scale", "tau_alpha")

out_single_station <- jags.model(data = jags.dat,
                                            file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_state/models/single.txt", 
                                            inits = inits, 
                                            n.chains = nc, 
                                            n.adapt = na)

update(out_single_station, n.iter = nb)

single_fit_station <- coda.samples(out_single_station,
                                              variable.names= params, 
                                              n.iter = ni, 
                                              thin = nt)

mcmcplot(single_fit_station)
summary(single_fit_station)
################################################################################



traj_dat_full_bay <- filter(traj_dat_full, subarea_1 == "Lower Bay" | subarea_1 == "Middle Bay" | subarea_1 == "West Arm")
traj_dat_full_strait <- filter(traj_dat_full, subarea_1 == "Chatham" | subarea_1 == "East Icy" | subarea_1 == "West Icy")