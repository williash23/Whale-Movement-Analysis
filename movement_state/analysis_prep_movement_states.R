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
locs_a <- arrange(dat, same_whale_ID, ob_order_time)
locs_b <- locs_a %>%
                dplyr::select(same_whale_ID, X_whale_UTM, Y_whale_UTM) %>%
                dplyr::rename(X = X_whale_UTM, Y = Y_whale_UTM)
locs_b$same_whale_ID <- droplevels(locs_b$same_whale_ID)

#   Create ltraj object        
whale_traj <- as.ltraj(xy = locs_b[,c("X","Y")], id = locs_b$same_whale_ID, typeII = FALSE)

#   Convert into dataframe
traj <- ld(whale_traj)
traj <- traj %>%
           dplyr::rename(same_whale_ID = id)

#  Add location where whale was spotted (as grid ID) to data set
gridID <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Sighting-Density/data/adjusted_density_by_first_sighting.csv")
gridcovs <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Sighting-Density/data/env_ship_covs_sighting_density_by_gridID.csv")
grid_dat <- left_join(gridID, gridcovs, by = "grid_ID")

#  Join covariates from grid cell to whale observations, by the grid that the whale was first sighted in.
traj_dat <-  left_join(traj, grid_dat, by = "same_whale_ID")
################################################################################

#  Sort only ID's that have 1+ observations. Used only "single", "double", 
#   "double covariate" models.

#  Select only needed variables and rename
tmp <- traj_dat %>%
         dplyr::select(same_whale_ID, x, y, dist, rel.angle, grid_ID, distance, ship_speed_scaled,
                                trk_length_sum_km, sst_clim, chlor_clim, bath, bath_buff_500, cell_area) 
names(tmp) <- c("ID", "X", "Y", "steps", "turns", "gridID",  "ship_dist", "ship_speed_scaled", 
                              "trk_length_sum_km", "sst_clim", "chlor_clim", "bath", "bath_buff_500", "cell_area")
tmp2 <- filter(tmp,!is.na(steps))
tmp3 <- filter(tmp2,!is.na(turns), !is.na(gridID))
tmp4 <- filter(tmp3, steps < 5000)

obs_1 <- tmp4 %>%
               group_by(ID) %>%
               mutate(occ = 1:n()) %>%
               ungroup() %>%
               mutate(ID_new = as.numeric(as.factor(as.character(ID)))) %>%
               arrange(ID_new) %>%
               as.data.frame()

#  Data for SINGLE model only
l_sing <- (obs_1$steps)/1000
theta_sing <- obs_1$turns

#   Indexing
npts_1 <- nrow(obs_1)
ind_1 <- obs_1$ID_new
nind_1 <- length(unique(obs_1$ID_new))
nocc_1 <- obs_1 %>%
                 group_by(ID_new) %>%
                 summarise(nocc = n()) %>%
                 .$nocc

#   Steps and turns - MAKE SQAURE!!! - For double, double cov models
l_1 <- obs_1 %>%
           mutate(steps_km = steps/1000) %>%
           dplyr::select(ID_new, occ, steps_km) %>%
           spread(occ, steps_km, fill = NA, convert = FALSE)

theta_1 <- obs_1 %>%
                  dplyr::select(ID_new, occ, turns) %>%
                  spread(occ, turns, fill = NA, convert = FALSE)

#  Covariates
ship_dist_1<- obs_1 %>%
                  dplyr::select(ID_new, occ, ship_dist) %>%
                  spread(occ, ship_dist, fill = NA, convert = FALSE)

sst_1 <- obs_1 %>%
              dplyr::select(ID_new, occ, sst_clim) %>%
              spread(occ, sst_clim, fill = NA, convert = FALSE)                  

bath_ave_1 <- obs_1 %>%
                         dplyr::select(ID_new, occ, bath_buff_500) %>%
                         spread(occ, bath_buff_500, fill = NA, convert = FALSE) 
################################################################################

#  Sort only ID's that have 2+ movement segments. Used in "double switch" and 
#   "double switch with covariate" models.
obs <- tmp4 %>%
           group_by(ID) %>%
           filter(n() > 1) %>%
           mutate(occ = 1:n()) %>%
           ungroup() %>%
           mutate(ID_new = as.numeric(as.factor(as.character(ID)))) %>%
           arrange(ID_new) %>%
           as.data.frame()

#   Indexing
npts <- nrow(obs)
ind <- obs$ID_new
nind <- length(unique(obs$ID_new))
nocc <- obs %>%
             group_by(ID_new) %>%
             summarise(nocc = n()) %>%
             .$nocc

#   Steps and turns - MAKE SQAURE!!!
l <- obs %>%
       mutate(steps_km = steps/1000) %>%
       dplyr::select(ID_new, occ, steps_km) %>%
       spread(occ, steps_km, fill = NA, convert = FALSE) 

theta <- obs %>%
              dplyr::select(ID_new, occ, turns) %>%
              spread(occ, turns, fill = NA, convert = FALSE)

#  Currently use covariates
ship_dist <- obs %>%
                    dplyr::select(ID_new, occ, ship_dist) %>%
                    spread(occ, ship_dist, fill = NA, convert = FALSE) 

 sst <- obs %>%
          dplyr::select(ID_new, occ, sst_clim) %>%
          spread(occ, sst_clim, fill = NA, convert = FALSE)              

bath_ave <- obs %>%
                     dplyr::select(ID_new, occ, bath_buff_500) %>%
                     mutate(bath_ave_scale = as.numeric(scale(obs$bath_buff_500))) %>%
                     dplyr::select(ID_new, occ, bath_ave_scale) %>%
                     spread(occ, bath_ave_scale, fill = NA, convert = FALSE) 
####################################################################################################

#  Behavioral assignment data subsetting for single movement mode model
#   Transit
tmp1_transit <- dat %>%
                           group_by(same_whale_ID) %>%
                           filter(last(whale_behavior)  == "DF-Dive-fluke-up") %>%
                           ungroup() %>%
                           as.data.frame()
 tmp2_transit <- tmp1_transit %>%
                            filter(whale_behavior  == "DF-Dive-fluke-up" | whale_behavior  == "DN-Dive-no-fluke" |
                            whale_behavior  == "BL-Blowing")
locs_a_transit <- arrange(tmp2_transit, same_whale_ID, ob_order_time)
locs_b_transit <- locs_a_transit %>%
                            dplyr::select(same_whale_ID, X_whale_UTM, Y_whale_UTM) %>%
                            dplyr::rename(X = X_whale_UTM, Y = Y_whale_UTM)
locs_b_transit$same_whale_ID <- droplevels(locs_b_transit$same_whale_ID)

#   Create ltraj object        
whale_traj_transit <- as.ltraj(xy = locs_b_transit[,c("X","Y")], id = locs_b_transit$same_whale_ID, typeII = FALSE)

#   Convert into dataframe
traj_transit <- ld(whale_traj_transit)
traj_transit <- traj_transit %>%
                        dplyr::rename(same_whale_ID = id)

#  Join covariates from grid cell to whale observations, by the grid that the whale was first sighted in.
traj_dat_transit <-  left_join(traj_transit, grid_dat, by = "same_whale_ID")

#  Remove erroneous steps and NAs
transit_dat <- traj_dat_transit %>%
                         dplyr::select(same_whale_ID, x, y, dist, rel.angle, grid_ID, distance, ship_speed_scaled,
                         trk_length_sum_km, sst_clim, chlor_clim, bath, bath_buff_500, cell_area) 
names(transit_dat) <- c("ID", "X", "Y", "steps", "turns", "gridID",  "ship_dist", "ship_speed_scaled", 
                                          "trk_length_sum_km", "sst_clim", "chlor_clim", "bath", "bath_buff_500", "cell_area")
transit_dat2 <- filter(transit_dat,!is.na(steps))
transit_dat3 <- filter(transit_dat2,!is.na(turns), !is.na(gridID))
transit_dat4 <- filter(transit_dat3, steps < 5000)

obs_1_transit <- transit_dat4 %>%
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

#  Data 
l_sing_transit <- (obs_1_transit$steps)/1000
theta_sing_transit <- obs_1_transit$turns

#  Stationary
 tmp1_station <- dat %>%
                            group_by(same_whale_ID) %>%
                            filter(last(whale_behavior) == "RE-Resting" | last(whale_behavior) == "SA- Surface-active"| 
                                     last(whale_behavior) == "LF-Lunge-feed") %>%
                            ungroup() %>%
                            as.data.frame()
locs_a_station <- arrange(tmp1_station, same_whale_ID, ob_order_time)
locs_b_station <- locs_a_station %>%
                            dplyr::select(same_whale_ID, X_whale_UTM, Y_whale_UTM) %>%
                            dplyr::rename(X = X_whale_UTM, Y = Y_whale_UTM)
locs_b_station$same_whale_ID <- droplevels(locs_b_station$same_whale_ID)

#   Create ltraj object        
whale_traj_station <- as.ltraj(xy = locs_b_station[,c("X","Y")], id = locs_b_station$same_whale_ID, typeII = FALSE)

#   Convert into dataframe
traj_station <- ld(whale_traj_station)
traj_station <- traj_station %>%
                        dplyr::rename(same_whale_ID = id)

#  Join covariates from grid cell to whale observations, by the grid that the whale was first sighted in.
traj_dat_station <-  left_join(traj_station, grid_dat, by = "same_whale_ID")

#  Remove erroneous steps and NAs
station_dat <- traj_dat_station %>%
                         dplyr::select(same_whale_ID, x, y, dist, rel.angle, grid_ID, distance, ship_speed_scaled,
                         trk_length_sum_km, sst_clim, chlor_clim, bath, bath_buff_500, cell_area) 
names(station_dat) <- c("ID", "X", "Y", "steps", "turns", "gridID",  "ship_dist", "ship_speed_scaled", 
                                          "trk_length_sum_km", "sst_clim", "chlor_clim", "bath", "bath_buff_500", "cell_area")
station_dat2 <- filter(station_dat,!is.na(steps))
station_dat3 <- filter(station_dat2,!is.na(turns), !is.na(gridID))
station_dat4 <- filter(station_dat3, steps < 5000)

obs_1_station <- station_dat4 %>%
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

#  Data
l_sing_station <- (obs_1_station$steps)/1000
theta_sing_station <- obs_1_station$turns
####################################################################################################

dat_cone <- dat %>%
                     filter(ship_whale_bearing <= 40 & ship_whale_bearing >= -40) %>%
                     as.data.frame()

#  Generate step lengths and turning angles using ADEpackage
locs_a_cone <- arrange(dat_cone, same_whale_ID, ob_order_time)
locs_b_cone <- locs_a_cone %>%
                          dplyr::select(same_whale_ID, X_whale_UTM, Y_whale_UTM) %>%
                          dplyr::rename(X = X_whale_UTM, Y = Y_whale_UTM)
locs_b_cone$same_whale_ID <- droplevels(locs_b_cone$same_whale_ID)

#   Create ltraj object        
whale_traj_cone <- as.ltraj(xy = locs_b_cone[,c("X","Y")], id = locs_b_cone$same_whale_ID, typeII = FALSE)

#   Convert into dataframe
traj_cone <- ld(whale_traj_cone)
traj_cone <- traj_cone %>%
                        dplyr::rename(same_whale_ID = id)

#  Join covariates from grid cell to whale observations, by the grid that the whale was first sighted in.
traj_dat_cone <-  left_join(traj_cone, grid_dat, by = "same_whale_ID")
################################################################################

#  Sort only ID's that have 1+ observations. Used only "single", "double", 
#   "double covariate" models.

#  Select only needed variables and rename
tmp_cone <- traj_dat_cone %>%
                      dplyr::select(same_whale_ID, x, y, dist, rel.angle, grid_ID, distance, ship_speed_scaled,
                                            trk_length_sum_km, sst_clim, chlor_clim, bath, bath_buff_500, cell_area) 
names(tmp_cone) <- c("ID", "X", "Y", "steps", "turns", "gridID",  "ship_dist", "ship_speed_scaled", 
                                        "trk_length_sum_km", "sst_clim", "chlor_clim", "bath", "bath_buff_500", "cell_area")
tmp2_cone <- filter(tmp_cone,!is.na(steps))
tmp3_cone <- filter(tmp2_cone,!is.na(turns), !is.na(gridID))
tmp4_cone <- filter(tmp3_cone, steps < 5000)

obs_1_cone <- tmp4_cone %>%
                         group_by(ID) %>%
                         mutate(occ = 1:n()) %>%
                         ungroup() %>%
                         mutate(ID_new = as.numeric(as.factor(as.character(ID)))) %>%
                         arrange(ID_new) %>%
                         as.data.frame()

#  Data for SINGLE model only
l_sing_cone <- (obs_1_cone$steps)/1000
theta_sing_cone <- obs_1_cone$turns

#   Indexing
npts_1_cone <- nrow(obs_1_cone)
ind_1_cone <- obs_1_cone$ID_new
nind_1_cone <- length(unique(obs_1_cone$ID_new))
nocc_1_cone <- obs_1_cone %>%
                           group_by(ID_new) %>%
                           summarise(nocc = n()) %>%
                           .$nocc





#  Other avilable covariates
# shore_dist_1 <-as.numeric(obs_1$shore_dist)
# ship_dist_1 <- as.numeric(obs_1$ship_dist)
# ship_speed_1 <-as.numeric(obs_1$ship_speed_scaled)
# ship_dens_1 <- as.numeric(obs_1$trk_length_sum_km)
# chlor_1 <- as.numeric(obs_1$chlor_clim)
# sst_1 <- as.numeric(obs_1$sst_clim)
# bath_1 <- as.numeric(obs_1$bath)
# bath_ave_1 <- as.numeric(obs_1$bath_buff_500)

# shore_dist <-as.numeric(obs$shore_dist)
# ship_dist <- as.numeric(obs$ship_dist)
# ship_speed <-as.numeric(obs$ship_speed_scaled)
# ship_dens <- as.numeric(obs$trk_length_sum_km)
# chlor <- as.numeric(obs$chlor_clim)
# sst <- as.numeric(obs$sst_clim)
# bath <- as.numeric(obs$bath) 
# bath_ave <- as.numeric(obs$bath_buff_500)

#   Scaled and centered
# shore_dist_1 <-as.numeric(scale(obs_1$shore_dist))
# ship_dist_1 <- as.numeric(scale(obs_1$ship_dist))
# ship_speed_1 <-as.numeric(scale(obs_1$ship_speed_scaled))
# ship_dens_1 <- as.numeric(scale(obs_1$trk_length_sum_km))
# chlor_1 <- as.numeric(scale(obs_1$chlor_clim)) 
# sst_1 <- as.numeric(scale(obs_1$sst_clim)) 
# bath_1 <- as.numeric(scale(obs_1$bath)) 
# bath_ave_1 <- as.numeric(scale(obs_1$bath_buff_500)) 

# shore_dist <-as.numeric(scale(obs$shore_dist))
# ship_dist <- as.numeric(scale(obs$ship_dist))
# ship_speed <-as.numeric(scale(obs$ship_speed_scaled))
# ship_dens <- as.numeric(scale(obs$trk_length_sum_km))
# chlor <- as.numeric(scale(obs$chlor_clim)) 
# sst <- as.numeric(scale(obs$sst_clim)) 
# bath <- as.numeric(scale(obs$bath)) 
# bath_ave <- as.numeric(scale(obs$bath_buff_500)) 
####################################################################################################

#  JAGS parameterization of weibull distribution.
#   dweib(v = shape, lambda)
#  R parameterization of weibull distribution.
#   rweibull(n, shape, scale)
#  Transform lambda to R scale parameter.
#    scale[1] <- (1/lambda[1])^ (1/v[1])
#  Find mean step length.
#   mean <- scale*gamma((1+1/shape))

