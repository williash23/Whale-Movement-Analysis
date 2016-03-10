# Sara Williams
# 12/8/2015; updated 2/1/2016, 3/9/2016
# Data prep for model running script.
################################################################################

#  Load packages
library(dplyr)
library(adehabitatLT)    
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
################################################################################

#  Sort only ID's that have 1+ observations. Used only "single", "double", 
#   "double covariate" models.

#  Select only needed variables and rename
tmp <- traj_dat_full %>%
         dplyr::select(id, x, y, dist, rel.angle, grid_ID, whale_dist_shore_m, ship_whale_dist, ship_speed_scaled,
                                trk_length_sum_km, sst_clim, chlor_clim, bath, bath_buff_500, cell_area) 
names(tmp) <- c("ID", "X", "Y", "steps", "turns", "gridID", "shore_dist", "ship_dist", "ship_speed_scaled", 
                              "trk_length_sum_km", "sst_clim", "chlor_clim", "bath", "bath_buff_500", "cell_area")
tmp2 <- filter(tmp,!is.na(steps))
tmp3 <- filter(tmp2,!is.na(turns), !is.na(gridID))
tmp4 <- filter(tmp3, steps < 5000)
tmp5 <- tmp4 %>%
               group_by(ID) %>%
               mutate(same_ID_indicator = ifelse(row_number() == 1, 1,0)) %>%
               as.data.frame()

#  Generate ordered ID variable
raw_ID_1 <- as.numeric(tmp4$ID)
ID_1 <- as.factor(raw_ID_1)
ID_new_1 <- as.numeric(ID_1)
tmp6 <- cbind(ID_new_1, tmp5)
obs_1 <- tmp6 %>%
                arrange(ID_new_1) %>%
                as.data.frame()

#   Indexing
npts_1 <- nrow(obs_1)
ID_1 <- obs_1$ID_new_1
same_1 <- obs_1$same_ID_indicator
nind_1 <- length(unique(ID_1))
#   Steps and turns
l_1 <-as.numeric((obs_1$steps)/1000)
theta_1 <- as.numeric(obs_1$turns)
#   Covariates
shore_dist_1 <-as.numeric(obs_1$shore_dist)
ship_dist_1 <- as.numeric(obs_1$ship_dist)
ship_speed_1 <-as.numeric(obs_1$ship_speed_scaled)
ship_dens_1 <- as.numeric(obs_1$trk_length_sum_km)
chlor_1 <- as.numeric(obs_1$chlor_clim)
sst_1 <- as.numeric(obs_1$sst_clim)
bath_1 <- as.numeric(obs_1$bath)
bath_ave_1 <- as.numeric(obs_1$bath_buff_500)
#   Scaled and centered covariates
# shore_dist_1 <-as.numeric(scale(obs_1$shore_dist))
# ship_dist_1 <- as.numeric(scale(obs_1$ship_dist))
# ship_speed_1 <-as.numeric(scale(obs_1$ship_speed_scaled))
# ship_dens_1 <- as.numeric(scale(obs_1$trk_length_sum_km))
# chlor_1 <- as.numeric(scale(obs_1$chlor_clim)) 
# sst_1 <- as.numeric(scale(obs_1$sst_clim)) 
# bath_1 <- as.numeric(scale(obs_1$bath)) 
# bath_ave_1 <- as.numeric(scale(obs_1$bath_buff_500)) 
################################################################################

#  Sort only ID's that have 2+ observations. Used in "double switch" and 
#   "double switch with covariate" models.
tmp7 <- tmp4 %>%
               group_by(ID) %>%
               filter(n() > 1) %>%
               mutate(same_ID_indicator = ifelse(row_number() == 1, 1,0)) %>%
               as.data.frame()

#  Generate ordered ID variable
raw_ID <- as.numeric(tmp7$ID)
ID <- as.factor(raw_ID)
ID_new <- as.numeric(ID)
tmp8 <- cbind(ID_new, tmp7)
obs <- tmp8 %>%
            arrange(ID_new) %>%
            as.data.frame()

 #   Indexing
npts <- nrow(obs)
ID <- obs$ID_new
same <- obs$same_ID_indicator
nind <- length(unique(ID))
#   Steps and turns
l <- as.numeric((obs$steps)/1000)
theta <- as.numeric(obs$turns)
#   Covariates
shore_dist <-as.numeric(obs$shore_dist)
ship_dist <- as.numeric(obs$ship_dist)
ship_speed <-as.numeric(obs$ship_speed_scaled)
ship_dens <- as.numeric(obs$trk_length_sum_km)
chlor <- as.numeric(obs$chlor_clim)
sst <- as.numeric(obs$sst_clim)
bath <- as.numeric(obs$bath) 
bath_ave <- as.numeric(obs$bath_buff_500)
#   Scaled and centered covariates
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

