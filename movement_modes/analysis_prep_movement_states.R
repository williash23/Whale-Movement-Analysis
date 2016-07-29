# Sara Williams
# 12/8/2015; updated 2/1/2016, 3/9/2016
# Data prep for model running script.
################################################################################

#  Load packages
library(dplyr)
library(adehabitatLT)    
library(tidyr)
library(stringr)
library(aspace)
################################################################################

#  Load data
dat_raw <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/Whales_0615_locations_clean.csv")
dat <- dat_raw %>%
            group_by(same_whale_ID) %>%
            filter(n() >1) %>%
            ungroup() %>%
            as.data.frame()


#  Generate step lengths and turning angles using ADEpackage
locs_a <- arrange(dat, same_whale_ID, ob_order_time)
locs_b <- locs_a %>%
                dplyr::select(same_whale_ID, unique_event_ID, X_whale_UTM, Y_whale_UTM) %>%
                dplyr::rename(X = X_whale_UTM, Y = Y_whale_UTM)
locs_b$same_whale_ID <- droplevels(locs_b$same_whale_ID)

# obs_path <- ggplot() + 
                     # geom_path(data = locs_b[1:1000,], aes(x = X, y = Y, group = same_whale_ID, colour = same_whale_ID),
                     # show.legend = FALSE) +
                     # #geom_path(data = df_XY2, aes(x = X, y = Y, group = walk_num)) +
                     # xlab("X (km)") +
                     # xlim(c(432500, 450000)) +
                     # ylab("Y (km)") +
                     # ylim(c(6460000, 6500000)) +
                     # theme_bw()
# obs_path
 
#   Create ltraj object        
whale_traj <- as.ltraj(xy = locs_b[,c("X","Y")], id = locs_b$same_whale_ID, typeII = FALSE)

#   Convert into dataframe
traj_dat <- ld(whale_traj)


#  Add location where whale was spotted (as grid ID) to data set
#gridID <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Sighting-Density/data/adjusted_density_by_first_sighting.csv")
#gridcovs <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Sighting-Density/data/env_ship_covs_sighting_density_by_gridID.csv")
#grid_dat <- left_join(gridID, gridcovs, by = "grid_ID")

#  Join covariates from grid cell to whale observations, by the grid that the whale was first sighted in.
#traj_dat <-  left_join(traj, grid_dat, by = "same_whale_ID")
# traj_dat <-  traj
################################################################################

#  Select only needed variables and rename
    # tmp <- traj_dat %>%
             # dplyr::select(same_whale_ID, x, y, dist, rel.angle, grid_ID, distance, ship_speed_scaled,
                                    # trk_length_sum_km, sst_clim, chlor_clim, bath, bath_buff_500, cell_area) 
    # names(tmp) <- c("ID", "X", "Y", "steps", "turns", "gridID",  "ship_dist", "ship_speed_scaled", 
                                  # "trk_length_sum_km", "sst_clim", "chlor_clim", "bath", "bath_buff_500", "cell_area")
tmp <- traj_dat %>%
            dplyr::select(id, x, y, dist, rel.angle, abs.angle) 
names(tmp) <- c("same_whale_ID", "X", "Y", "steps", "turns", "abs_angle")

#  Add time difference between successive observations to observations.
        time_1 <- dat %>%
                         dplyr::select(TimeTxt, X_whale_UTM, Y_whale_UTM, ob_order_time, ship_whale_dist) %>%
                         as.data.frame()
        time_2 <- as.data.frame(str_split_fixed(time_1$TimeTxt, ":", 3))
        time_3 <- cbind(time_2, time_1$X_whale_UTM, time_1$Y_whale_UTM, time_1$ob_order_time,
                                   time_1$ship_whale_dist)
        names(time_3)[4] <- "X"
        names(time_3)[5] <- "Y"
        names(time_3)[6] <- "ob_order_time"
        names(time_3)[7] <- "ship_whale_dist"
        time_3$V1 <- as.character(time_3$V1)
        time_3$V2 <- as.character(time_3$V2)
        time_3$V3 <- as.character(time_3$V3)
        time_3$V1 <- as.integer(time_3$V1)
        time_3$V2 <- as.integer(time_3$V2)
        time_3$V3 <- as.integer(time_3$V3)
        time_3$V1 <- as.numeric(time_3$V1)
        time_3$V2 <- as.numeric(time_3$V2)
        time_3$V3 <- as.numeric(time_3$V3)

        time_4 <- time_3 %>%
                         dplyr::rename (hr = V1, min = V2, sec = V3) %>%
                         dplyr::mutate(sec_frac = sec/60, min_sec = sec_frac+min) %>%
                         dplyr::mutate(min_frac = min_sec/60, time = hr+min_frac)

tmp2 <- cbind(tmp, time_4$time, time_4$ob_order_time, time_4$ship_whale_dist)
names(tmp2)[7] <- "time_hrs"
names(tmp2)[8] <- "ob_order_time"
names(tmp2)[9] <- "ship_whale_dist"
tmp3 <- tmp2 %>%
              group_by(same_whale_ID) %>%
              mutate(time_diff = lead(time_hrs) - time_hrs) %>%
              mutate(time_diff_sec = time_diff*3600) 
tmp3$time_diff_sec[tmp3$time_diff_sec == 0] <- 0.01
tmp4 <- tmp3%>%
              mutate(velocity_m_s = steps/time_diff_sec) %>%
              as.data.frame()

#  Remove data points based on obvious breaks or errors
input_dat <- filter(tmp4, steps < 10000) #, time_diff_sec < 750, time_diff_sec > 6, velocity_m_s < 50)
#### 2088 m and less is broken at the 95 percentile of step lengths
 
# #  Remove observations that do not have a turning angle
# tmp5 <- filter(tmp4, steps < 6000, time_diff_sec < 1200, time_diff_sec > 6, velocity_m_s < 50)
# input_dat <-  filter(tmp5,!is.na(turns))
#tmp3 <- filter(tmp2,!is.na(turns), !is.na(gridID))
################################################################################

#  Organize data and give sequential IDs: used in "single", "double", "double covariate" models.
obs_1 <- input_dat %>%
               group_by(same_whale_ID) %>%
               mutate(occ = 1:n()) %>%
               ungroup() %>%
               mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
               as.data.frame()

#  Data for SINGLE model only
l_single <- (obs_1$steps)/1000
theta_single <- obs_1$turns
# ship <- as.numeric(scale(obs_1$ship_dist))

#  Indexing
npts_1 <- nrow(obs_1)
ind_1 <- obs_1$ID_new
nind_1 <- length(unique(obs_1$ID_new))
nocc_1 <- obs_1 %>%
                 group_by(ID_new) %>%
                 summarise(nocc = n()) %>%
                 .$nocc

#  Make data SQAURE - For double, double cov models
l_double <- obs_1 %>%
                    mutate(steps_km = steps/1000) %>%
                    dplyr::select(ID_new, occ, steps_km) %>%
                    spread(occ, steps_km, fill = NA, convert = FALSE)

theta_double <- obs_1 %>%
                           dplyr::select(ID_new, occ, turns) %>%
                           spread(occ, turns, fill = NA, convert = FALSE)

#  Covariates
ship_dist_double<- obs_1 %>%
                                dplyr::select(ID_new, occ, ship_whale_dist) %>%
                                spread(occ, ship_whale_dist, fill = NA, convert = FALSE)

# sst_1 <- obs_1 %>%
              # dplyr::select(ID_new, occ, sst_clim) %>%
              # spread(occ, sst_clim, fill = NA, convert = FALSE)                  

# bath_ave_1 <- obs_1 %>%
                         # dplyr::select(ID_new, occ, bath_buff_500) %>%
                         # spread(occ, bath_buff_500, fill = NA, convert = FALSE) 
################################################################################

#  Behavioral assignment data subsetting for single movement mode model

#   Transit
tmp1_transit <- dat %>%
                           group_by(same_whale_ID) %>%
                           filter(n() >1) %>%
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
traj_dat_transit <- ld(whale_traj_transit)

# #  Join covariates from grid cell to whale observations, by the grid that the whale was first sighted in.
# traj_dat_transit <-  left_join(traj_transit, grid_dat, by = "same_whale_ID")
###################################################

#  Select only needed variables and rename
tmp_transit <- traj_dat_transit %>%
                        dplyr::select(id, x, y, dist, rel.angle, abs.angle) 
names(tmp_transit) <- c("same_whale_ID", "X", "Y", "steps", "turns", "abs_angle")

#  Add time difference between successive observations to observations.
        time_1_transit <- locs_a_transit %>%
                                     dplyr::select(TimeTxt, X_whale_UTM, Y_whale_UTM, ob_order_time, whale_behavior) %>%
                                     as.data.frame()
        time_2_transit <- as.data.frame(str_split_fixed(time_1_transit$TimeTxt, ":", 3))
        time_3_transit <- cbind(time_2_transit, time_1_transit$X_whale_UTM, time_1_transit$Y_whale_UTM, 
                                                time_1_transit$ob_order_time, time_1_transit$whale_behavior)
        names(time_3_transit)[4] <- "X"
        names(time_3_transit)[5] <- "Y"
        names(time_3_transit)[6] <- "ob_order_time"
        names(time_3_transit)[7] <- "whale_behavior"
        time_3_transit$V1 <- as.character(time_3_transit$V1)
        time_3_transit$V2 <- as.character(time_3_transit$V2)
        time_3_transit$V3 <- as.character(time_3_transit$V3)
        time_3_transit$V1 <- as.integer(time_3_transit$V1)
        time_3_transit$V2 <- as.integer(time_3_transit$V2)
        time_3_transit$V3 <- as.integer(time_3_transit$V3)
        time_3_transit$V1 <- as.numeric(time_3_transit$V1)
        time_3_transit$V2 <- as.numeric(time_3_transit$V2)
        time_3_transit$V3 <- as.numeric(time_3_transit$V3)

        time_4_transit <- time_3_transit %>%
                                      dplyr::rename (hr = V1, min = V2, sec = V3) %>%
                                      dplyr::mutate(sec_frac = sec/60, min_sec = sec_frac+min) %>%
                                      dplyr::mutate(min_frac = min_sec/60, time = hr+min_frac)

tmp2_transit <- cbind(tmp_transit, time_4_transit$time, time_4_transit$ob_order_time, time_4_transit$whale_behavior)
names(tmp2_transit)[7] <- "time_hrs"
names(tmp2_transit)[8] <- "ob_order_time"
names(tmp2_transit)[9] <- "whale_behavior"
tmp3_transit <- tmp2_transit %>%
                           group_by(same_whale_ID) %>%
                           dplyr::mutate(time_diff = lead(time_hrs) - time_hrs) %>%
                           dplyr::mutate(time_diff_sec = time_diff*3600) 
tmp3_transit$time_diff_sec[tmp3_transit$time_diff_sec == 0] <- 0.01
tmp4_transit <- tmp3_transit%>%
                           dplyr::mutate(velocity_m_s = steps/time_diff_sec) %>%
                           as.data.frame()

#  Remove data points based on obvious breaks or errors
input_dat_transit <- filter(tmp4_transit, steps < 10000) #, time_diff_sec < 750, time_diff_sec > 6, velocity_m_s < 50)
# #  Remove observations that do not have a turning angle
# tmp5_transit <- filter(tmp4_transit, steps < 6000, time_diff_sec < 1200, time_diff_sec > 6, velocity_m_s < 50)
# input_dat_transit <-  filter(tmp5_transit,!is.na(turns))
###################################################

#  Organize data and give sequential IDs: used in "single", "double", "double covariate" models.
obs_1_transit <- input_dat_transit %>%
                            group_by(same_whale_ID) %>%
                            dplyr::mutate(occ = 1:n()) %>%
                            ungroup() %>%
                            dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
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
l_single_transit <- (obs_1_transit$steps)/1000
theta_single_transit <- obs_1_transit$turns
###################################################

#  Stationary
tmp1_station <- dat %>%
                            group_by(same_whale_ID) %>%
                            filter(n() >1) %>%
                            filter(first(whale_behavior) == "RE-Resting" | first(whale_behavior) == "SA- Surface-active"| 
                                     first(whale_behavior) == "LF-Lunge-feed") %>%
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
traj_dat_station <- ld(whale_traj_station)

# #  Join covariates from grid cell to whale observations, by the grid that the whale was first sighted in.
# traj_dat_station <-  left_join(traj_station, grid_dat, by = "same_whale_ID")
###################################################

#  Select only needed variables and rename
tmp_station <- traj_dat_station %>%
                         dplyr::select(id, x, y, dist, rel.angle, abs.angle) 
names(tmp_station) <- c("same_whale_ID", "X", "Y", "steps", "turns", "abs_angle")

#  Add time difference between successive observations to observations.
        time_1_station <- locs_a_station %>%
                                      dplyr::select(TimeTxt, X_whale_UTM, Y_whale_UTM, ob_order_time, whale_behavior) %>%
                                      as.data.frame()
        time_2_station <- as.data.frame(str_split_fixed(time_1_station$TimeTxt, ":", 3))
        time_3_station <- cbind(time_2_station, time_1_station$X_whale_UTM, time_1_station$Y_whale_UTM, 
                                                time_1_station$ob_order_time, time_1_station$whale_behavior)
        names(time_3_station)[4] <- "X"
        names(time_3_station)[5] <- "Y"
        names(time_3_station)[6] <- "ob_order_time"
        names(time_3_station)[7] <- "whale_behavior"
        time_3_station$V1 <- as.character(time_3_station$V1)
        time_3_station$V2 <- as.character(time_3_station$V2)
        time_3_station$V3 <- as.character(time_3_station$V3)
        time_3_station$V1 <- as.integer(time_3_station$V1)
        time_3_station$V2 <- as.integer(time_3_station$V2)
        time_3_station$V3 <- as.integer(time_3_station$V3)
        time_3_station$V1 <- as.numeric(time_3_station$V1)
        time_3_station$V2 <- as.numeric(time_3_station$V2)
        time_3_station$V3 <- as.numeric(time_3_station$V3)

        time_4_station <- time_3_station %>%
                                      dplyr::rename (hr = V1, min = V2, sec = V3) %>%
                                      dplyr::mutate(sec_frac = sec/60, min_sec = sec_frac+min) %>%
                                      dplyr::mutate(min_frac = min_sec/60, time = hr+min_frac)

tmp2_station <- cbind(tmp_station, time_4_station$time, time_4_station$ob_order_time, time_4_station$whale_behavior)
names(tmp2_station)[7] <- "time_hrs"
names(tmp2_station)[8] <- "ob_order_time"
names(tmp2_station)[9] <- "whale_behavior"
tmp3_station <- tmp2_station %>%
                           group_by(same_whale_ID) %>%
                           dplyr::mutate(time_diff = lead(time_hrs) - time_hrs) %>%
                           dplyr::mutate(time_diff_sec = time_diff*3600) 
tmp3_station$time_diff_sec[tmp3_station$time_diff_sec == 0] <- 0.01
tmp4_station <- tmp3_station%>%
                           dplyr::mutate(velocity_m_s = steps/time_diff_sec) %>%
                           as.data.frame()

#  Remove data points based on obvious breaks or errors
input_dat_station <- filter(tmp4_station, steps < 10000) #, time_diff_sec < 750, time_diff_sec > 6, velocity_m_s < 50)
# #  Remove observations that do not have a turning angle
# tmp5_station <- filter(tmp4_station, steps < 6000, time_diff_sec < 1200, time_diff_sec > 6, velocity_m_s < 50)
# input_dat_station <-  filter(tmp5_station,!is.na(turns))
###################################################

#  Organize data and give sequential IDs: used in "single", "double", "double covariate" models.
obs_1_station <- input_dat_station %>%
                            group_by(same_whale_ID) %>%
                            dplyr::mutate(occ = 1:n()) %>%
                            ungroup() %>%
                            dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
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
l_single_station <- (obs_1_station$steps)/1000
theta_single_station <- obs_1_station$turns
####################################################################################################

#  Data subsetting within "cone of concern"
tmp1_cone <- dat %>%
                        group_by(same_whale_ID) %>%
                        filter(n() > 1) %>%
                        filter(first(ship_whale_bearing) <= 40 & first(ship_whale_bearing) >= -40) %>%
                        ungroup() %>%
                       as.data.frame()
locs_a_cone <- arrange(tmp1_cone, same_whale_ID, ob_order_time)
locs_b_cone <- locs_a_cone %>%
                            dplyr::select(same_whale_ID, X_whale_UTM, Y_whale_UTM) %>%
                            dplyr::rename(X = X_whale_UTM, Y = Y_whale_UTM)
locs_b_cone$same_whale_ID <- droplevels(locs_b_cone$same_whale_ID)

#   Create ltraj object
whale_traj_cone <- as.ltraj(xy = locs_b_cone[,c("X","Y")], id = locs_b_cone$same_whale_ID, typeII = FALSE)

#   Convert into dataframe
traj_dat_cone <- ld(whale_traj_cone)

# #  Join covariates from grid cell to whale observations, by the grid that the whale was first sighted in.
# traj_dat_cone <-  left_join(traj_cone, grid_dat, by = "same_whale_ID")
###################################################

#  Select only needed variables and rename
tmp_cone <- traj_dat_cone %>%
                        dplyr::select(id, x, y, dist, rel.angle, abs.angle) 
names(tmp_cone) <- c("same_whale_ID", "X", "Y", "steps", "turns", "abs_angle")

#  Add time difference between successive observations to observations.
        time_1_cone <- locs_a_cone %>%
                                     dplyr::select(TimeTxt, X_whale_UTM, Y_whale_UTM, ob_order_time) %>%
                                     as.data.frame()
        time_2_cone <- as.data.frame(str_split_fixed(time_1_cone$TimeTxt, ":", 3))
        time_3_cone <- cbind(time_2_cone, time_1_cone$X_whale_UTM, time_1_cone$Y_whale_UTM, 
                                                time_1_cone$ob_order_time)
        names(time_3_cone)[4] <- "X"
        names(time_3_cone)[5] <- "Y"
        names(time_3_cone)[6] <- "ob_order_time"
        time_3_cone$V1 <- as.character(time_3_cone$V1)
        time_3_cone$V2 <- as.character(time_3_cone$V2)
        time_3_cone$V3 <- as.character(time_3_cone$V3)
        time_3_cone$V1 <- as.integer(time_3_cone$V1)
        time_3_cone$V2 <- as.integer(time_3_cone$V2)
        time_3_cone$V3 <- as.integer(time_3_cone$V3)
        time_3_cone$V1 <- as.numeric(time_3_cone$V1)
        time_3_cone$V2 <- as.numeric(time_3_cone$V2)
        time_3_cone$V3 <- as.numeric(time_3_cone$V3)

        time_4_cone <- time_3_cone %>%
                                      dplyr::rename (hr = V1, min = V2, sec = V3) %>%
                                      dplyr::mutate(sec_frac = sec/60, min_sec = sec_frac+min) %>%
                                      dplyr::mutate(min_frac = min_sec/60, time = hr+min_frac)

tmp2_cone <- cbind(tmp_cone, time_4_cone$time, time_4_cone$ob_order_time)
names(tmp2_cone)[7] <- "time_hrs"
names(tmp2_cone)[8] <- "ob_order_time"
tmp3_cone <- tmp2_cone %>%
                           group_by(same_whale_ID) %>%
                           dplyr::mutate(time_diff = lead(time_hrs) - time_hrs) %>%
                           dplyr::mutate(time_diff_sec = time_diff*3600) 
tmp3_cone$time_diff_sec[tmp3_cone$time_diff_sec == 0] <- 0.01
tmp4_cone <- tmp3_cone%>%
                           dplyr::mutate(velocity_m_s = steps/time_diff_sec) %>%
                           as.data.frame()

#  Remove data points based on obvious breaks or errors
input_dat_cone <- filter(tmp4_cone, steps < 10000) #, time_diff_sec < 750, time_diff_sec > 6, velocity_m_s < 50)
# #  Remove observations that do not have a turning angle
# tmp5_cone <- filter(tmp4_cone, steps < 6000, time_diff_sec < 1200, time_diff_sec > 6, velocity_m_s < 50)
# input_dat_cone <-  filter(tmp5_cone,!is.na(turns))
###################################################

#  Organize data and give sequential IDs: used in "single", "double", "double covariate" models.
obs_1_cone <- input_dat_cone %>%
                            group_by(same_whale_ID) %>%
                            dplyr::mutate(occ = 1:n()) %>%
                            ungroup() %>%
                            dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
                            arrange(ID_new) %>%
                            as.data.frame()

#   Indexing
npts_1_cone <- nrow(obs_1_cone)
ind_1_cone <- obs_1_cone$ID_new
nind_1_cone <- length(unique(obs_1_cone$ID_new))
nocc_1_cone <- obs_1_cone %>%
                              group_by(ID_new) %>%
                              summarise(nocc = n()) %>%
                              .$nocc

#  Data 
l_single_cone <- (obs_1_cone$steps)/1000
theta_single_cone <- obs_1_cone$turns
################################################################################

#  Sort only ID's that have 2+ movement segments. Used in "double switch" and 
#   "double switch with covariate" models.
obs <- input_dat %>%
           group_by(same_whale_ID) %>%
           filter(n() > 1) %>%
           dplyr::mutate(occ = 1:n()) %>%
           ungroup() %>%
           dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
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
       dplyr::mutate(steps_km = steps/1000) %>%
       dplyr::select(ID_new, occ, steps_km) %>%
       spread(occ, steps_km, fill = NA, convert = FALSE) 

theta <- obs %>%
              dplyr::select(ID_new, occ, turns) %>%
              spread(occ, turns, fill = NA, convert = FALSE)

#  Currently use covariates
ship_dist <- obs %>%
                    dplyr::select(ID_new, occ, ship_whale_dist) %>%
                    spread(occ, ship_whale_dist, fill = NA, convert = FALSE) 

 # sst <- obs %>%
          # dplyr::select(ID_new, occ, sst_clim) %>%
          # spread(occ, sst_clim, fill = NA, convert = FALSE)              

# bath_ave <- obs %>%
                     # dplyr::select(ID_new, occ, bath_buff_500) %>%
                     # mutate(bath_ave_scale = as.numeric(scale(obs$bath_buff_500))) %>%
                     # dplyr::select(ID_new, occ, bath_ave_scale) %>%
                     # spread(occ, bath_ave_scale, fill = NA, convert = FALSE) 
################################################################################

# #  Using observations for step only models - do not need to remove turn NAs
# obs_step <- input_dat_steps %>%
                    # group_by(same_whale_ID) %>%
                    # mutate(occ = 1:n()) %>%
                    # ungroup() %>%
                    # mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
                    # arrange(ID_new) %>%
                    # as.data.frame()

# #  Indexing for step only models
# npts_step <- nrow(obs_step)
# ind_step <- obs_step$ID_new
# nind_step <- length(unique(obs_step$ID_new))
# nocc_step <- obs_step %>%
                      # group_by(ID_new) %>%
                      # summarise(nocc = n()) %>%
                      # .$nocc

# #  Data for single step only model
# l_single_step <-  (obs_step$steps)/1000

# #  Data for double and double cov step only model
# l_double_step <- obs_step %>%
                            # mutate(steps_km = steps/1000) %>%
                            # dplyr::select(ID_new, occ, steps_km) %>%
                            # spread(occ, steps_km, fill = NA, convert = FALSE)
################################################################################


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

