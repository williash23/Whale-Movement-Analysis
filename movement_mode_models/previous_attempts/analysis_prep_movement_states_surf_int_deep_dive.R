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
            filter(count == 1) %>%
            filter(year > 2009) %>% ## per info from Karin, PTB (bearing) was kind of weird before 2010
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
            dplyr::select(id, x, y, dist, rel.angle, abs.angle, pkey) 
names(tmp) <- c("same_whale_ID", "X", "Y", "steps", "turns", "abs_angle", "pkey")

#  Prep data - this is what I used for my original model runs.
#  Add time difference between successive observations to observations.
        time_1 <- dat %>%
                         dplyr::select(TimeTxt, X_whale_UTM, Y_whale_UTM, ob_order_time, 
                         ship_whale_dist, ship_whale_bearing, whale_behavior) %>%
                         as.data.frame()
        time_2 <- as.data.frame(str_split_fixed(time_1$TimeTxt, ":", 3))
        time_3 <- cbind(time_2, time_1$X_whale_UTM, time_1$Y_whale_UTM, time_1$ob_order_time,
                                   time_1$ship_whale_dist, time_1$ship_whale_bearing, time_1$whale_behavior, 
                                   time_1$TimeTxt)
        names(time_3)[4] <- "X"
        names(time_3)[5] <- "Y"
        names(time_3)[6] <- "ob_order_time"
        names(time_3)[7] <- "ship_whale_dist"
        names(time_3)[8] <- "ship_whale_bearing"
        names(time_3)[9] <- "whale_behavior"
        names(time_3)[10] <- "time_txt"
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

tmp2 <- cbind(tmp, time_4$time, time_4$ob_order_time, time_4$ship_whale_dist, 
                        time_4$ship_whale_bearing, time_4$whale_behavior, time_4$time_txt)
names(tmp2)[8] <- "time_hrs"
names(tmp2)[9] <- "ob_order_time"
names(tmp2)[10] <- "ship_whale_dist"
names(tmp2)[11] <- "ship_whale_bearing"
names(tmp2)[12] <- "whale_behavior"
names(tmp2)[13] <- "time_txt"

#  Prep data - this is for new sepeartions based on deep dive vs surface interval!!!!!!!!!!!!!!!
input_dat <- tmp2 %>%
              group_by(same_whale_ID) %>%
              mutate(time_diff = lead(time_hrs) - time_hrs) %>%
              mutate(time_diff_sec = time_diff*3600) %>%
              ungroup(.)%>%
              mutate(type1 = ifelse(time_diff_sec > 180, "deep_dive", "")) %>%
              mutate(type2 = ifelse(time_diff_sec < 30, "surface_interval", "")) %>%
              mutate(type = paste(type1, type2)) %>%
              dplyr::select(steps, turns, time_diff_sec, type, same_whale_ID, X, Y, time_txt, ob_order_time, 
              ship_whale_dist, ship_whale_bearing, whale_behavior) %>%
              filter(steps < 6001) %>%
              as.data.frame()

input_dat$type[input_dat$type==" deep_dive"] <- "deep_dive" 
input_dat$type[input_dat$type=="deep_dive "] <- "deep_dive"
input_dat$type[input_dat$type==" surface_interval"] <- "surface_interval"
input_dat$type[input_dat$type=="surface_interval "] <- "surface_interval"
input_dat$type <- as.factor(input_dat$type)
sep_input_dat <- filter(input_dat, type == "deep_dive" | type == "surface_interval")

#tst <- filter(sep_input_dat, type == "deep_dive" & whale_behavior != "DF-Dive-fluke-up")

#  Select only steps where time difference is small and at natural break - would indicate that whale did not 
#  initiate a deep dive
surf_int <- filter(sep_input_dat, type=="surface_interval")
#  Select only steps where time difference is great - indicating initiaion of deep dive
deep_dive_obs <- filter(input_dat, type=="deep_dive") 
#  Find whale IDs that have a deep dive
deep_dive_whales <- deep_dive_obs %>%
                                    dplyr::select(same_whale_ID) %>%
                                    distinct()
#  Get entire group of sightings for those whales
deep_dive_full_path <- semi_join(input_dat, deep_dive_whales, by = "same_whale_ID")
deep_dive_full_path_long <- deep_dive_full_path %>%
                                                group_by(same_whale_ID) %>%
                                                filter(time_diff_sec > 180) %>%
                                                as.data.frame()
#  Plot for an idea of what we're looking at...
surf_steps_tmp <- surf_int %>% 
                       dplyr::select(steps) %>%
                       filter(steps < 6001)
surf_steps <- surf_steps_tmp %>%
                       mutate(iter_num = 1:nrow(surf_steps_tmp)) %>%
                       mutate(type = "surfacing interval steps")
dive_steps_tmp <- deep_dive_full_path_long %>% 
                               dplyr::select(steps) %>%
                               filter(steps < 6001)
dive_steps <- dive_steps_tmp %>%
                       mutate(iter_num = 1:nrow(dive_steps_tmp)) %>%
                       mutate(type = "deep dive steps")
both <- rbind(surf_steps, dive_steps)
steps_both_plot <- ggplot(both, aes(steps, colour = type, fill = type)) + 
                                geom_density(alpha= 0.7) +
                                xlab("Step length (km)") +
                                ylab("Frequency") +
                                xlim(c(0, 6000)) +
                                ylim(c(0, 0.005)) +
                                theme_bw() +
                                theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
steps_both_plot

surf_turns_tmp <- surf_int %>% 
                               dplyr::select(turns, steps) %>%
                               filter(steps < 6001)
surf_turns <- surf_turns_tmp %>%
                       dplyr::select(turns) %>%
                       mutate(iter_num = 1:nrow(surf_turns_tmp)) %>%
                       mutate(type = "surfacing interval turns")
surf_turns$turns_deg <- ((surf_turns$turns)*180)/3.14159265359
dive_turns_tmp <- deep_dive_full_path_long %>% 
                               dplyr::select(turns, steps) %>%
                               filter(steps < 6001)
dive_turns <- dive_turns_tmp %>%
                       dplyr::select(turns) %>%
                       mutate(iter_num = 1:nrow(dive_turns_tmp)) %>%
                       mutate(type = "deep dive turns")
dive_turns$turns_deg <- ((dive_turns$turns)*180)/3.14159265359
both <- rbind(surf_turns, dive_turns)
turns_both_plot <- ggplot(both, aes(turns_deg, colour = type, fill = type)) + 
                                geom_density(alpha= 0.7) +
                                xlab("Turn angle (deg)") +
                                ylab("Frequency") +
                                #xlim(c(-3.15, 3.15)) +
                                #ylim(c(0, 0.003)) +
                                theme_bw() +
                                theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
turns_both_plot
################################################################################


#  Generate surface interval only data and deep dive only data for models
#  Organize data and give sequential IDs: 

#   Surface interval
obs_1_surf <- surf_int %>%
                       filter(steps < 6001) %>%
                       dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
                       arrange(ID_new) %>%
                       as.data.frame()

#   Indexing
npts_1_surf <- nrow(obs_1_surf)
ind_1_surf <- obs_1_surf$ID_new
nind_1_surf <- length(unique(obs_1_surf$ID_new))
nocc_1_surf <- obs_1_surf %>%
                              group_by(ID_new) %>%
                              summarise(nocc = n()) %>%
                              .$nocc

#  Data 
l_single_surf <- obs_1_surf$steps #/1000
theta_single_surf <- obs_1_surf$turns


#  Deep dive
obs_1_dive <- deep_dive_full_path_long %>%
                        filter(steps < 6001) %>%
                        dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
                        arrange(ID_new) %>%
                        as.data.frame()

#   Indexing
npts_1_dive <- nrow(obs_1_dive)
ind_1_dive <- obs_1_dive$ID_new
nind_1_dive <- length(unique(obs_1_dive$ID_new))
nocc_1_dive <- obs_1_dive %>%
                              group_by(ID_new) %>%
                              summarise(nocc = n()) %>%
                              .$nocc

#  Data 
l_single_dive <- obs_1_dive$steps #/1000
theta_single_dive <- obs_1_dive$turns
################################################################################




#  Select only steps where whale was first sighted close by (<1000m)
close_obs <- sep_input_dat %>%
                      filter(steps < 6001) %>%
                      filter(ship_whale_dist < 1000)
#  Select only steps where whale was first sighted very far away (> 3000m)
far_obs <- sep_input_dat %>%
                  filter(steps < 6001) %>%
                  filter(ship_whale_dist > 3000)
#  Plot for an idea of what we're looking at...
close_obs_steps_tmp <- close_obs %>% 
                                        dplyr::select(steps) %>%
                                        filter(steps < 6001)
close_obs_steps <- close_obs_steps_tmp %>%
                                mutate(iter_num = 1:nrow(close_obs_steps_tmp)) %>%
                                mutate(type = "first sighting close steps")
far_obs_steps_tmp <- far_obs %>% 
                                     dplyr::select(steps) %>%
                                     filter(steps < 6001)
far_obs_steps <- far_obs_steps_tmp %>%
                            mutate(iter_num = 1:nrow(far_obs_steps_tmp)) %>%
                            mutate(type = "first sighting far steps")
both <- rbind(close_obs_steps, far_obs_steps)
steps_both_plot <- ggplot(both, aes(steps, colour = type, fill = type)) + 
                                geom_density(alpha= 0.7) +
                                xlab("Step length (km)") +
                                ylab("Frequency") +
                                xlim(c(0, 6000)) +
                                ylim(c(0, 0.006)) +
                                theme_bw() +
                                theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
steps_both_plot

close_obs_turns_tmp <- close_obs %>% 
                                        dplyr::select(turns, steps) %>%
                                        filter(steps < 6001)
close_obs_turns <- close_obs_turns_tmp %>%
                               dplyr::select(turns) %>%
                               mutate(iter_num = 1:nrow(close_obs_turns_tmp)) %>%
                               mutate(type = "first sighting close turns")
close_obs_turns$turns_deg <- ((close_obs_turns$turns)*180)/3.14159265359
far_obs_turns_tmp <- far_obs %>% 
                                     dplyr::select(turns, steps) %>%
                                     filter(steps < 6001)
far_obs_turns <- far_obs_turns_tmp %>%
                            dplyr::select(turns) %>%
                            mutate(iter_num = 1:nrow(far_obs_turns_tmp)) %>%
                            mutate(type = "first sighting far turns")
far_obs_turns$turns_deg <- ((far_obs_turns$turns)*180)/3.14159265359
both <- rbind(close_obs_turns, far_obs_turns)
turns_both_plot <- ggplot(both, aes(turns_deg, colour = type, fill = type)) + 
                                geom_density(alpha= 0.7) +
                                xlab("Turn angle (deg)") +
                                ylab("Frequency") +
                                #xlim(c(-3.15, 3.15)) +
                                #ylim(c(0, 0.003)) +
                                theme_bw() +
                                theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
turns_both_plot
################################################################################


#  Generate first sighting close only data and first sighting far only data for models
#  Organize data and give sequential IDs: 

#   First sighting close (<1000m)
obs_1_close <- close_obs %>%
                         dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
                         arrange(ID_new) %>%
                         as.data.frame()

#   Indexing
npts_1_close <- nrow(obs_1_close)
ind_1_close <- obs_1_close$ID_new
nind_1_close <- length(unique(obs_1_close$ID_new))
nocc_1_close <- obs_1_close %>%
                            group_by(ID_new) %>%
                            summarise(nocc = n()) %>%
                            .$nocc

#  Data 
l_single_close <- obs_1_close$steps #/1000
theta_single_close <- obs_1_close$turns


#  First sighting far ( > 3000m)
obs_1_far <- far_obs %>%
                      dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
                      arrange(ID_new) %>%
                      as.data.frame()

#   Indexing
npts_1_far <- nrow(obs_1_far)
ind_1_far <- obs_1_far$ID_new
nind_1_far <- length(unique(obs_1_far$ID_new))
nocc_1_far <- obs_1_far %>%
                        group_by(ID_new) %>%
                        summarise(nocc = n()) %>%
                        .$nocc

#  Data 
l_single_far <- obs_1_far$steps #/1000
theta_single_far <- obs_1_far$turns
################################################################################





#  Select only steps where whale was first sighted appraoching ship from the side
side_obs <- sep_input_dat %>%
                    filter(steps < 6001) %>%
                    filter(ship_whale_bearing < -60 | ship_whale_bearing > 60)
#  Select only steps where whale was first sighted appraoching ship from the front
front_obs <- sep_input_dat %>%
                     filter(steps < 6001) %>%
                     filter(ship_whale_bearing < 20 | ship_whale_bearing < -20)
#  Plot for an idea of what we're looking at...
side_obs_steps_tmp <- side_obs %>% 
                                        dplyr::select(steps) %>%
                                        filter(steps < 6001)
side_obs_steps <- side_obs_steps_tmp %>%
                                mutate(iter_num = 1:nrow(side_obs_steps_tmp)) %>%
                                mutate(type = "first sighting side steps")
front_obs_steps_tmp <- front_obs %>% 
                                     dplyr::select(steps) %>%
                                     filter(steps < 6001)
front_obs_steps <- front_obs_steps_tmp %>%
                                mutate(iter_num = 1:nrow(front_obs_steps_tmp)) %>%
                                mutate(type = "first sighting front steps")
both <- rbind(side_obs_steps, front_obs_steps)
steps_both_plot <- ggplot(both, aes(steps, colour = type, fill = type)) + 
                                geom_density(alpha= 0.7) +
                                xlab("Step length (km)") +
                                ylab("Frequency") +
                                xlim(c(0, 6000)) +
                                ylim(c(0, 0.006)) +
                                theme_bw() +
                                theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
steps_both_plot

side_obs_turns_tmp <- side_obs %>% 
                                      dplyr::select(turns, steps) %>%
                                      filter(steps < 6001)
side_obs_turns <- side_obs_turns_tmp %>%
                             dplyr::select(turns) %>%
                             mutate(iter_num = 1:nrow(side_obs_turns_tmp)) %>%
                             mutate(type = "first sighting side turns")
side_obs_turns$turns_deg <- ((side_obs_turns$turns)*180)/3.14159265359
front_obs_turns_tmp <- front_obs %>% 
                                       dplyr::select(turns, steps) %>%
                                       filter(steps < 6001)
front_obs_turns <- front_obs_turns_tmp %>%
                                dplyr::select(turns) %>%
                                mutate(iter_num = 1:nrow(front_obs_turns_tmp)) %>%
                                mutate(type = "first sighting front turns")
front_obs_turns$turns_deg <- ((front_obs_turns$turns)*180)/3.14159265359
both <- rbind(side_obs_turns, front_obs_turns)
turns_both_plot <- ggplot(both, aes(turns_deg, colour = type, fill = type)) + 
                                geom_density(alpha= 0.7) +
                                xlab("Turn angle (deg)") +
                                ylab("Frequency") +
                                #xlim(c(-3.15, 3.15)) +
                                #ylim(c(0, 0.003)) +
                                theme_bw() +
                                theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
turns_both_plot
################################################################################


#  Generate first sighting side only data and first sighting front only data for models
#  Organize data and give sequential IDs: 

#   First sighting side (60 deg or greater bearing from ship)
obs_1_side <- side_obs %>%
                        dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
                        arrange(ID_new) %>%
                        as.data.frame()

#   Indexing
npts_1_side <- nrow(obs_1_side)
ind_1_side <- obs_1_side$ID_new
nind_1_side <- length(unique(obs_1_side$ID_new))
nocc_1_side <- obs_1_side %>%
                          group_by(ID_new) %>%
                          summarise(nocc = n()) %>%
                          .$nocc

#  Data 
l_single_side <- obs_1_side$steps #/1000
theta_single_side <- obs_1_side$turns


#  First sighting front (20 deg or less bearing from ship)
obs_1_front <- front_obs %>%
                         dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
                         arrange(ID_new) %>%
                         as.data.frame()

#   Indexing
npts_1_front <- nrow(obs_1_front)
ind_1_front <- obs_1_front$ID_new
nind_1_front <- length(unique(obs_1_front$ID_new))
nocc_1_front <- obs_1_front %>%
                           group_by(ID_new) %>%
                           summarise(nocc = n()) %>%
                           .$nocc

#  Data 
l_single_front <- obs_1_front$steps #/1000
theta_single_front <- obs_1_front$turns
################################################################################


#  Using only data where first observation was a blow and subsequent obs were blow or dive.
tmp1_blow <- dat %>%
                             group_by(same_whale_ID) %>%
                             filter(first(whale_behavior)  == "BL-Blowing") %>%
                             ungroup() %>%
                             as.data.frame()
tmp2_blow <- tmp1_blow %>%
                            filter(whale_behavior  == "DF-Dive-fluke-up" | whale_behavior  == "DN-Dive-no-fluke" |
                            whale_behavior  == "BL-Blowing")
tmp2_blow$whale_behavior <- droplevels(tmp2_blow$whale_behavior)
locs_a_blow <- arrange(tmp2_blow, same_whale_ID, ob_order_time)
locs_b_blow <- locs_a_blow %>%
                            dplyr::select(same_whale_ID, X_whale_UTM, Y_whale_UTM) %>%
                            dplyr::rename(X = X_whale_UTM, Y = Y_whale_UTM)
locs_b_blow$same_whale_ID <- droplevels(locs_b_blow$same_whale_ID)

#   Create ltraj object
whale_traj_blow <- as.ltraj(xy = locs_b_blow[,c("X","Y")], id = locs_b_blow$same_whale_ID, typeII = FALSE)

#   Convert into dataframe
traj_dat_blow <- ld(whale_traj_blow)

#  Generate first sighting close only data and first sighting far only data for models
#  Organize data and give sequential IDs: 
input_dat_blow <- filter(traj_dat_blow, dist < 5000) 

#   First sighting is blow
obs_1_blow <- input_dat_blow %>%
                         dplyr::mutate(ID_new = as.numeric(as.factor(as.character(id)))) %>%
                         arrange(ID_new) %>%
                         as.data.frame()

#   Indexing
npts_1_blow <- nrow(obs_1_blow)
ind_1_blow <- obs_1_blow$ID_new
nind_1_blow <- length(unique(obs_1_blow$ID_new))
nocc_1_blow <- obs_1_blow %>%
                            group_by(ID_new) %>%
                            summarise(nocc = n()) %>%
                            .$nocc

#  Data 
l_single_blow <- obs_1_blow$dist #/1000
theta_single_blow <- obs_1_blow$rel.angle




























sm <- sep_input_dat %>%
          group_by(same_whale_ID) %>%
          filter(n() > 2) %>%
          as.data.frame()
sm_err <- sm %>%
                  mutate(min_step = steps - ()

                  
                  
sqrt((Y_whale_UTM - lag(Y_whale_UTM))²+(X_whale_UTM - lag(X_whale_UTM))²)