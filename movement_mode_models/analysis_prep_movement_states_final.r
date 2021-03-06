# Sara Williams
# 12/8/2015; updated 2/1/2016, 3/9/2016, 12/16/2016; 1/17/2017; 9/4/2019
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
dat_raw <- read.csv("C:/Users/saraw/Documents/Whales/data/Whales_0615_general_clean.csv")

#  Data clean up and manipulation
tmp1 <- dat_raw %>%
	filter(year > 2009) %>% ## per info from Karin, PTB (bearing) was kind of weird before 2010
	#filter(ship_whale_dist < 8000) %>%
	filter(count == 1) %>%            
	group_by(same_whale_ID) %>%
	filter(n() > 1) %>%
	ungroup() %>%
	#dplyr::filter(same_whale_ID != "2015-08-21-K-034") %>%
	as.data.frame()
tmp2 <- arrange(tmp1, same_whale_ID, ob_order_time)
tmp3 <- tmp2 %>%
	dplyr::select(same_whale_ID, X_whale_UTM, Y_whale_UTM, ship_whale_dist, 
		ship_whale_bearing, whale_behavior, ob_order_time, TimeTxt, same_whale_ID) %>%
	dplyr::rename(X = X_whale_UTM, Y = Y_whale_UTM)

#  Add time difference between successive observations to observations.
t1 <- tmp3
t2 <- as.data.frame(str_split_fixed(t1$TimeTxt, ":", 3))
t3 <- cbind(t2, t1$X, t1$Y, t1$ob_order_time, t1$ship_whale_dist, 
                    t1$ship_whale_bearing, t1$whale_behavior, t1$TimeTxt, t1$same_whale_ID)
names(t3)[4] <- "X"
names(t3)[5] <- "Y"
names(t3)[6] <- "ob_order_time"
names(t3)[7] <- "ship_whale_dist"
names(t3)[8] <- "ship_whale_bearing"
names(t3)[9] <- "whale_behavior"
names(t3)[10] <- "time_txt"
names(t3)[11] <- "same_whale_ID"
t3$V1 <- as.numeric(as.integer(as.character(t3$V1)))
t3$V2 <- as.numeric(as.integer(as.character(t3$V2)))
t3$V3 <- as.numeric(as.integer(as.character(t3$V3)))
#   time_diff is the difference in time between a sighting and the sighting before it
t4 <- t3 %>%
	dplyr::rename (hr = V1, min = V2, sec = V3) %>%
	dplyr::mutate(sec_frac = sec/60, min_sec = sec_frac+min) %>%
	dplyr::mutate(min_frac = min_sec/60, time_dec = hr+min_frac) %>%
	group_by(same_whale_ID) %>%
	mutate(time_diff = lead(time_dec) - time_dec) %>%
	mutate(time_diff_sec = time_diff*3600) %>%
	ungroup() %>%
	as.data.frame()
		 
tot_time_obs <- t4 %>%
	group_by(same_whale_ID) %>%
	mutate(n_obs = max(ob_order_time)) %>%
	mutate(surface_event_time = sum(time_diff_sec, na.rm = TRUE)) %>%
	slice(1) %>%
	as.data.frame(.)
tot_time_obs95 <- quantile(tot_time_obs$surface_event_time, probs = c(0.95))
tot_time_obs90 <- quantile(tot_time_obs$surface_event_time, probs = c(0.90))
tot_time_obs75 <- quantile(tot_time_obs$surface_event_time, probs = c(0.75))
	
################################################################################

#  Create trajectory (ltraj) object for steps and turns
tmp4 <- t4 %>%
	dplyr::select(X, Y, ob_order_time, time_txt, time_dec, ship_whale_dist,
		ship_whale_bearing, whale_behavior, same_whale_ID, time_diff_sec)
tmp4$same_whale_ID <- droplevels(tmp4$same_whale_ID)
traj <- as.ltraj(xy = tmp4[,c("X","Y")], id = tmp4$same_whale_ID, typeII = FALSE)

	#  With dates dummy coded for seconds between obs
	#tmp5 <- tmp4 %>%
	#	group_by(same_whale_ID) %>%
	#	mutate(row_num = row_number()) %>%
	#	mutate(dummy_date_idx =ifelse(row_num == 1, 1, 
	#		ifelse(row_num == 2, 1 + lag(time_diff_sec), 
	#		ifelse(row_num == 3, 1 + lag(time_diff_sec) + lag(time_diff_sec, 2), 
	#		ifelse(row_num == 4, 1 + lag(time_diff_sec) + lag(time_diff_sec, 2) + lag(time_diff_sec, 3), 
	#		ifelse(row_num == 5, 1 + lag(time_diff_sec) + lag(time_diff_sec, 2) + lag(time_diff_sec, 3) + lag(time_diff_sec, 4), 
	#		ifelse(row_num == 6, 1 + lag(time_diff_sec) + lag(time_diff_sec, 2) + lag(time_diff_sec, 3) + 
	#			lag(time_diff_sec, 4) + lag(time_diff_sec, 5), 
	#		ifelse(row_num == 7, 1 + lag(time_diff_sec) + lag(time_diff_sec, 2) + lag(time_diff_sec, 3) + 
	#			lag(time_diff_sec, 4) + lag(time_diff_sec, 5) + lag(time_diff_sec, 6),
	#		ifelse(row_num == 8, 1 + lag(time_diff_sec) + lag(time_diff_sec, 2) + lag(time_diff_sec, 3) + 
	#			lag(time_diff_sec, 4) + lag(time_diff_sec, 5) + lag(time_diff_sec, 6)  + lag(time_diff_sec, 7),
	#		ifelse(row_num == 9, 1 + lag(time_diff_sec) + lag(time_diff_sec, 2) + lag(time_diff_sec, 3) + 
	#			lag(time_diff_sec, 4) + lag(time_diff_sec, 5) + lag(time_diff_sec, 6)  + lag(time_diff_sec, 7) + lag(time_diff_sec, 8),
	#		ifelse(row_num == 10, 1 + lag(time_diff_sec) + lag(time_diff_sec, 2) + lag(time_diff_sec, 3) + 
	#			lag(time_diff_sec, 4) + lag(time_diff_sec, 5) + lag(time_diff_sec, 6)  + lag(time_diff_sec, 7) + lag(time_diff_sec, 8) +
	#			lag(time_diff_sec, 9),
	#		ifelse(row_num == 11, 1 + lag(time_diff_sec) + lag(time_diff_sec, 2) + lag(time_diff_sec, 3) + 
	#			lag(time_diff_sec, 4) + lag(time_diff_sec, 5) + lag(time_diff_sec, 6)  + lag(time_diff_sec, 7) + lag(time_diff_sec, 8) +
	#			lag(time_diff_sec, 9) + lag(time_diff_sec, 10), 1 + lag(time_diff_sec) + lag(time_diff_sec, 2) + lag(time_diff_sec, 3) + 
	#			lag(time_diff_sec, 4) + lag(time_diff_sec, 5) + lag(time_diff_sec, 6)  + lag(time_diff_sec, 7) + lag(time_diff_sec, 8) +
	#			lag(time_diff_sec, 9) + lag(time_diff_sec, 10) + lag(time_diff_sec, 11))))))))))))) %>%
	#	as.data.frame()
	#	
	#needed_dates <- seq(from = as.POSIXct("2010/1/1"), 
	#	by = "day", length.out = max(tmp5$dummy_date_idx, na.rm = TRUE))
	#dummy_dates <- as.data.frame(seq(1:max(tmp5$dummy_date_idx, na.rm = TRUE))) %>%
	#	bind_cols(as.data.frame(needed_dates))
	#colnames(dummy_dates) <- c("dummy_date_idx", "date")	
	#
	#tmp5$dummy_date_idx <- round(tmp5$dummy_date_idx)
	#tmp5$dummy_date_idx <- as.character(round(tmp5$dummy_date_idx))
	#dummy_dates$dummy_date_idx <- as.character(dummy_dates$dummy_date_idx)
	#
	#tmp6 <- tmp5 %>%
	#	left_join(dummy_dates, by = "dummy_date_idx")
	#
	#tmp7 <- tmp6 %>%
	#  group_by(same_whale_ID) %>%
	#  mutate(count = n_distinct(date)) %>%
	#  mutate(n_obs = n()) %>%
	#  as.data.frame() %>%
	#  mutate(date_obs_diff = ifelse(count == n_obs, 0, 1)) %>%
	#  filter(date_obs_diff == 1)
	#  
	#traj <- as.ltraj(xy = tmp6[,c("X","Y")], id = tmp6$same_whale_ID, burst = tmp6$same_whale_ID,
		#typeII = TRUE, date = tmp6$date)






#  Convert into dataframe
traj_df_tmp <- ld(traj)
traj_df <- traj_df_tmp %>%
                 dplyr::select(id, x, y, dist, rel.angle, abs.angle, dx, dy) 
names(traj_df) <- c("same_whale_ID", "X", "Y", "step", "turn", "abs_angle", "dx", "dy")

#  Join trajectory dataframe with dataframe holding location variables
mod_dat_tmp <- full_join(traj_df, tmp4, by = c("same_whale_ID", "X", "Y")) %>%
	mutate(ave_swim_spd = step/time_diff_sec)
#  Remove erroneous steps based on 90th or 95th quantile
quantile(mod_dat_tmp$step, probs = c(0.90, 0.95), na.rm = TRUE)
#quantile(mod_dat_tmp$ave_swim_spd, probs = c(0.90, 0.95), na.rm = TRUE)
#mod_dat <- filter(mod_dat_tmp, step < 2086 | is.na(step))
mod_dat <- mod_dat_tmp %>%
	dplyr::filter(ave_swim_spd < 6.6877) %>% #4.361111
	dplyr::filter(!is.na(step))


#mod_dat <- filter(mod_dat, ave_swim_spd < 4.361111)


################################################################################

#  Subset data based on wanting to assess influence of ship and whale activity

#  Whale activity
#   Deep dive > 120s, surface interval dive < 50
surf <- filter(mod_dat, time_diff_sec < 50)
#surf <- filter(mod_dat, time_diff_sec < 30)
deep <- filter(mod_dat, time_diff_sec > 120)

#  Whale activity 2
#   Fluke dive or surface activity/feeding behavior described by visual assessment
#fluke <- mod_dat %>%
#	group_by(same_whale_ID) %>%
#	filter(any(whale_behavior == "DF-Dive-fluke-up")) %>%
#	filter(!any(whale_behavior == "SA-Surface-active")) %>%
#	as.data.frame()
#surf_act <- mod_dat %>%
#	group_by(same_whale_ID) %>%
#	filter(any(whale_behavior == "SA-Surface-active")) %>%
#	filter(!any(whale_behavior == "DF-Dive-fluke-up")) %>%
#	as.data.frame()
	

#  Distance to ship at first sighting
#   Near (<1000m), far away (> 3000m)
near_ids <- mod_dat %>%
	group_by(same_whale_ID) %>%
	slice(1) %>%
	filter(ship_whale_dist < 1000) %>%
	ungroup() %>%
	as.data.frame(.)
near <- semi_join(mod_dat, near_ids, by = c("same_whale_ID"))
#near <- semi_join(mod_dat, near_ids, by = c("same_whale_ID", "X", "Y")  ?? correct join?
far_ids <- mod_dat %>%
	group_by(same_whale_ID) %>%
	slice(1) %>%
	filter(ship_whale_dist > 3000) %>%
	ungroup() %>%
	as.data.frame(.)
far <- semi_join(mod_dat, far_ids, by = "same_whale_ID") 

#  Relative bearing (approach) at first sighting
#   Side approach( greater than +/- 40 deg relative to bow), front approach (less than +/- 20)
side_ids <- mod_dat %>%
	group_by(same_whale_ID) %>%
	slice(1) %>%
	filter(ship_whale_bearing < -40 | ship_whale_bearing > 40) %>%
	ungroup() %>%
	as.data.frame(.)
side <- semi_join(mod_dat, side_ids, by = c("same_whale_ID"))
front_ids <- mod_dat %>%
	group_by(same_whale_ID) %>%
	slice(1) %>%
	filter(ship_whale_bearing < 20 | ship_whale_bearing < -20) %>%
	ungroup() %>%
	as.data.frame(.)
front <- semi_join(mod_dat, front_ids, by = c("same_whale_ID"))


#  Relative bearing (approach) at first sighting AND first sighting within 1k
#   Side approach( greater than +/- 40 deg relative to bow), front approach (less than +/- 20)
side_1k_ids <- near %>%
                    group_by(same_whale_ID) %>%
                    slice(1) %>%
                    filter(ship_whale_bearing < -40 | ship_whale_bearing > 40) %>%
					#filter(ship_whale_dist < 1001) %>%
                    ungroup() %>%
                    as.data.frame(.)
side_1k <- semi_join(mod_dat, side_1k_ids, by = c("same_whale_ID"))
front_1k_ids <- near %>%
                    group_by(same_whale_ID) %>%
                    slice(1) %>%
                    filter(ship_whale_bearing < 20 | ship_whale_bearing < -20) %>%
					#filter(ship_whale_dist < 1001) %>%
                    ungroup() %>%
                    as.data.frame(.)
front_1k <- semi_join(mod_dat, front_1k_ids, by = c("same_whale_ID"))
################################################################################

#  Generate data for use in model runs

#   All data single
obs_1 <- mod_dat %>%
	group_by(same_whale_ID) %>%
	dplyr::mutate(occ = 1:n()) %>%
	ungroup() %>%
	dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
	arrange(ID_new) %>%
	as.data.frame()
npts_1 <- nrow(obs_1)
ind_1 <- obs_1$ID_new
nind_1 <- length(unique(obs_1$ID_new))
nocc_1 <- obs_1 %>%
	group_by(ID_new) %>%
	summarise(nocc = n()) %>%
	.$nocc
l <- (obs_1$step)/1000
theta <- obs_1$turn

# DREM
obs_1 <- mod_dat %>%
	group_by(same_whale_ID) %>%
	dplyr::mutate(occ = 1:n()) %>%
	ungroup() %>%
	dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
	arrange(ID_new) %>%
	as.data.frame() %>%
	#mutate(mean_ship_whale_dist = mean(ship_whale_dist, na.rm = TRUE)) %>%
	#mutate(mean_abs_rel_bearing = mean(abs(ship_whale_bearing), na.rm = TRUE)) %>%
	#mutate(DREM_dist_ind = ifelse(ship_whale_dist < mean_ship_whale_dist, 1, 2)) %>%
	#mutate(DREM_bear_ind = ifelse(ship_whale_bearing < mean_abs_rel_bearing, 1, 2)) %>%
	mutate(DREM_dist_ind = ifelse(ship_whale_dist < 1001, 1, 
		ifelse(ship_whale_dist >= 1001 & ship_whale_dist < 3001, 2, 3))) %>%
	mutate(DREM_bear_ind = ifelse(abs(ship_whale_bearing) < 22.56, 1,
		ifelse(abs(ship_whale_bearing) >= 22.56 & abs(ship_whale_bearing) < 45, 2, 3))) %>%
	mutate(DREM_dur_ind = ifelse(time_diff_sec < 50, 1,
		ifelse(time_diff_sec >= 50 & time_diff_sec < 120, 2, 3))) %>%
	dplyr::filter(!is.na(X))
	
npts_1 <- nrow(obs_1)
ind_1 <- obs_1$ID_new
nind_1 <- length(unique(obs_1$ID_new))
nocc_1 <- obs_1 %>%
	group_by(ID_new) %>%
	summarise(nocc = n()) %>%
	.$nocc
l <- (obs_1$step)/1000
theta <- obs_1$turn
dist_ind <- obs_1$DREM_dist_ind
bear_ind <- obs_1$DREM_bear_ind
dur_ind <- obs_1$DREM_dur_ind

#  All data double (made square)
l_double <- obs_1 %>%
	mutate(step_km = (step/1000)) %>%
	dplyr::select(ID_new, occ, step_km) %>%
	spread(occ, step_km, fill = NA, convert = FALSE)
theta_double <- obs_1 %>%
	dplyr::select(ID_new, occ, turn) %>%
	spread(occ, turn, fill = NA, convert = FALSE)
ship_dist_double <- obs_1 %>%
	mutate(ship_dist_scale = scale(ship_whale_dist)) %>%
	dplyr::select(ID_new, occ, ship_dist_scale) %>%
	spread(occ, ship_dist_scale, fill = NA, convert = FALSE)
                                
#   Surface interval
obs_1_surf <- surf %>%
	dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
	arrange(ID_new) %>%
	as.data.frame()
npts_1_surf <- nrow(obs_1_surf)
ind_1_surf <- obs_1_surf$ID_new
nind_1_surf <- length(unique(obs_1_surf$ID_new))
nocc_1_surf <- obs_1_surf %>%
	group_by(ID_new) %>%
	summarise(nocc = n()) %>%
	.$nocc
l_surf <- (obs_1_surf$step)/1000
theta_surf <- obs_1_surf$turn

#   Deep dive
obs_1_dive <- deep %>%
	dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
	arrange(ID_new) %>%
	as.data.frame()
npts_1_dive <- nrow(obs_1_dive)
ind_1_dive <- obs_1_dive$ID_new
nind_1_dive <- length(unique(obs_1_dive$ID_new))
nocc_1_dive <- obs_1_dive %>%
	group_by(ID_new) %>%
	summarise(nocc = n()) %>%
	.$nocc
l_dive <- (obs_1_dive$step)/1000
theta_dive <- obs_1_dive$turn


##   Surface active or feeding behavior
#obs_1_surf_act <- surf_act %>%
#	dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
#	arrange(ID_new) %>%
#	as.data.frame()
#npts_1_surf_act <- nrow(obs_1_surf_act)
#ind_1_surf_act <- obs_1_surf_act$ID_new
#nind_1_surf_act <- length(unique(obs_1_surf_act$ID_new))
#nocc_1_surf_act <- obs_1_surf_act %>%
#	group_by(ID_new) %>%
#	summarise(nocc = n()) %>%
#	.$nocc
#l_surf_act <- (obs_1_surf_act$step)/1000
#theta_surf_act <- obs_1_surf_act$turn
#
##   Fluke dive
#obs_1_fluke <- fluke %>%
#	dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
#	arrange(ID_new) %>%
#	as.data.frame()
#npts_1_fluke <- nrow(obs_1_fluke)
#ind_1_fluke <- obs_1_fluke$ID_new
#nind_1_fluke <- length(unique(obs_1_fluke$ID_new))
#nocc_1_fluke <- obs_1_fluke %>%
#	group_by(ID_new) %>%
#	summarise(nocc = n()) %>%
#	.$nocc
#l_fluke <- (obs_1_fluke$step)/1000
#theta_fluke <- obs_1_fluke$turn
#
#

#   First sighting near (<1000m)
obs_1_near <- near %>%
	dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
	arrange(ID_new) %>%
	as.data.frame()
npts_1_near <- nrow(obs_1_near)
ind_1_near <- obs_1_near$ID_new
nind_1_near <- length(unique(obs_1_near$ID_new))
nocc_1_near <- obs_1_near %>%
	group_by(ID_new) %>%
	summarise(nocc = n()) %>%
	.$nocc
l_near <- (obs_1_near$step)/1000
theta_near <- obs_1_near$turn

#   First sighting far ( > 3000m)
obs_1_far <- far %>%
	dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
	arrange(ID_new) %>%
	as.data.frame()
npts_1_far <- nrow(obs_1_far)
ind_1_far <- obs_1_far$ID_new
nind_1_far <- length(unique(obs_1_far$ID_new))
nocc_1_far <- obs_1_far %>%
	group_by(ID_new) %>%
	summarise(nocc = n()) %>%
	.$nocc
l_far <- (obs_1_far$step)/1000
theta_far <- obs_1_far$turn

#   First sighting approaching from side 
obs_1_side <- side %>%
	dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
	arrange(ID_new) %>%
	as.data.frame()
npts_1_side <- nrow(obs_1_side)
ind_1_side <- obs_1_side$ID_new
nind_1_side <- length(unique(obs_1_side$ID_new))
nocc_1_side <- obs_1_side %>%
	group_by(ID_new) %>%
	summarise(nocc = n()) %>%
	.$nocc
l_side <- (obs_1_side$step)/1000
theta_side <- obs_1_side$turn

#   First sighting approaching from front
obs_1_front <- front %>%
	dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
	arrange(ID_new) %>%
	as.data.frame()
npts_1_front <- nrow(obs_1_front)
ind_1_front <- obs_1_front$ID_new
nind_1_front <- length(unique(obs_1_front$ID_new))
nocc_1_front <- obs_1_front %>%
	group_by(ID_new) %>%
	summarise(nocc = n()) %>%
	.$nocc
l_front <- (obs_1_front$step)/1000
theta_front <- obs_1_front$turn


#
##   First sighting approaching from side AND first sighting within 1k
#obs_1_side_1k <- side_1k %>%
#	dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
#	arrange(ID_new) %>%
#	as.data.frame()
#npts_1_side_1k <- nrow(obs_1_side_1k)
#ind_1_side_1k <- obs_1_side_1k$ID_new
#nind_1_side_1k <- length(unique(obs_1_side_1k$ID_new))
#nocc_1_side_1k <- obs_1_side_1k %>%
#	group_by(ID_new) %>%
#	summarise(nocc = n()) %>%
#	.$nocc
#l_side_1k <- (obs_1_side_1k$step)/1000
#theta_side_1k <- obs_1_side_1k$turn

##   First sighting approaching from front AND first sighting within 1k
#obs_1_front_1k <- front_1k %>%
#	dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
#	arrange(ID_new) %>%
#	as.data.frame()
#npts_1_front_1k <- nrow(obs_1_front_1k)
#ind_1_front_1k <- obs_1_front_1k$ID_new
#nind_1_front_1k <- length(unique(obs_1_front_1k$ID_new))
#nocc_1_front_1k <- obs_1_front_1k %>%
#	group_by(ID_new) %>%
#	summarise(nocc = n()) %>%
#	.$nocc
#l_front_1k <- (obs_1_front_1k$step)/1000
#theta_front_1k <- obs_1_front_1k$turn
################################################################################





















####### OTHER IDEAS AND SUBSETS....Not currently used

#  Histogram for number of cues per individual
tmp5 <- tmp4 %>%
              group_by(same_whale_ID) %>%
              mutate(n_cue = n()) %>%
              ungroup() %>%
              as.data.frame()



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

#  Generate first sighting near only data and first sighting far only data for models
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
























################################################################################
#  Plot for an idea of what we're looking at...
surf_steps_tmp <- short_int %>% 
                       dplyr::select(step) %>%
                       filter(step < 6001)
surf_steps <- surf_steps_tmp %>%
                        mutate(iter_num = 1:nrow(surf_steps_tmp)) %>%
                        mutate(type = "surfacing interval step")
dive_steps_tmp <- long_int_full_path_long %>% 
                                dplyr::select(step) %>%
                                filter(step < 6001)
dive_steps <- dive_steps_tmp %>%
                        mutate(iter_num = 1:nrow(dive_steps_tmp)) %>%
                        mutate(type = "deep dive step")
both <- rbind(surf_steps, dive_steps)
steps_both_plot <- ggplot(both, aes(step, colour = type, fill = type)) + 
                                  geom_density(alpha= 0.7) +
                                  xlab("Step length (km)") +
                                  ylab("Frequency") +
                                  xlim(c(0, 6000)) +
                                  ylim(c(0, 0.005)) +
                                  theme_bw() +
                                  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                  axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
steps_both_plot

surf_turns_tmp <- short_int %>% 
                                dplyr::select(turn, step) %>%
                                filter(step < 6001)
surf_turns <- surf_turns_tmp %>%
                        dplyr::select(turn) %>%
                        mutate(iter_num = 1:nrow(surf_turns_tmp)) %>%
                        mutate(type = "surfacing interval turn")
surf_turns$turns_deg <- ((surf_turns$turn)*180)/3.14159265359
dive_turns_tmp <- long_int_full_path_long %>% 
                                 dplyr::select(turn, step) %>%
                                 filter(step < 6001)
dive_turns <- dive_turns_tmp %>%
                        dplyr::select(turn) %>%
                        mutate(iter_num = 1:nrow(dive_turns_tmp)) %>%
                        mutate(type = "deep dive turn")
dive_turns$turns_deg <- ((dive_turns$turn)*180)/3.14159265359
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


#  Plot for an idea of what we're looking at...
side_obs_steps_tmp <- side_obs %>% 
                                        dplyr::select(step) %>%
                                        filter(step < 6001)
side_obs_steps <- side_obs_steps_tmp %>%
                                mutate(iter_num = 1:nrow(side_obs_steps_tmp)) %>%
                                mutate(type = "first sighting side step")
front_obs_steps_tmp <- front_obs %>% 
                                     dplyr::select(step) %>%
                                     filter(step < 6001)
front_obs_steps <- front_obs_steps_tmp %>%
                                mutate(iter_num = 1:nrow(front_obs_steps_tmp)) %>%
                                mutate(type = "first sighting front step")
both <- rbind(side_obs_steps, front_obs_steps)
steps_both_plot <- ggplot(both, aes(step, colour = type, fill = type)) + 
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
                                      dplyr::select(turn, step) %>%
                                      filter(step < 6001)
side_obs_turns <- side_obs_turns_tmp %>%
                             dplyr::select(turn) %>%
                             mutate(iter_num = 1:nrow(side_obs_turns_tmp)) %>%
                             mutate(type = "first sighting side turn")
side_obs_turns$turns_deg <- ((side_obs_turns$turn)*180)/3.14159265359
front_obs_turns_tmp <- front_obs %>% 
                                       dplyr::select(turn, step) %>%
                                       filter(step < 6001)
front_obs_turns <- front_obs_turns_tmp %>%
                                dplyr::select(turn) %>%
                                mutate(iter_num = 1:nrow(front_obs_turns_tmp)) %>%
                                mutate(type = "first sighting front turn")
front_obs_turns$turns_deg <- ((front_obs_turns$turn)*180)/3.14159265359
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



#  Plot for an idea of what we're looking at...
near_obs_steps_tmp <- near_obs_full_path %>% 
                                          dplyr::select(step) %>%
                                          filter(step < 6001)
near_obs_steps <- near_obs_steps_tmp %>%
                                  mutate(iter_num = 1:nrow(near_obs_steps_tmp)) %>%
                                  mutate(type = "first sighting near step")
far_obs_steps_tmp <- far_obs_full_path %>% 
                                      dplyr::select(step) %>%
                                      filter(step < 6001)
far_obs_steps <- far_obs_steps_tmp %>%
                              mutate(iter_num = 1:nrow(far_obs_steps_tmp)) %>%
                              mutate(type = "first sighting far step")
both <- rbind(near_obs_steps, far_obs_steps)
steps_both_plot <- ggplot(both, aes(step, colour = type, fill = type)) + 
                                  geom_density(alpha= 0.7) +
                                  xlab("Step length (km)") +
                                  ylab("Frequency") +
                                  xlim(c(0, 6000)) +
                                  ylim(c(0, 0.006)) +
                                  theme_bw() +
                                  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                  axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
steps_both_plot

near_obs_turns_tmp <- near_obs_full_path %>% 
                                        dplyr::select(turn, step) %>%
                                        filter(step < 6001)
near_obs_turns <- near_obs_turns_tmp %>%
                               dplyr::select(turn) %>%
                               mutate(iter_num = 1:nrow(near_obs_turns_tmp)) %>%
                               mutate(type = "first sighting near turn")
near_obs_turns$turns_deg <- ((near_obs_turns$turn)*180)/3.14159265359
far_obs_turns_tmp <- far_obs_full_path %>% 
                                     dplyr::select(turn, step) %>%
                                     filter(step < 6001)
far_obs_turns <- far_obs_turns_tmp %>%
                            dplyr::select(turn) %>%
                            mutate(iter_num = 1:nrow(far_obs_turns_tmp)) %>%
                            mutate(type = "first sighting far turn")
far_obs_turns$turns_deg <- ((far_obs_turns$turn)*180)/3.14159265359
both <- rbind(near_obs_turns, far_obs_turns)
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















# # obs_path <- ggplot() + 
                     # # geom_path(data = locs_b[1:1000,], aes(x = X, y = Y, group = same_whale_ID, colour = same_whale_ID),
                     # # show.legend = FALSE) +
                     # # #geom_path(data = df_XY2, aes(x = X, y = Y, group = walk_num)) +
                     # # xlab("X (km)") +
                     # # xlim(c(432500, 450000)) +
                     # # ylab("Y (km)") +
                     # # ylim(c(6460000, 6500000)) +
                     # # theme_bw()
# # obs_path