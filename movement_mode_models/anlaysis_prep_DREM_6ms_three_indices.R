# Sara Williams
# 9/18/2019; 2/2/2020; 10/23/2020
# Data prep 
################################################################################

#  Load packages
library(dplyr)
library(adehabitatLT)    
library(tidyr)
library(stringr)
library(aspace)
library(sf)

`%!in%` = function(x,y) !(x %in% y)

################################################################################

#  Load data
dat_raw <- read.csv("C:/Users/saraw/Documents/Whales/data/Whales_0615_general_clean.csv")
sa <- st_read("C:/Users/saraw/Documents/Whales/Tracks_0615/survey_area_buff15k.shp")

#  Data clean up and manipulation
tmp1 <- dat_raw %>%
	dplyr::filter(year > 2009) %>% ## per info from Karin, PTB (bearing) was kind of weird before 2010
	dplyr::filter(!is.na(X_whale_UTM)) %>%
	group_by(same_whale_ID) %>%
	filter(n() > 1) %>%
	as.data.frame()  %>%
	dplyr::filter(ship_whale_dist < 5001) %>%
	arrange(same_whale_ID, ob_order_time) 

ids_only1count <- tmp1 %>% 
	group_by(same_whale_ID) %>%
	mutate(max_count = max(count)) %>%
	as.data.frame() %>%
	dplyr::filter(max_count < 2)
ids_only1count <- as.character(unique(ids_only1count$same_whale_ID))

tmp2 <- tmp1 %>%
	dplyr::filter(same_whale_ID %in% ids_only1count) 

tmp3 <- tmp2 %>% 
	tidyr::separate(TimeTxt, c("hrs", "mins", "sec"), sep = ":", remove = FALSE) %>%
	mutate_if(is.character, as.numeric) %>%
	mutate(time_sec = sec + mins*60 + hrs*60*60) %>%
	group_by(same_whale_ID) %>%
	mutate(tdiff_s = abs(time_sec - lead(time_sec))) %>%
	as.data.frame() 
	
tmp4a <- tmp3 %>%
	dplyr::select(-tdiff_s) %>%
	mutate(X = X_whale_UTM, Y = Y_whale_UTM) %>%
	st_as_sf(coords = c("X", "Y"), crs = 32608) %>%
	group_by(same_whale_ID) %>%
	mutate(lead_tmp = geometry[row_number() + 1]) %>%
    mutate(dist_m = as.numeric(st_distance(geometry, lead_tmp, by_element = T))) %>%
    mutate(tdiff_s = abs(time_sec - lead(time_sec))) %>%
	mutate(spd_mps = dist_m/tdiff_s) %>%
	as.data.frame() %>% 
	dplyr::select(unique_event_ID, same_whale_ID, count,
		X_whale_UTM, Y_whale_UTM, 
		dist_m, tdiff_s, spd_mps, time_sec, ob_order_time, 
		ship_whale_dist,  ship_whale_bearing, 
		whale_behavior, visibility, TimeTxt) 

too_fast1 <- tmp4a %>% dplyr::filter(spd_mps > 6.7)		
rem_ids1 <- as.character(unique(too_fast1$unique_event_ID))
tmp4 <- tmp4a %>% dplyr::filter(unique_event_ID %!in% rem_ids1) 


tmp5a <- tmp4 %>%
	dplyr::select(unique_event_ID, same_whale_ID, count,
		X_whale_UTM, Y_whale_UTM, 
		time_sec, ob_order_time, 
		ship_whale_dist,  ship_whale_bearing, 
		whale_behavior, visibility, TimeTxt) %>%
	mutate(X = X_whale_UTM, Y = Y_whale_UTM) %>%
	st_as_sf(coords = c("X", "Y"), crs = 32608) %>%
	group_by(same_whale_ID) %>%
	filter(n() > 1) %>%
	mutate(lead_tmp = geometry[row_number() + 1]) %>%
    mutate(dist_m = as.numeric(st_distance(geometry, lead_tmp, by_element = T))) %>%
	mutate(tdiff_s = abs(time_sec - lead(time_sec))) %>%
	mutate(spd_mps = dist_m/tdiff_s) %>%
	as.data.frame() 

too_fast2 <- tmp5a %>% dplyr::filter(spd_mps > 6.7)
rem_ids2 <- as.character(unique(too_fast2$unique_event_ID))
tmp5 <- tmp5a %>% dplyr::filter(unique_event_ID %!in% rem_ids2) 


tmp6a <- tmp5 %>%
	dplyr::select(unique_event_ID, same_whale_ID, count,
		X_whale_UTM, Y_whale_UTM, 
		time_sec, ob_order_time, 
		ship_whale_dist,  ship_whale_bearing, 
		whale_behavior, visibility, TimeTxt) %>%
	mutate(X = X_whale_UTM, Y = Y_whale_UTM) %>%
	st_as_sf(coords = c("X", "Y"), crs = 32608) %>%
	group_by(same_whale_ID) %>%
	filter(n() > 1) %>%
	mutate(lead_tmp = geometry[row_number() + 1]) %>%
    mutate(dist_m = as.numeric(st_distance(geometry, lead_tmp, by_element = T))) %>%
	mutate(tdiff_s = abs(time_sec - lead(time_sec))) %>%
	mutate(spd_mps = dist_m/tdiff_s) %>%
	as.data.frame() 
	
too_fast3 <- tmp6a %>% dplyr::filter(spd_mps > 6.7)	
rem_ids3 <- as.character(unique(too_fast3$unique_event_ID))
tmp6 <- tmp6a %>% dplyr::filter(unique_event_ID %!in% rem_ids3) 


tmp7a <- tmp6 %>%
	dplyr::select(unique_event_ID, same_whale_ID, count,
		X_whale_UTM, Y_whale_UTM, 
		time_sec, ob_order_time, 
		ship_whale_dist,  ship_whale_bearing, 
		whale_behavior, visibility, TimeTxt) %>%
	mutate(X = X_whale_UTM, Y = Y_whale_UTM) %>%
	st_as_sf(coords = c("X", "Y"), crs = 32608) %>%
	group_by(same_whale_ID) %>%
	filter(n() > 1) %>%
	mutate(lead_tmp = geometry[row_number() + 1]) %>%
    mutate(dist_m = as.numeric(st_distance(geometry, lead_tmp, by_element = T))) %>%
	mutate(tdiff_s = abs(time_sec - lead(time_sec))) %>%
	mutate(spd_mps = dist_m/tdiff_s) %>%
	as.data.frame() 
	
too_fast4 <- tmp7a %>% dplyr::filter(spd_mps > 6.7)
rem_ids4 <- as.character(unique(too_fast4$unique_event_ID))
tmp7 <- tmp7a %>% dplyr::filter(unique_event_ID %!in% rem_ids4) 

### All movements that are greater than 6.7 m/s removed by this point

tmp8 <- tmp7 %>%
	dplyr::select(unique_event_ID, same_whale_ID, count,
		X_whale_UTM, Y_whale_UTM, 
		time_sec, ob_order_time, 
		ship_whale_dist,  ship_whale_bearing, 
		whale_behavior, visibility, TimeTxt) %>%
	mutate(X = X_whale_UTM, Y = Y_whale_UTM) %>%
	st_as_sf(coords = c("X", "Y"), crs = 32608) %>%
	group_by(same_whale_ID) %>%
	filter(n() > 1) %>%
	mutate(lead_tmp = geometry[row_number() + 1]) %>%
    mutate(dist_m = as.numeric(st_distance(geometry, lead_tmp, by_element = T))) %>%
	mutate(tdiff_s = abs(time_sec - lead(time_sec))) %>%
	mutate(spd_mps = dist_m/tdiff_s) %>%
	as.data.frame() 
tmp8_sf <- tmp8 %>%
	mutate(X = X_whale_UTM, Y = Y_whale_UTM) %>%
	st_as_sf(coords = c("X", "Y"), crs = 32608) 


#  Create trajectory (ltraj) object for steps and turns
tmp8$same_whale_ID <- droplevels(tmp8$same_whale_ID)
tmp9 <- tmp8 %>%
	dplyr::select(X_whale_UTM, Y_whale_UTM, 
		ob_order_time, TimeTxt, tdiff_s, 
		ship_whale_dist, ship_whale_bearing, whale_behavior, 
		same_whale_ID, unique_event_ID, time_sec) %>%
	arrange(same_whale_ID, ob_order_time)
traj <- as.ltraj(xy = tmp9[,c("X_whale_UTM","Y_whale_UTM")], 
	id = tmp9$same_whale_ID, typeII = FALSE)

#  Convert into dataframe
traj_df_tmp <- ld(traj)
traj_df <- traj_df_tmp %>%
	dplyr::select(id, x, y, dist, rel.angle, abs.angle, dx, dy) 
names(traj_df) <- c("same_whale_ID", "X_whale_UTM", "Y_whale_UTM", 
	"step", "turn", "abs_angle", "dx", "dy")

#  Join trajectory dataframe with dataframe holding location variables
dat <- full_join(traj_df, tmp9, by = c("same_whale_ID", "X_whale_UTM", "Y_whale_UTM")) %>%
	mutate(spd_mps = step/tdiff_s) %>%
	group_by(same_whale_ID) %>%
	dplyr::mutate(occ = 1:n()) %>%
	dplyr::mutate(n_occ = n()) %>%
	as.data.frame()


#  Data summary
# of surfacings
nrow(dat)
# of unique surfacing bouts/whales
length(unique(dat$same_whale_ID))

# Constrain to 99% quantiles for surfacing lengths 
sl_99 <- as.numeric(quantile(dat$step[!is.na(dat$step)], probs = c(0.99)))
# sa_99 <- as.numeric(quantile(dat$turn[!is.na(dat$turn)], probs = c(0.99)))

kp_sl <- which(dat$step <= sl_99) 
dat99 <- dat[kp_sl,]


#  Format data for use in model with parameters estaimted per index
mod_dat <- dat99 %>%
	dplyr::filter(!is.na(step)) %>%
	mutate(dist_ind = ifelse(ship_whale_dist < 1001, 1, 
		ifelse(ship_whale_dist >= 1001 & ship_whale_dist < 3001, 2, 3))) %>%
	mutate(bear_ind = ifelse(abs(ship_whale_bearing) < 22.56, 1, 
		ifelse(abs(ship_whale_bearing) >= 22.56 & abs(ship_whale_bearing) < 45, 2, 3))) %>%
	mutate(beh_ind = ifelse(whale_behavior == "BL-Blowing" | whale_behavior == "DF-Dive0fluke-up" |
		whale_behavior == "DN-Dive-no-fluke", 1, 2))  #%>%
	#mutate(ind = ifelse(dist_ind == 1 & bear_ind == 1, 1, # close and front
		#ifelse(dist_ind == 1 & bear_ind == 3, 2, # close and side
		#ifelse(dist_ind == 3 & bear_ind == 1, 3, # far and front
		#ifelse(dist_ind == 3 & bear_ind == 3, 4,  NA))))) # far and side
	


#  Format data for use in model with parameters estaimted per index
npts <- nrow(mod_dat)
ind <- mod_dat$ID_new
nind <- length(unique(mod_dat$ID_new))

l <- (mod_dat$step)/1000
theta <- mod_dat$turn
dist_ind <- mod_dat$dist_ind
bear_ind <- mod_dat$bear_ind
#ind <- mod_dat$ind

