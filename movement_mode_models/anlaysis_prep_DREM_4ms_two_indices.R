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
	#filter(approx == "n") %>%
	as.data.frame()  %>%
	arrange(same_whale_ID, ob_order_time) 

ids_only1count <- tmp1 %>% 
	group_by(same_whale_ID) %>%
	mutate(max_count = max(count)) %>%
	as.data.frame() %>%
	dplyr::filter(max_count < 2)
	#dplyr::filter(max_count < 4)
ids_only1count <- as.character(unique(ids_only1count$same_whale_ID))

tmp2a <- tmp1 %>%
	dplyr::filter(same_whale_ID %in% ids_only1count) 

#  Whale IDs for those that were within 3000m at first sighting
kp_close_IDs <- tmp2a %>%
	group_by(same_whale_ID) %>%
	arrange(ob_order_time) %>%
	slice(1) %>%
	as.data.frame() %>%
	#dplyr::filter(ship_whale_dist < 3001)
	dplyr::filter(ship_whale_dist < 2001)
kp_close_IDs <- as.character(kp_close_IDs$same_whale_ID)

tmp2b <- tmp2a %>%
	dplyr::filter(same_whale_ID %in% kp_close_IDs)

kp_cone_IDs <- tmp2b %>%
	#mutate(abs_bear = abs(ship_whale_bearing)) %>%
	group_by(same_whale_ID) %>%
	#arrange(ob_order_time) %>%
	#mutate(min_abs_bear = min(abs_bear)) %>%
	#mutate(ever_in_cone = ifelse(min_abs_bear < 46, 1, 0)) %>%
	slice(1) %>%
	as.data.frame() %>%
	dplyr::filter(ship_whale_bearing < 45 & ship_whale_bearing > -45) 
	#dplyr::filter(ever_in_cone == 1) 
kp_cone_IDs <- as.character(kp_cone_IDs$same_whale_ID)

tmp2c <- tmp2b %>%
	dplyr::filter(same_whale_ID %in% kp_cone_IDs) %>%
	mutate(ob_time = mdy_hms(paste(DateD, TimeTxt)))

tmp3 <- tmp2c %>% 
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

too_fast1 <- tmp4a %>% dplyr::filter(spd_mps > 4.8) #7.5)		
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

too_fast2 <- tmp5a %>% dplyr::filter(spd_mps > 4.8) #7.5)		
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
	
too_fast3 <- tmp6a %>% dplyr::filter(spd_mps > 4.8) #7.5)		
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
	
too_fast4 <- tmp7a %>% dplyr::filter(spd_mps > 4.8) #7.5)		
rem_ids4 <- as.character(unique(too_fast4$unique_event_ID))
tmp7 <- tmp7a %>% dplyr::filter(unique_event_ID %!in% rem_ids4) 

### All movements that are greater than 4.8 m/s removed by this point

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
traj <- as.ltraj(xy = tmp8[,c("X_whale_UTM","Y_whale_UTM")], 
	id = tmp8$same_whale_ID, typeII = FALSE)

#  Convert into dataframe
traj_df_tmp <- ld(traj)
traj_df <- traj_df_tmp %>%
	dplyr::select(id, x, y, dist, rel.angle, abs.angle, dx, dy) 
names(traj_df) <- c("same_whale_ID", "X_whale_UTM", "Y_whale_UTM", 
	"step", "turn", "abs_angle", "dx", "dy")

#  Join trajectory dataframe with dataframe holding location variables
dat <- full_join(traj_df, tmp8, by = c("same_whale_ID", "X_whale_UTM", "Y_whale_UTM")) %>%
	mutate(spd_mps = step/tdiff_s) %>%
	group_by(same_whale_ID) %>%
	dplyr::mutate(occ = 1:n()) %>%
	dplyr::mutate(n_occ = n()) %>%
	as.data.frame()


# Constrain to 99% quantiles for surfacing lengths 
sl_99 <- as.numeric(quantile(dat$step[!is.na(dat$step)], probs = c(0.99)))
# sa_99 <- as.numeric(quantile(dat$turn[!is.na(dat$turn)], probs = c(0.99)))
kp_sl  <- which(dat$step <= sl_99) 
dat99 <- dat[kp_sl,]

#  Data summary
# of surfacings
nrow(dat)
nrow(dat99)
# of unique surfacing bouts/whales
length(unique(dat$same_whale_ID))
length(unique(dat99$same_whale_ID))
# of steps
length(dat$step[!is.na(dat$step)])
length(dat99$step[!is.na(dat99$step)])
# of turns
length(dat$turn[!is.na(dat$turn)])
length(dat99$turn[!is.na(dat99$turn)])



# Number of surfacings
nrow(dat)

# Number of unique surfacing bouts/whales
length(unique(dat$same_whale_ID))

# Number of surfacing lengths
length(dat$step[!is.na(dat$step)])

# Number of surfacing angles
length(dat$turn[!is.na(dat$turn)])



rem_sl <- which(dat$step > sl_99) 
dat_rem <- dat[rem_sl,]



#  Format data for use in model with parameters estaimted per index
mod_dat <- dat %>%
	dplyr::filter(!is.na(step)) %>%
	mutate(bin100_ship_whale_dist = round(ship_whale_dist, -2)) %>%
	mutate(bin1000_ship_whale_dist = round(ship_whale_dist, -3)) %>%
	# mutate(dist_ind = ifelse(ship_whale_dist < 667, 1, 
		# ifelse(ship_whale_dist >= 667 & ship_whale_dist < 1001, 2, 3))) %>%
	# mutate(bear_ind = ifelse(abs(ship_whale_bearing) < 15, 1, 
		# ifelse(abs(ship_whale_bearing) >= 15 & abs(ship_whale_bearing) < 30, 2, 3))) %>%
	# mutate(beh_ind = ifelse(whale_behavior == "BL-Blowing" | whale_behavior == "DF-Dive0fluke-up" |
		# whale_behavior == "DN-Dive-no-fluke", 1, 2))  %>%
	mutate(dist_ind = ifelse(ship_whale_dist < 1001, 1, 2)) %>%
	mutate(bear_ind = ifelse(abs(ship_whale_bearing) < 22.6, 1, 2)) %>%
	mutate(beh_ind = ifelse(whale_behavior == "BL-Blowing" | whale_behavior == "DF-Dive0fluke-up" |
		whale_behavior == "DN-Dive-no-fluke", 1, 2))  #%>%
	# mutate(comb_ind = ifelse(dist_ind == 1 & bear_ind == 1, 1, # close and front
		# ifelse(dist_ind == 1 & bear_ind == 3, 2, # close and side
		# ifelse(dist_ind == 3 & bear_ind == 1, 3, # far and front
		# ifelse(dist_ind == 3 & bear_ind == 3, 4,  NA))))) # far and side
	


n_pts <- nrow(mod_dat)
ind <- mod_dat$ID_new
nind <- length(unique(mod_dat$ID_new))

l <- (mod_dat$step)/1000
theta <- mod_dat$turn
dist_ind <- mod_dat$dist_ind
bear_ind <- mod_dat$bear_ind
comb_ind <- mod_dat$comb_ind
beh_ind <- mod_dat$beh_ind


kp_dist <- which(dist_ind == 1 | dist_ind == 3)
dist_ind_sm <- dist_ind[kp_dist]
dist_ind_sm[dist_ind_sm == 3] <- 2
l_dist <- l[kp_dist]
theta_dist <- theta[kp_dist]
n_pts_dist <- length(dist_ind_sm)


kp_bear <- which(bear_ind !=2)
bear_ind_sm <- bear_ind[kp_bear]
bear_ind_sm[bear_ind_sm == 3] <- 2
l_bear <- l[kp_bear]
theta_bear <- theta[kp_bear]
n_pts_bear <- length(bear_ind_sm)



kp_comb <- which(!is.na(comb_ind))
comb_ind_sm <-comb_ind[kp_comb]
l_comb <- l[kp_comb]
theta_comb <- theta[kp_comb]
n_pts_comb <- length(comb_ind_sm)





l_dist <- l
theta_dist <- theta
n_pts_dist <- length(dist_ind_sm)


kp_bear <- which(bear_ind !=2)
bear_ind_sm <- bear_ind[kp_bear]
bear_ind_sm[bear_ind_sm == 3] <- 2
l_bear <- l[kp_bear]
theta_bear <- theta[kp_bear]
n_pts_bear <- length(bear_ind_sm)
