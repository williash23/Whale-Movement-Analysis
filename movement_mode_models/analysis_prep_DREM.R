# Sara Williams
# 9/18/2019
# Data prep for DREM 
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
	dplyr::rename(X = X_whale_UTM, Y = Y_whale_UTM) %>%
	dplyr::filter(ship_whale_dist < 5001)
	
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

#  Convert into dataframe
traj_df_tmp <- ld(traj)
traj_df <- traj_df_tmp %>%
                 dplyr::select(id, x, y, dist, rel.angle, abs.angle, dx, dy) 
names(traj_df) <- c("same_whale_ID", "X", "Y", "step", "turn", "abs_angle", "dx", "dy")

#  Join trajectory dataframe with dataframe holding location variables
mod_dat_tmp <- full_join(traj_df, tmp4, by = c("same_whale_ID", "X", "Y")) %>%
	mutate(ave_swim_spd = step/time_diff_sec)
#  Remove erroneous steps based on 90th or 95th quantile
#quantile(mod_dat_tmp$step, probs = c(0.90, 0.95), na.rm = TRUE)
#quantile(mod_dat_tmp$ave_swim_spd, probs = c(0.90, 0.95), na.rm = TRUE)
#mod_dat <- filter(mod_dat_tmp, step < 2086 | is.na(step))
mod_dat <- mod_dat_tmp %>%
	dplyr::filter(ave_swim_spd < 6.6877) %>% #4.361111
	dplyr::filter(!is.na(step))



################################################################################

#  Generate data for use in DREM 
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
l_double <- as.matrix(l_double)
l_double <- l_double[,-1]	

theta_double <- obs_1 %>%
	dplyr::select(ID_new, occ, turn) %>%
	spread(occ, turn, fill = NA, convert = FALSE)
theta_double <- as.matrix(theta_double)
theta_double <- theta_double[,-1]	

dist_ind <- obs_1 %>%
	mutate(DREM_dist_ind = ifelse(ship_whale_dist < 1001, 1, 
		ifelse(ship_whale_dist >= 1001 & ship_whale_dist < 3001, 2, 3))) %>%
	dplyr::select(ID_new, occ, DREM_dist_ind) %>%
	spread(occ, DREM_dist_ind, fill = NA, convert = FALSE)
dist_ind <- as.matrix(dist_ind)
dist_ind <- dist_ind[,-1]

bear_ind <- obs_1 %>%
	mutate(DREM_bear_ind = ifelse(abs(ship_whale_bearing) < 22.56, 1,
		ifelse(abs(ship_whale_bearing) >= 22.56 & abs(ship_whale_bearing) < 45, 2, 3))) %>%
	dplyr::select(ID_new, occ, DREM_bear_ind) %>%
	spread(occ, DREM_bear_ind, fill = NA, convert = FALSE)
bear_ind <- as.matrix(bear_ind)
bear_ind <- bear_ind[,-1]

dur_ind <- obs_1 %>%
	mutate(DREM_dur_ind = ifelse(time_diff_sec < 50, 1,
		ifelse(time_diff_sec >= 50 & time_diff_sec < 120, 2, 3))) %>%
	dplyr::select(ID_new, occ, DREM_dur_ind) %>%
	spread(occ, DREM_dur_ind, fill = NA, convert = FALSE)
dur_ind <- as.matrix(dur_ind)
dur_ind <- dur_ind[,-1]



#  Plotting
# load API

#  Load data
dat_raw <- read.csv("C:/Users/saraw/Documents/Whales/data/Whales_0615_general_clean.csv")

#  Data clean up and manipulation
tmp1 <- dat_raw %>%
	filter(!is.na(X)) %>%
	as.data.frame()
tmp2 <- arrange(tmp1, same_whale_ID, ob_order_time)
tmp3 <- tmp2 %>%
	dplyr::select(same_whale_ID, X_whale_UTM, Y_whale_UTM, ship_whale_dist, 
		ship_whale_bearing, whale_behavior, ob_order_time, TimeTxt, same_whale_ID,
		year, month) %>%
	dplyr::rename(X = X_whale_UTM, Y = Y_whale_UTM) %>%
	filter(!is.na(X)) %>%
	filter(!is.na(Y))

dat <- st_as_sf(tmp3, coords = c("X", "Y"), crs = 32608) %>%
	st_transform(4326) 
lon <- data.frame(lon = st_coordinates(dat)[,1])
lat <- data.frame(lat = st_coordinates(dat)[,2])
dat <- cbind(dat, lon$lon)
dat <- cbind(dat, lat$lat)

center = c(lon = -135.9865722, lat = 58.5854358)


p <- qmap(center, zoom = 8, source = "stamen", maptype = "toner-lite")  +
	geom_point(data = dat, size = 0.75,
		aes(x = lon.lon, y = lat.lat, colour = ship_whale_dist)) +
	scale_colour_gradient(name = "Dist. to Ship", low = "#F3DF6C", high = "#B40F20") +
	theme_bw() +
	xlab("Longitude") +
	ylab("Latitude") +
	theme(text = element_text(family = "Calibri Light", size = 12),
		axis.title.y = element_text(face = "bold", angle = 90, vjust = 2),
		axis.title.x = element_text(face = "bold", vjust = -0.2),
		axis.text.y = element_text(face = "bold"), 
		axis.text.x = element_text(angle = 90, vjust = 0.5), 
		axis.ticks = element_line(),
		panel.grid.major = element_line(colour="#f0f0f0"),
		panel.grid.minor = element_blank(),
		legend.title = element_text(face = "bold"),
		plot.margin=unit(c(10,5,5,5),"mm"),
		panel.border = element_rect(colour = "black", fill = NA, size = 1),
		legend.background = element_blank(),
		legend.box.background = element_rect(colour = "black", fill = "white"),
		legend.position = c(0.125, 0.175))
p


p <- ggmap(get_googlemap(center = c(lon = -135.9865722, lat = 58.5854358),
	zoom = 8, scale = 2,
	maptype ='terrain',
	color = 'bw'))
p + geom_point(aes(x = lon.lon, y = lat..1.,  colour = Initial.Type.Group), data = i2, size = 0.5) + 
  theme(legend.position="bottom")


p2 <- ggplot() + 
	geom_histogram(data = dat, bins = 1500, aes(ship_whale_dist, fill = factor(year)))
p2


, position = "dodge"

p2 <- ggplot() + 
	geom_histogram(data = dat, aes(x = factor(year), y = ship_whale_dist, fill = factor(month)),
	stat = "identity", position = "dodge")
p2










