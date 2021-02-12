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
librayr(lubridate)

`%!in%` = function(x,y) !(x %in% y)

################################################################################

#  Load data
dat_raw <- read.csv("C:/Users/saraw/Documents/Whales/data/Whales_0615_general_clean.csv")
sa <- st_read("C:/Users/saraw/Documents/Whales/Tracks_0615/survey_area_buff15k.shp")

#  Data clean up and manipulation
tmp <- dat_raw %>%
	dplyr::filter(year > 2009) %>% ## per info from Karin, PTB (bearing) was kind of weird before 2010
	dplyr::filter(!is.na(X_whale_UTM)) %>%
	mutate(exact = ifelse(approx == "n", 1, 0)) %>%
	group_by(same_whale_ID) %>%
	filter(n() > 1) %>%
	mutate(ever_exact = max(exact)) %>%
	#filter(approx == "n") %>%
	as.data.frame()  %>%
	arrange(same_whale_ID, ob_order_time)  %>%
	mutate(ob_time = mdy_hms(paste(DateD, TimeTxt))) %>% 
	dplyr::select(unique_event_ID, same_whale_ID, count,
		X_whale_UTM, Y_whale_UTM, approx,
		ob_order_time, ob_time, 
		ship_whale_dist,  ship_whale_bearing, 
		whale_behavior, visibility) 
		
#  Create trajectory (ltraj) object for steps and turns
tmp$same_whale_ID <- droplevels(tmp$same_whale_ID)
traj <- as.ltraj(xy = tmp[,c("X_whale_UTM","Y_whale_UTM")], 
	id = tmp$same_whale_ID, typeII = FALSE)

#  Convert into dataframe
traj_df_tmp <- ld(traj)
traj_df <- traj_df_tmp %>%
	dplyr::select(id, x, y, dist, rel.angle, abs.angle, dx, dy) 
names(traj_df) <- c("same_whale_ID", "X_whale_UTM", "Y_whale_UTM", 
	"step", "turn", "abs_angle", "dx", "dy")

#  Join trajectory dataframe with dataframe holding location variables
dat1 <- full_join(traj_df, tmp, by = c("same_whale_ID", "X_whale_UTM", "Y_whale_UTM")) %>%
	group_by(same_whale_ID) %>%
	dplyr::mutate(occ = 1:n()) %>%
	dplyr::mutate(n_occ = n()) %>%
	mutate(tdiff_s_tmp = lead(ob_time) - ob_time) %>%
	mutate(tdiff_s = ifelse(tdiff_s_tmp < 1, 1, tdiff_s_tmp)) %>%
	mutate(spd_mps = as.numeric(step)/as.numeric(tdiff_s)) %>%
	as.data.frame() %>%
	dplyr::select(-tdiff_s_tmp)

dat2 <- dat1 %>%
	dplyr::filter(count == 1) %>%
	dplyr::filter(ship_whale_dist <= 3000) %>%
	dplyr::filter(ship_whale_bearing <= 45 & ship_whale_bearing >= -45) %>%
	dplyr::filter(spd_mps <= 4.8) %>%
	dplyr::filter(tdiff_s <= 300)

# Constrain to 99% quantiles for surfacing lengths 
t95 <- as.numeric(quantile(as.numeric(dat2$tdiff_s[!is.na(dat2$step)]), probs = c(0.95)))
t99 <- as.numeric(quantile(as.numeric(dat2$tdiff_s[!is.na(dat2$step)]), probs = c(0.99)))
sl95 <- as.numeric(quantile(dat2$step[!is.na(dat2$step)], probs = c(0.95)))
sl99 <- as.numeric(quantile(dat2$step[!is.na(dat2$step)], probs = c(0.99)))

kp_sl  <- which(dat2$step <= sl95)
dat95 <- dat2[kp_sl,]



# Number of surfacings
nrow(dat2)
nrow(dat95)

# Number of unique surfacing bouts/whales
length(unique(dat2$same_whale_ID))
length(unique(dat95$same_whale_ID))

# Number of surfacing lengths
length(dat2$step[!is.na(dat2$step)])
length(dat95$step[!is.na(dat95$step)])

# Number of surfacing angles
length(dat2$turn[!is.na(dat2$turn)])
length(dat95$turn[!is.na(dat95$turn)])

t95
t99
sl95
sl99


#  Format data for use in model with parameters estaimted per index
mod_dat <- dat2 %>%
	dplyr::filter(!is.na(step)) %>%
	mutate(bin100_ship_whale_dist = round(ship_whale_dist, -2)) %>%
	mutate(bin10_ship_whale_dist = round(ship_whale_dist, -1)) %>%
	mutate(dist_ind = ifelse(ship_whale_dist <= 1500, 1, 2)) %>%
	mutate(bear_ind = ifelse(abs(ship_whale_bearing) <= 22.6, 1, 2)) %>%
	mutate(beh_ind = ifelse(whale_behavior == "BL-Blowing" | whale_behavior == "DF-Dive-fluke-up" |
		whale_behavior == "DN-Dive-no-fluke", 1, 2))  
		
n_pts <- nrow(mod_dat)
nind <- length(unique(mod_dat$same_whale_ID))

l <- (mod_dat$step)/1000
theta <- mod_dat$turn
dist_ind <- mod_dat$dist_ind
bear_ind <- mod_dat$bear_ind
comb_ind <- mod_dat$comb_ind
beh_ind <- mod_dat$beh_ind


