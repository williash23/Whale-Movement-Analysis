
#  Load packages
library(dplyr)
library(adehabitatLT)    
library(tidyr)
library(stringr)
library(aspace)
################################################################################

#  Load data
dat_raw <- read.csv("C:/Users/saraw/Documents/UM_current_work/GitHub/Whale-Movement-Analysis/data/Whales_0615_locations_clean.csv")

#  Data clean up and manipulation
tmp1 <- dat_raw %>%
	filter(year > 2009) %>% ## per info from Karin, PTB (bearing) was kind of weird before 2010
	filter(ship_whale_dist < 8000) %>%
	filter(count == 1) %>%            
	group_by(same_whale_ID) %>%
	filter(n() > 1) %>%
	ungroup() %>%
	as.data.frame()
tmp2 <- arrange(tmp1, same_whale_ID, ob_order_time)
tmp3 <- tmp2 %>%
	dplyr::select(same_whale_ID, X_whale_UTM, Y_whale_UTM, ship_whale_dist, 
		ship_whale_bearing, whale_orientation, whale_direction, ob_order_time, TimeTxt) %>%
	dplyr::rename(X = X_whale_UTM, Y = Y_whale_UTM) %>%
	mutate(sight_side = ifelse(ship_whale_bearing > 0, "starb", "port"))
tmp4 <- tmp3 %>%
	group_by(same_whale_ID) %>%
	mutate(switch_sides = ifelse(lag(sight_side) == sight_side, 0, 1)) %>%
	as.data.frame()
switchers_id <- tmp4 %>%
	group_by(same_whale_ID) %>%
	summarise(switch_sum = sum(switch_sides, na.rm = TRUE)) %>%
	as.data.frame()  %>%
	dplyr::filter(switch_sum > 0)
switchers_dat <- inner_join(tmp4, switchers_id, by = "same_whale_ID")
non_switcher_dat <- anti_join(tmp4, switchers_id, by = "same_whale_ID")
port_dat_tmp <- non_switcher_dat %>%
	dplyr::filter(sight_side == "port") %>%
	dplyr::select(-whale_orientation)
starb_dat_tmp <- non_switcher_dat %>%
	dplyr::filter(sight_side == "starb") %>%
	dplyr::select(-whale_orientation)

port_dat_rep <- port_dat_tmp
port_dat_rep$whale_direction <- as.character(port_dat_rep $whale_direction)
port_dat_rep[port_dat_rep == "U-Unsure"] <- NA 
port_dat_rep[port_dat_rep == "U - Unsure"] <- NA
port_dat_rep[port_dat_rep == "W-With"] <- 0
port_dat_rep[port_dat_rep == "W/A-Wi/Aw"] <- 45
port_dat_rep[port_dat_rep == "WA-With/Aw"] <- 45
port_dat_rep[port_dat_rep == "A-Away"] <- 90
port_dat_rep[port_dat_rep == "GA-Ag/Away"] <- 135
port_dat_rep[port_dat_rep == "G/A-Ag/Aw"] <- 135
port_dat_rep[port_dat_rep == "G-Against"] <- 180
port_dat_rep[port_dat_rep == "G/T-Ag/to"] <- 225
port_dat_rep[port_dat_rep == "GT-Ag/To"] <- 225
port_dat_rep[port_dat_rep == "T-Toward"] <- 270
port_dat_rep[port_dat_rep == "W/T-Wi/T0"] <- 315
port_dat_rep[port_dat_rep == "W/T-Wit/To"] <- 315
port_dat_rep[port_dat_rep == "WT-With/To"] <- 315
port_dat <- port_dat_rep
port_dat$whale_direction <- as.numeric(port_dat$whale_direction)

starb_dat_rep <- starb_dat_tmp
starb_dat_rep$whale_direction <- as.character(starb_dat_rep $whale_direction)
starb_dat_rep[starb_dat_rep == "U-Unsure"] <- NA 
starb_dat_rep[starb_dat_rep == "U - Unsure"] <- NA
starb_dat_rep[starb_dat_rep == "W-With"] <- 0
starb_dat_rep[starb_dat_rep == "W/A-Wi/Aw"] <- 45
starb_dat_rep[starb_dat_rep == "WA-With/Aw"] <- 45
starb_dat_rep[starb_dat_rep == "A-Away"] <- 90
starb_dat_rep[starb_dat_rep == "GA-Ag/Away"] <- 135
starb_dat_rep[starb_dat_rep == "G/A-Ag/Aw"] <- 135
starb_dat_rep[starb_dat_rep == "G-Against"] <- 180
starb_dat_rep[starb_dat_rep == "G/T-Ag/to"] <- 225
starb_dat_rep[starb_dat_rep == "GT-Ag/To"] <- 225
starb_dat_rep[starb_dat_rep == "T-Toward"] <- 270
starb_dat_rep[starb_dat_rep == "W/T-Wi/T0"] <- 315
starb_dat_rep[starb_dat_rep == "W/T-Wit/To"] <- 315
starb_dat_rep[starb_dat_rep == "WT-With/To"] <- 315
starb_dat <- starb_dat_rep
starb_dat$whale_direction <- as.numeric(starb_dat$whale_direction)

port_dat <- port_dat %>%
	group_by(same_whale_ID) %>%
	mutate(sight_num = 1:n()) %>%
	mutate(turn = abs(whale_direction - lag(whale_direction))) %>%
	mutate(step = abs(
	as.data.frame()
starb_dat <- starb_dat %>%
	group_by(same_whale_ID) %>%
	mutate(sight_num = 1:n()) %>%
	mutate(turn = abs(whale_direction - lag(whale_direction))) %>%
	as.data.frame()
################################################################################

#  Generate data for use in model runs
#   Port
obs_1 <- port_dat %>%
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
theta <- obs_1$turn

#   Starboard
obs_1 <- starb_dat %>%
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
theta <- obs_1$turn