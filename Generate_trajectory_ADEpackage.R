#  Sara Williams
#  7/14/2015
#  Script to get trajectories from repeated observations of same whale and to summarize
#  step lengths and turning angles 
################################################################################

#  Load packages
library(stats)
library(reshape)
library(plyr)
library(dplyr)
library(adehabitatLT)			
library(CircStats)
library(ggplot2)			

#  Read in whale location data (all whale observations with duplicated sightings at behvaior change or regular dataset).
temp <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/all_locs_w_dup_last_row_rem.csv")
temp <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/all_locs_w_dup.csv")
temp <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/all_locs_not_dup.csv")

#  Only use observations where there is more than one observation for comparison of behaviors
dat <- temp %>%
			 group_by(SwB_Wpt_ID)%>%
			 filter(n() > 1) %>% 
			 arrange(SwB_Wpt_ID, ObOrder_Time) %>%
			 as.data.frame(.)
			 
#  Create "transit" and "stationary" data sets.
target <- c(1, 2, 3)
transit <- filter(dat, new_beh %in% target)

target <- c(4, 5, 6)
station <- filter(dat, new_beh %in% target)
 
#  Create data sets for each behavior.
blow <- filter(dat, new_beh == 1)

target <- c(2)
dive <-filter(dat, new_beh %in% target)

target <- c(4, 5, 6)
surf <- filter(dat, new_beh %in% target)

################################################################################

########## CHOOSE one of below options for desired subset of data 
#  Arrange locations in order of the same whale (so in a group of multiple observations)
#   and then arrange also by observation order
locs2 <- arrange(dat, SwB_Wpt_ID, ObOrder_Time)
locs2 <- arrange(transit, SwB_Wpt_ID, ObOrder_Time)
locs2 <- arrange(station, SwB_Wpt_ID, ObOrder_Time)
locs2 <- arrange(dive, SwB_Wpt_ID, ObOrder_Time)
locs2 <- arrange(surf, SwB_Wpt_ID, ObOrder_Time)

locs3<- dplyr::select(locs2, SwB_Wpt_ID, whale_easting, whale_northing)
n_uni <- length(unique(locs3$SwB_Wpt_ID))
locs_tmp <- locs3 %>% 
						distinct(SwB_Wpt_ID) %>%
						mutate(Name = seq(1,n_uni))
locs_tmp2 <- full_join(locs3, locs_tmp, by = "SwB_Wpt_ID")				
locs <- locs_tmp2 %>%
			  dplyr::select(Name, whale_easting.x, whale_northing.x) %>%
			  dplyr::rename(X = whale_easting.x, Y = whale_northing.x)
################################################################################
			
########## CHOOSE one of below options for desired subset of data
#   Create ltraj object				
whale_traj <- as.ltraj(xy = locs[,c("X","Y")], id = locs$Name, typeII = FALSE)
whale_traj_transit <- as.ltraj(xy = locs[,c("X","Y")], id = locs$Name, typeII = FALSE)
whale_traj_station <- as.ltraj(xy = locs[,c("X","Y")], id = locs$Name, typeII = FALSE)
whale_traj_dive <- as.ltraj(xy = locs[,c("X","Y")], id = locs$Name, typeII = FALSE)
whale_traj_surf <- as.ltraj(xy = locs[,c("X","Y")], id = locs$Name, typeII = FALSE)
#   Convert into dataframe
traj_dat <- ld(whale_traj)
traj_dat_transit <- ld(whale_traj_transit)
traj_dat_station <- ld(whale_traj_station)
traj_dat_dive <- ld(whale_traj_dive)
traj_dat_surf <- ld(whale_traj_surf)
#   Connect traj_dat to dataframe with SwB_Wpt_ID
traj_dat_full <- cbind(traj_dat, locs_tmp2)
################################################################################

#  Create vectors for data remove NA's
#   Step lengths
steps <- as.vector(traj_dat$dist)
steps <- na.omit(steps)
steps <- as.numeric(steps)
		
steps.transit <- as.vector(traj_dat_transit$dist)
steps.transit <- na.omit(steps.transit)
steps.transit <- as.numeric(steps.transit)
		
steps.station <- as.vector(traj_dat_station$dist)
steps.station <- na.omit(steps.station)
steps.station <- as.numeric(steps.station)
		
steps.dive <- as.vector(traj_dat_dive$dist)
steps.dive <- na.omit(steps.dive)
steps.dive <- as.numeric(steps.dive)
		
steps.surf <- as.vector(traj_dat_surf$dist)
steps.surf <- na.omit(steps.surf)
steps.surf <- as.numeric(steps.surf)
		
#   Turn angles
turns <- as.vector(traj_dat$rel.angle)
turns <- na.omit(turns)
turns <- as.numeric(turns)
		
turns.transit <- as.vector(traj_dat_transit$rel.angle)
turns.transit <- na.omit(turns.transit)
turns.transit <- as.numeric(turns.transit)
		
turns.station <- as.vector(traj_dat_station$rel.angle)
turns.station <- na.omit(turns.station)
turns.station <- as.numeric(turns.station)
		
turns.dive <- as.vector(traj_dat_dive$rel.angle)
turns.dive <- na.omit(turns.dive)
turns.dive <- as.numeric(turns.dive)
		
turns.surf <- as.vector(traj_dat_surf$rel.angle)
turns.surf <- na.omit(turns.surf)
turns.surf <- as.numeric(turns.surf)

turns.station.deg <- deg(turns.station)
turns.station.deg <- turns.station.deg[turns.station.deg < 0]  + 360
turns.transit.deg <- deg(turns.transit)
turns.transit.deg <- turns.transit.deg[turns.transit.deg < 0]  + 360
turns.dive.deg <- deg(turns.dive)
turns.dive.deg <- turns.dive.deg[turns.dive.deg < 0]  + 360
turns.surf.deg <-deg(turns.surf)
turns.surf.deg <- turns.surf.deg[turns.surf.deg < 0]  + 360


#  Displacement
dis_transit <- cbind(transit, traj_dat_transit)
dis_transit <- dis_transit %>%
							dplyr::group_by(SwB_Wpt_ID) %>%
							summarise(max_dis = max(R2n)) %>%
							mutate(displace = sqrt(max_dis)/1000)

dis_station <- cbind(station, traj_dat_station)
dis_station <- dis_station %>%
							dplyr::group_by(SwB_Wpt_ID) %>%
							summarise(max_dis = max(R2n)) %>%
							mutate(displace = sqrt(max_dis)/1000)
		
dis.transit <- as.vector(dis_transit$displace)
dis.transit <- na.omit(dis.transit)
dis.transit <- as.numeric(dis.transit)
		
dis.station <- as.vector(dis_station$displace)
dis.station <- na.omit(dis.station)
dis.station <- as.numeric(dis.station)

dis_dive <- cbind(dive, traj_dat_dive)
dis_dive <- dis_dive %>%
						dplyr::group_by(SwB_Wpt_ID) %>%
						summarise(max_dis = max(R2n)) %>%
						mutate(displace = sqrt(max_dis)/1000)

dis_surf <- cbind(surf, traj_dat_surf)
dis_surf <- dis_surf %>%
					  dplyr::group_by(SwB_Wpt_ID) %>%
					  summarise(max_dis = max(R2n)) %>%
					  mutate(displace = sqrt(max_dis)/1000)

dis.dive <- as.vector(dis_dive$displace)
dis.dive <- na.omit(dis.dive)
dis.dive <- as.numeric(dis.dive)
		
dis.surf <- as.vector(dis_surf$displace)
dis.surf <- na.omit(dis.surf)
dis.surf <- as.numeric(dis.surf)
		
#   Number of "relocations"
n_locs_transit <- cbind(transit, traj_dat_transit)
n_locs_transit <- n_locs_transit %>%
								   dplyr::group_by(SwB_Wpt_ID) %>%
								   summarise(num_locs = max(ObOrder_Time)) %>%
								   filter(num_locs > 1)

n_locs_station <- cbind(station, traj_dat_station)
n_locs_station <- n_locs_station %>%
									dplyr::group_by(SwB_Wpt_ID) %>%
									summarise(num_locs = max(ObOrder_Time)) %>%
									filter(num_locs > 1)
										
nlocs.transit <- as.vector(n_locs_transit$num_locs)
nlocs.transit <- na.omit(nlocs.transit)
nlocs.transit <- as.numeric(nlocs.transit)
		
nlocs.station <- as.vector(n_locs_station$num_locs)
nlocs.station <- na.omit(nlocs.station)
nlocs.station <- as.numeric(nlocs.station)

n_locs_dive <- cbind(dive, traj_dat_dive)
n_locs_dive <- n_locs_dive %>%
							  dplyr::group_by(SwB_Wpt_ID) %>%
							  summarise(num_locs = max(ObOrder_Time)) %>%
							  filter(num_locs > 1)

n_locs_surf <- cbind(surf, traj_dat_surf)
n_locs_surf <- n_locs_surf %>%
							 dplyr::group_by(SwB_Wpt_ID) %>%
							 summarise(num_locs = max(ObOrder_Time)) %>%
							 filter(num_locs > 1)
										
nlocs.dive <- as.vector(n_locs_dive$num_locs)
nlocs.dive <- na.omit(nlocs.dive)
nlocs.dive <- as.numeric(nlocs.dive)
		
nlocs.surf <- as.vector(n_locs_surf$num_locs)
nlocs.surf <- na.omit(nlocs.surf)
nlocs.surf <- as.numeric(nlocs.surf)
################################################################################

#  Plots
#  Step lengths
steps_hist_transit <- ggplot(traj_dat_transit, aes(x=dist)) + 
										   geom_histogram(aes(y=..density..),
										   color="black", fill="grey", 
										   binwidth=100) +
										   geom_density(alpha=.2, fill="#FF6666") +
										   xlim(0, 4000) +
										   ylim(0,0.003) +
										   theme_bw()
steps_hist_transit	

steps_hist_station <- ggplot(traj_dat_station, aes(x=dist)) + 
											geom_histogram(aes(y=..density..),
											color="black", fill="grey", 
											binwidth=100) +
											geom_density(alpha=.2, fill="#FF6666") +
											xlim(0, 4000) +
											ylim(0,0.003) +
											theme_bw()
steps_hist_station

steps_hist_dive <- ggplot(traj_dat_dive, aes(x=dist)) + 
									  geom_histogram(aes(y=..density..),
									  color="black", fill="grey", 
									  binwidth=100) +
									  geom_density(alpha=.2, fill="#FF6666") +
									  xlim(0, 4000) +
									  ylim(0,0.003) +
									  theme_bw()
steps_hist_dive	

steps_hist_surf <- ggplot(traj_dat_surf, aes(x=dist)) + 
									 geom_histogram(aes(y=..density..),
									 color="black", fill="grey", 
									 binwidth=100) +
									 geom_density(alpha=.2, fill="#FF6666") +
									 xlim(0, 4000) +
									 ylim(0,0.003) +
									 theme_bw()
steps_hist_surf
	
#   Create dataframe of observed dive steps and observed stationary steps together.
#   Have to have objects of same length, so add on NAs to shorter vector.
max.len <- max(length(steps.dive), length(steps.station))
steps.station <- c(steps.station, rep(NA, max.len - length(steps.station)))
obs_steps <- as.data.frame(cbind(steps.dive, steps.station))
names(obs_steps)[1] <- "dist_dive"
names(obs_steps)[2] <- "dist_station"

#   Look at density plots of both on same figure
df_steps_obs <- melt(obs_steps)
compare_steps_obs <- ggplot(df_steps_obs) + 
											   geom_density(aes(x = value, colour = variable)) +
											   xlim(0, 4000) +
											   theme_bw()
compare_steps_obs

#   Create dataframe of observed dive steps and observed surface steps
#   Have to have objects of same length, so add on NAs to shorter vector.
max.len <- max(length(steps.dive), length(steps.surf))
steps.surf<- c(steps.surf, rep(NA, max.len - length(steps.surf)))
obs_steps_dive_surf <- as.data.frame(cbind(steps.dive, steps.surf))
names(obs_steps_dive_surf )[1] <- "dist_dive"
names(obs_steps_dive_surf )[2] <- "dist_surf"

#   Look at density plots of both on same figure
df_steps_obs_dive_surf <- melt(obs_steps_dive_surf)
compare_steps_obs_dive_surf <- ggplot(df_steps_obs_dive_surf ) + 
																	 geom_density(aes(x = value, colour = variable)) +
																	 xlim(0, 4000) +
																	 theme_bw()
compare_steps_obs_dive_surf 

#  Turn angles
turns_hist_transit <- ggplot(traj_dat_transit, aes(x=rel.angle)) + 
										   geom_histogram(aes(y=..density..),
										   color="black", fill="grey", 
										   binwidth=.1) +
										   geom_density(alpha=.2, fill="#FF6666") +
										   xlim(-3.5, 3.5) +
										   ylim(0,1.0) +
										   theme_bw()
turns_hist_transit	

turns_hist_station <- ggplot(traj_dat_station, aes(x=rel.angle)) + 
										   geom_histogram(aes(y=..density..),
										   color="black", fill="grey", 
										   binwidth=.1) +
										   geom_density(alpha=.2, fill="#FF6666") +
										   xlim(-3.5, 3.5) +
										   ylim(0,1.0) +
										   theme_bw()
turns_hist_station	

turns_hist_dive <- ggplot(traj_dat_dive, aes(x=rel.angle)) + 
									  geom_histogram(aes(y=..density..),
									  color="black", fill="grey", 
									  binwidth=.1) +
									  geom_density(alpha=.2, fill="#FF6666") +
									  xlim(-3.5, 3.5) +
									  ylim(0,1.0) +
									  theme_bw()
turns_hist_dive	

turns_hist_surf <- ggplot(traj_dat_surf, aes(x=rel.angle)) + 
									 geom_histogram(aes(y=..density..),
									 color="black", fill="grey", 
									 binwidth=.1) +
									 geom_density(alpha=.2, fill="#FF6666") +
									 xlim(-3.5, 3.5) +
									 ylim(0,1.0) +
									 theme_bw()
turns_hist_surf	

#   Create dataframe of observed dive turns and observed stationary turns
#   Have to have objects of same length, so add on NAs to shorter vector.
max.len <- max(length(turns.dive), length(turns.station))
turns.station<- c(turns.station, rep(NA, max.len - length(turns.station)))
obs_turns <- as.data.frame(cbind(turns.dive, turns.station))
names(obs_turns)[1] <- "turn_dive"
names(obs_turns)[2] <- "turn_station"
	
#   Look at density plots of both on same figure
df_turns_obs <- melt(obs_turns)
compare_turns_obs <- ggplot(df_turns_obs) + 
											   geom_density(aes(x = value, colour = variable)) +
										       xlim(-4, 4) +
											   theme_bw()
compare_turns_obs	

#   Create dataframe of observed transit turns and observed stationary turns IN DEGRESS
#   Have to have objects of same length, so add on NAs to shorter vector.
max.len <- max(length(turns.dive.deg), length(turns.station.deg))
turns.station.deg<- c(turns.station.deg, rep(NA, max.len - length(turns.station.deg)))
obs_turns_deg <- as.data.frame(cbind(turns.dive.deg, turns.station.deg))
names(obs_turns_deg)[1] <- "turn_dive"
names(obs_turns_deg)[2] <- "turn_station"
	
#   Look at density plots of both on same figure
df_turns_obs_deg <- melt(obs_turns_deg)
compare_turns_obs_deg <- ggplot(df_turns_obs_deg) + 
														geom_density(aes(x = value, colour = variable)) +
														xlim(0, 360) +
														theme_bw()
compare_turns_obs_deg	
	
#   Create dataframe of observed dive turns and observed surface turns
#   Have to have objects of same length, so add on NAs to shorter vector.
max.len <- max(length(turns.dive), length(turns.surf))
turns.surf<- c(turns.surf, rep(NA, max.len - length(turns.surf)))
obs_turns_dive_surf <- as.data.frame(cbind(turns.dive, turns.surf))
names(obs_turns_dive_surf)[1] <- "turn_dive"
names(obs_turns_dive_surf)[2] <- "turn_surf"

#   Look at density plots of both on same figure
df_turns_obs_dive_surf <- melt(obs_turns_dive_surf)
compare_turns_obs_dive_surf <- ggplot(df_turns_obs_dive_surf) + 
																	geom_density(aes(x = value, colour = variable)) +
																	xlim(-3.2, 3.2) +
																	theme_bw()
compare_turns_obs_dive_surf

#  Displacement
#   Create dataframe of observed transit total displacement and observed stationary displacement
#   Have to have objects of same length, so add on NAs to shorter vector.
max.len <- max(length(dis.transit), length(dis.station))
dis.station <- c(dis.station, rep(NA, max.len - length(dis.station)))
obs_dis <- as.data.frame(cbind(dis.transit, dis.station))
names(obs_dis)[1] <- "displace_transit"
names(obs_dis)[2] <- "dispalce_station"
#  Look at density plots of both on same figure
df_dis_obs <- melt(obs_dis)
compare_dis_obs <- ggplot(df_dis_obs) + 
										  geom_density(aes(x = value, colour = variable)) +
										  xlim(0, 10) +
										  theme_bw()
compare_dis_obs		

#   Create dataframe of observed dive total displacement and observed surface displacement
#   Have to have objects of same length, so add on NAs to shorter vector.
max.len <- max(length(dis.dive), length(dis.surf))
dis.surf <- c(dis.surf, rep(NA, max.len - length(dis.surf)))
obs_dis_dive_surf <- as.data.frame(cbind(dis.dive, dis.surf))
names(obs_dis_dive_surf)[1] <- "displace_dive"
names(obs_dis_dive_surf)[2] <- "dispalce_surf"

#   Look at density plots of both on same figure
df_dis_obs_dive_surf <- melt(obs_dis_dive_surf)
compare_dis_obs_dive_surf <- ggplot(df_dis_obs_dive_surf) + 
															   geom_density(aes(x = value, colour = variable)) +
															   xlim(0, 6) +
															   theme_bw()
compare_dis_obs_dive_surf		

#  Number of "relocations"
#   Create dataframe of observed transit relocations and observed stationary relocations
#   Have to have objects of same length, so add on NAs to shorter vector.
max.len <- max(length(nlocs.transit), length(nlocs.station))
nlocs.station <- c(nlocs.station, rep(NA, max.len - length(nlocs.station)))
obs_nlocs <- as.data.frame(cbind(nlocs.transit, nlocs.station))
names(obs_nlocs)[1] <- "nlocs_transit"
names(obs_nlocs)[2] <- "nlocs_station"

#   Look at density plots of both on same figure
df_nlocs_obs <- melt(obs_nlocs)
compare_nlocs_obs <- ggplot(df_nlocs_obs) + 
											   geom_density(aes(x = value, colour = variable)) +
										       xlim(0, 10) +
											   theme_bw()
compare_nlocs_obs		
		
#   Create dataframe of observed dive relocations and observed surf relocations
#   Have to have objects of same length, so add on NAs to shorter vector.
max.len <- max(length(nlocs.dive), length(nlocs.surf))
nlocs.surf <- c(nlocs.surf, rep(NA, max.len - length(nlocs.surf)))
obs_nlocs_dive_surf <- as.data.frame(cbind(nlocs.dive, nlocs.surf))
names(obs_nlocs_dive_surf )[1] <- "nlocs_dive"
names(obs_nlocs_dive_surf )[2] <- "nlocs_surf"

#   Look at density plots of both on same figure
df_nlocs_obs_dive_surf  <- melt(obs_nlocs_dive_surf )
compare_nlocs_obs_dive_surf  <- ggplot(df_nlocs_obs_dive_surf ) + 
																	  geom_density(aes(x = value, colour = variable)) +
																	  xlim(0, 10) +
																	  theme_bw()
compare_nlocs_obs_dive_surf 				

# # ##################################################################################################
# # ######## Example code from adehabitatLT vignette ####################################################
# # data(puechabonsp)
# # head(puechabonsp)
# # locs <- puechabonsp$relocs
# # locs <- as.data.frame(locs)
# # head(locs)

# # puech <- as.ltraj(xy = locs[,c("X","Y")], date = da, id = locs$Name)
# # puech
# # head(puech[[1]])
# # plot(puech)
# # ##################################################################################################