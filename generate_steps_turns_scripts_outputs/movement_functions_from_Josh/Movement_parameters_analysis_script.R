# Sara Williams
# 7/14/2015
# Analysis script to run movement parameter generation functions
################################################################################

#  Load packages
library(stats)
library(reshape)
library(plyr)
library(dplyr)
library(adehabitatLT)			
library(CircStats)
library(ggplot2)			

#  Read in whale location data (all whale observations with duplicated sightings at behvaior change).
temp <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/Whales_0615_w_dup_last_row_rem.csv")

#  Arrange locations in order of the same whale (so in a group of multiple observations)
#   and by observation order. 
#  Create dataset of all points with no duplicated points. Duplication is unncessary when not dividing by behavior.
#  Use only whales that were observed at least 2 times.
dat <- temp %>%
			 group_by(same_whale_ID) %>%
			 filter(n() > 1) %>%
			 ungroup(.) %>%
			 arrange(same_whale_ID, ob_order_time) %>%
			 as.data.frame(.)
		 
#  Create behavior type data sets.
target <- c(1, 2, 3)
transit <- filter(dat, new_beh %in% target)

target <- c(4, 5, 6)
station <- filter(dat, new_beh %in% target)

blow <- filter(dat, new_beh == 1)

target <- c(2)
dive <-filter(dat, new_beh %in% target)

target <- c(4, 5, 6)
surf <- filter(dat, new_beh %in% target)
#############################################################################
########## CHOOSE one of below options for desired subset of data #################

pts1 <- dat
pts1 <- transit
pts1 <- station

pts1 <- blow
pts1 <- dive
pts1 <- surf

# #   Instead, using only points that are, before, at ,or after CPA
# source(file.path("C:\Users\sara.williams\Documents\GitHub\Whale-Behavior-Analysis\exploratory_behavior_change.R"))
# pts1 <- before
# pts1 <- cpa
# pts1 <- after

#  Source Movement_parameters_fun
source(file.path("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/Move_param_fun.R"))

#  Prep locs.data: this is one of the inputs for the move_param_fun
#   Arrange by same whale ID and observation order. Also remove pts where distance between 
#   whale and shore is less than 0 (erroneous?)
pts2 <- pts1 %>% 							
			   dplyr::select(same_whale_ID, X_whale_UTM, Y_whale_UTM) %>%
			   dplyr::rename(x = X_whale_UTM , y = Y_whale_UTM) 
							
loc.data <- as.data.frame(pts2)
loc.data.transit <- as.data.frame(pts2)
loc.data.station <- as.data.frame(pts2)
loc.data.dive <- as.data.frame(pts2)
loc.data.surf <- as.data.frame(pts2)
		
		#  Location data processing: go through data set by "same_whale_ID", select on whales that have more than 
		#   2 sightings (need at least 3 sightings for step length and turning angle).		
		#  Run move.fun!!!
		data.out <- loc.data %>%
								group_by(SwB_Wpt_ID)%>%
								filter(n()>2) %>%
								do(move.fun(data.frame(x = .$x, y = .$y))) %>%
								as.data.frame(.)
									
		data.out.transit <- loc.data.transit %>%
												group_by(SwB_Wpt_ID)%>%
												#filter(n()>2) %>%
												do(move.fun(data.frame(x = .$x, y = .$y))) %>%
												as.data.frame(.)

		data.out.station <- loc.data.station %>%
												group_by(SwB_Wpt_ID)%>%
												filter(n()>2) %>%
												do(move.fun(data.frame(x = .$x, y = .$y))) %>%
												as.data.frame(.)
												
		data.out.dive <- loc.data.dive %>%
												group_by(same_whale_ID)%>%
												filter(n()>2) %>%
												do(move.fun(data.frame(x = .$x, y = .$y))) %>%
												as.data.frame(.)
												
		data.out.surf <- loc.data.surf %>%
												group_by(same_whale_ID)%>%
												filter(n()>2) %>%
												do(move.fun(data.frame(x = .$x, y = .$y))) %>%
												as.data.frame(.)
												
		# #  Compare histograms of step lengths and turning angles of all sightings vs sightings just within 2000m
		library(ggplot2)
		# #   Step lengths
		steps_hist_transit <- ggplot(data.out.transit, aes(x=dist)) + 
													geom_histogram(aes(y=..density..),
													color="black", fill="grey", 
													binwidth=50) +
													geom_density(alpha=.2, fill="#FF6666") +
													xlim(0, 4000) +
													ylim(0,0.003) +
													theme_bw()
		steps_hist_transit	

		steps_hist_station <- ggplot(data.out.station, aes(x=dist)) + 
													geom_histogram(aes(y=..density..),
													color="black", fill="grey", 
													binwidth=50) +
													geom_density(alpha=.2, fill="#FF6666") +
													xlim(0, 4000) +
													ylim(0,0.003) +
													theme_bw()
		steps_hist_station

		# #  Turn angles
		turns_hist_transit <- ggplot(data.out.transit, aes(x=turn.angle)) + 
									geom_histogram(aes(y=..density..),
									color="black", fill="grey", 
									binwidth=.1) +
									geom_density(alpha=.2, fill="#FF6666") +
									xlim(-3.5, 3.5) +
									ylim(0,1.0) +
									theme_bw()
		turns_hist_transit	

		turns_hist_station <- ggplot(data.out.station, aes(x=turn.angle)) + 
									geom_histogram(aes(y=..density..),
									color="black", fill="grey", 
									binwidth=.1) +
									geom_density(alpha=.2, fill="#FF6666") +
									xlim(-3.5, 3.5) +
									ylim(0,1.0) +
									theme_bw()
		turns_hist_station			

		#  Create vectors of just turns and step lengths and remove NA's
		#   Step lengths
		steps <- as.vector(data.out$dist)
		steps <- na.omit(steps)
	    steps <- as.numeric(steps)
		
		steps.transit <- as.vector(data.out.transit$dist)
		steps.transit <- na.omit(steps.transit)
	    steps.transit <- as.numeric(steps.transit)
		
		steps.station <- as.vector(data.out.station$dist)
		steps.station <- na.omit(steps.station)
	    steps.station <- as.numeric(steps.station)
		
		steps.dive <- as.vector(data.out.dive$dist)
		steps.dive <- na.omit(steps.dive)
	    steps.dive <- as.numeric(steps.dive)
		
		steps.surf <- as.vector(data.out.surf$dist)
		steps.surf <- na.omit(steps.surf)
	    steps.surf <- as.numeric(steps.surf)
		
		#   Turn angles
		turns <- as.vector(data.out$turn.angle)
		turns <- na.omit(turns)
		turns <- as.numeric(turns)
		
		turns.transit <- as.vector(data.out.transit$turn.angle)
		turns.transit <- na.omit(turns.transit)
		turns.transit <- as.numeric(turns.transit)
		
		turns.station <- as.vector(data.out.station$turn.angle)
		turns.station <- na.omit(turns.station)
		turns.station <- as.numeric(turns.station)
		
		turns.dive <- as.vector(data.out.dive$turn.angle)
		turns.dive <- na.omit(turns.dive)
		turns.dive <- as.numeric(turns.dive)
		
		turns.surf <- as.vector(data.out.surf$turn.angle)
		turns.surf <- na.omit(turns.surf)
		turns.surf <- as.numeric(turns.surf)
		
		length(steps.transit)
		range(steps.transit)
		length(turns.transit)
		range(turns.transit)
		
		length(steps.station)
		range(steps.station)
		length(turns.station)
		range(turns.station)
#############################################################################	
