# Sara Williams
# 7/14/2015
# Analysis script to run movement parameter generation functions
################################################################################
		#  Load packages
		library(plyr)
		library(dplyr)
		
		#  Read in whale observation location data
		#  Raw data
		temp <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/Whale_Pts.csv")
		head(temp)
		table(temp$whale_behavior)

		temp2 <-filter(temp, temp$whale_behavior=="BL-Blowing" | temp$whale_behavior=="DF-Dive-fluke-up" | 
											temp$whale_behavior=="DN-Dive-no-fluke" | temp$whale_behavior=="LF-Lunge-feed" | 
											temp$whale_behavior=="RE-Resting" | temp$whale_behavior=="SA-Surface-active")
		temp2$whale_behavior <- droplevels(temp2$whale_behavior)

		#  Only use observations where there is more than one observation for comparison of behaviors
		temp2 <- filter(temp2, ObType=="MultiOb")

		#  Make a new whale_behavior category that is numeric and attach to dataframe
		#   1: transit type behaviors (blowing/dive with no fluke, dive with fluke up) and 2: stationary type behaviors 
		#   (lunge feed, resting, surface active)
		new_beh <- revalue(temp2$whale_behavior, c("BL-Blowing"=1, "DF-Dive-fluke-up"=1, "DN-Dive-no-fluke"=1, 
													"LF-Lunge-feed"=2, "RE-Resting"=2, "SA-Surface-active"=2))
		dat <- cbind(temp2,new_beh)
				
		#  Create "transit" and "stationary" data sets.
		transit <- filter(dat, new_beh == 1)
		station <- filter(dat, new_beh == 2)
 
		#  Create data sets for each behavior.
		blow <- filter(dat, whale_behavior == "BL-Blowing")

		target <- c("DF-Dive-fluke-up", "DN-Dive-no-fluke")
		dive <-filter(dat, whale_behavior %in% target)

		target <- c("SA-Surface-active", "LF-Lunge-feed")
		surf <- filter(dat, whale_behavior %in% target)

		rest <-  filter(dat, whale_behavior == "RE-Resting")
		
		#############################################################################
		########## CHOOSE one of below options for desired subset of data #################
		
		pts1 <- dat
		pts1 <- transit
		pts1 <- station
		
		# Instead, using groups of observationally assigned more specific behavior categories
		# pts1 <- blow
		pts1 <- dive
		pts1 <- surf
		pts1 <- rest
		
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
		pts2 <- arrange(pts1, SwB_Wpt_ID, ObOrder_Time)
		pts3 <- pts2 %>%
					   filter( whale_dist_to_shore_m > 3) %>%
					   filter(SwB_Wpt_ID != "2008-05-28-N-032") %>% 
						filter(SwB_Wpt_ID != "2008-09-15-K-011" ) %>%
						filter(SwB_Wpt_ID != "2010-05-21-N-039") %>%
						filter(SwB_Wpt_ID != "2010-06-25-N-003") %>%
						filter(SwB_Wpt_ID != "2013-07-14-S-002")					
		pts4 <- pts3 %>% 							
					   dplyr::select(SwB_Wpt_ID, whale_easting, whale_northing) %>%
					   dplyr::rename( x = whale_easting , y = whale_northing) 
									
		loc.data <- as.data.frame(pts4)
		loc.data.transit <- as.data.frame(pts4)
		loc.data.station <- as.data.frame(pts4)
		loc.data.dive <- as.data.frame(pts4)
		loc.data.surf <- as.data.frame(pts4)
		
		#  Location data processing: go through data set by "SwB_Wpt_ID", select on whales that have more than 
		#   2 sightings (need at least 3 sightings for step length and turning angle).		
		#  Run move.fun!!!
		data.out <- loc.data %>%
								group_by(SwB_Wpt_ID)%>%
								filter(n()>2) %>%
								do(move.fun(data.frame(x = .$x, y = .$y))) %>%
								as.data.frame(.)
									
		data.out.transit <- loc.data.transit %>%
												group_by(SwB_Wpt_ID)%>%
												filter(n()>2) %>%
												do(move.fun(data.frame(x = .$x, y = .$y))) %>%
												as.data.frame(.)

		data.out.station <- loc.data.station %>%
												group_by(SwB_Wpt_ID)%>%
												filter(n()>2) %>%
												do(move.fun(data.frame(x = .$x, y = .$y))) %>%
												as.data.frame(.)
												
		data.out.dive <- loc.data.dive %>%
												group_by(SwB_Wpt_ID)%>%
												filter(n()>2) %>%
												do(move.fun(data.frame(x = .$x, y = .$y))) %>%
												as.data.frame(.)
												
		data.out.surf <- loc.data.surf %>%
												group_by(SwB_Wpt_ID)%>%
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
		all.steps <- as.vector(data.out$dist)
		all.steps <- na.omit(all.steps)
	    all.steps <- as.numeric(all.steps)
		
		all.steps.transit <- as.vector(data.out.transit$dist)
		all.steps.transit <- na.omit(all.steps.transit)
	    all.steps.transit <- as.numeric(all.steps.transit)
		
		all.steps.station <- as.vector(data.out.station$dist)
		all.steps.station <- na.omit(all.steps.station)
	    all.steps.station <- as.numeric(all.steps.station)
		
		all.steps.dive <- as.vector(data.out.dive$dist)
		all.steps.dive <- na.omit(all.steps.dive)
	    all.steps.dive <- as.numeric(all.steps.dive)
		
		all.steps.surf <- as.vector(data.out.surf$dist)
		all.steps.surf <- na.omit(all.steps.surf)
	    all.steps.surf <- as.numeric(all.steps.surf)
		
		#   Turn angles
		all.turns <- as.vector(data.out$turn.angle)
		all.turns <- na.omit(all.turns)
		all.turns <- as.numeric(all.turns)
		
		all.turns.transit <- as.vector(data.out.transit$turn.angle)
		all.turns.transit <- na.omit(all.turns.transit)
		all.turns.transit <- as.numeric(all.turns.transit)
		
		all.turns.station <- as.vector(data.out.station$turn.angle)
		all.turns.station <- na.omit(all.turns.station)
		all.turns.station <- as.numeric(all.turns.station)
		
		all.turns.dive <- as.vector(data.out.dive$turn.angle)
		all.turns.dive <- na.omit(all.turns.dive)
		all.turns.dive <- as.numeric(all.turns.dive)
		
		all.turns.surf <- as.vector(data.out.surf$turn.angle)
		all.turns.surf <- na.omit(all.turns.surf)
		all.turns.surf <- as.numeric(all.turns.surf)
		
		length(all.steps.transit)
		range(all.steps.transit)
		length(all.turns.transit)
		range(all.turns.transit)
		
		length(all.steps.station)
		range(all.steps.station)
		length(all.turns.station)
		range(all.turns.station)
#############################################################################	
#############################################################################	

####  Looking for erroneous points - crazy step lenghts. Based on detection probability, 
####   removing step lengths greater than 5000m.											
# filter(data.out, dist > 5000)
####  Output: 
		#   SwB_Wpt_ID  	     		 x       			    y 					turn.angle    	heading      		dist
		#  2008-05-28-N-032 		408077.7 		6525132 		-2.8700604  	1.9158865 		18365.832
		#  2008-09-15-K-011 		438057.9 		6488211 		-3.0638468 		0.3093324  		5893.211
		#  2010-05-21-N-039 	    383546.1 		6531663        NA 					-1.2978190 		30239.332
		#  2010-05-21-N-039 		412428.1 		6523589  		3.1412104  		1.8433914 		29989.382
		#  2010-06-25-N-003 		439855.7 		6481496  		0.7511768  		2.9004882  		7495.769
		#  2010-06-25-N-003 		438163.6 		6488788 		-3.1284931 		-0.2280049  	7485.660
		#  2013-07-14-S-002 		442848.6 		6467393        NA 					-1.2835693  	5678.349
####  Remove entire tracks from these points. Back up in pts2 object, 
####   then re run data.out.						
#############################################################################	
		
	# #  Run move.fun for sightings that are within 2000m
		# pts5 <- pts3 %>% 
					   # filter(ship_bb_to_whale_dist < 2001) %>%
					   # dplyr::select(SwB_Wpt_ID, whale_easting, whale_northing) %>%
					   # dplyr::rename( x = whale_easting , y = whale_northing) 
									
		# loc.data.close <- as.data.frame(pts5)
		
		# data.out.close <- loc.data.close %>%
										# group_by(SwB_Wpt_ID)%>%
										# filter(n()>2) %>%
										# do(move.fun(data.frame(x = .$x, y = .$y))) %>%
										# as.data.frame(.)
	
		# #  Run move.fun for sightings that are outside 2000m
		# pts6 <- pts3 %>% 
					   # filter(ship_bb_to_whale_dist > 2000) %>%
					   # dplyr::select(SwB_Wpt_ID, whale_easting, whale_northing) %>%
					   # dplyr::rename( x = whale_easting , y = whale_northing) 
									
		# loc.data.far <- as.data.frame(pts6)
		
		# data.out.far <-  loc.data.far %>%
										# group_by(SwB_Wpt_ID)%>%
										# filter(n()>2) %>%
										# do(move.fun(data.frame(x = .$x, y = .$y))) %>%
										# as.data.frame(.)
	
		
	# close.steps <- ggplot(data.out.close, aes(x=dist)) + 
								# geom_histogram(color="black", fill="grey", binwidth=50) +
								# theme_bw()
	# close.steps
	
	# far.steps <- ggplot(data.out.far, aes(x=dist)) + 
								# geom_histogram(color="black", fill="grey", binwidth=50) +
								# theme_bw()
	# far.steps
		
	# close.turns <- ggplot(data.out.close, aes(x=turn.angle)) + 
								# geom_histogram(color="black", fill="grey", binwidth=0.25) +
								# theme_bw()
	# close.turns
	
	# far.turns<- ggplot(data.out.far, aes(x=turn.angle)) + 
								# geom_histogram(color="black", fill="grey", binwidth=0.25) +
								# theme_bw()
	# far.turns
	
	
	# close.steps <- as.vector(data.out.close$dist)
	# close.steps <- na.omit(close.steps)
	
	# far.steps <- as.vector(data.out.far$dist)
	# far.steps <- na.omit(far.steps)
	
	# close.turns <- as.vector(data.out.close$turn.angle)
	# close.turns <- na.omit(close.turns)
	
	# far.turns <- as.vector(data.out.far$turn.angle)
	# far.turns <- na.omit(far.turns)
#############################################################################		

		#  Manual trouble shooting of code to pull in data and run function										
		# tmp_d <- loc.data %>% group_by(SwB_Wpt_ID) %>% filter(n() > 2)
		# d_in <- split(as.data.frame(tmp_d), tmp_d$SwB_Wpt_ID, drop = T)
		# out_d <- bind_rows(lapply(1:length(d_in), function(i){
			# print(i)
			# move.fun(d_in[[i]])
		# }))	
#############################################################################		

#  Create dummy "dates" vector: this is the other input for the move_param_fun 
		#   Do not need to make dates - can run through move_param_fun without it.
		# dates <- (seq(as.Date("2011/5/22"), as.Date("1995/1/1"), by=-1))
				# If you want to create a column in the locs dataframe for the "dummy" dates column:
				# locs3 <- mutate(locs2, dates = (seq(as.Date("2011/5/22"), as.Date("1995/1/1"), by=-1)))