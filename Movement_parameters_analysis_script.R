# Sara Williams
# 7/14/2015
# Analysis script to run movement parameter generation functions
################################################################################
		#  Load packages
		require(sp)
		require(raster)
		library(rgdal)
		library(dplyr)
		
		#  Read in whale observation location data
		pts1<- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/Whale_Pts_SST.csv")
		
		#  Source Movement_parameters_fun
		source(file.path("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/Move_param_fun.R"))
		
		#  Prep locs.data: this is one of the inputs for the move_param_fun
		#   Arrange by same whale ID and observation order. Also remove pts where distance between whale and shore is less than 0 (erroneous?)
		pts2 <- arrange(pts1, SwB_Wpt_ID, ObOrder_Time)
		pts3 <- pts2 %>%
					   filter( whale_dist_to_shore_m > 0) %>%
					   filter(SwB_Wpt_ID != "2008-05-28-N-032") %>% 
						filter(SwB_Wpt_ID != "2008-09-15-K-011" ) %>%
						filter(SwB_Wpt_ID != "2010-05-21-N-039") %>%
						filter(SwB_Wpt_ID != "2010-06-25-N-003") %>%
						filter(SwB_Wpt_ID != "2013-07-14-S-002")					
		pts4 <- pts3 %>% 							
					   dplyr::select(SwB_Wpt_ID, whale_easting, whale_northing) %>%
					   dplyr::rename( x = whale_easting , y = whale_northing) 
									
		loc.data <- as.data.frame(pts4)
									
		#  Create dummy "dates" vector: this is the other input for the move_param_fun 
		#   Do not need to make dates - can run through move_param_fun without it.
		# dates <- (seq(as.Date("2011/5/22"), as.Date("1995/1/1"), by=-1))
				# If you want to create a column in the locs dataframe for the "dummy" dates column:
				# locs3 <- mutate(locs2, dates = (seq(as.Date("2011/5/22"), as.Date("1995/1/1"), by=-1)))
		
		#  Location data processing: go through data set by "SwB_Wpt_ID", select on whales that have more than 
		#   2 sightings, and do movement parameter creation function.		
		#  Run move.fun!!!
		data.out <- loc.data %>%
								group_by(SwB_Wpt_ID)%>%
								filter(n()>2) %>%
								do(move.fun(data.frame(x = .$x, y = .$y))) %>%
								as.data.frame(.)
		
		#  Run move.fun for sightings that are within 2000m
		pts5 <- pts3 %>% 
					   filter(ship_bb_to_whale_dist < 2001) %>%
					   dplyr::select(SwB_Wpt_ID, whale_easting, whale_northing) %>%
					   dplyr::rename( x = whale_easting , y = whale_northing) 
									
		loc.data.close <- as.data.frame(pts5)
		
		data.out.close <- loc.data.close %>%
										group_by(SwB_Wpt_ID)%>%
										filter(n()>2) %>%
										do(move.fun(data.frame(x = .$x, y = .$y))) %>%
										as.data.frame(.)
	
		#  Run move.fun for sightings that are outside 2000m
		pts6 <- pts3 %>% 
					   filter(ship_bb_to_whale_dist > 2000) %>%
					   dplyr::select(SwB_Wpt_ID, whale_easting, whale_northing) %>%
					   dplyr::rename( x = whale_easting , y = whale_northing) 
									
		loc.data.far <- as.data.frame(pts6)
		
		data.out.far <-  loc.data.far %>%
										group_by(SwB_Wpt_ID)%>%
										filter(n()>2) %>%
										do(move.fun(data.frame(x = .$x, y = .$y))) %>%
										as.data.frame(.)
		
		
		####  Looking for erroneous points - crazy step lenghts. Based on detection probability, removing step lengths greater than 5000m.														
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
	####  Remove entire tracks from these points. Back up in pts2 object, then re run data.out.						
	
	#  Compare histograms of step lengths and turning angles of all sightings vs sightings just within 2000m
	library(ggplot2)
	#   Step lengths
	all.steps <- ggplot(data.out, aes(x=dist)) + 
							geom_histogram(color="black", fill="grey", binwidth=50) +
							theme_bw()
	all.steps					
	
	close.steps <- ggplot(data.out.close, aes(x=dist)) + 
								geom_histogram(color="black", fill="grey", binwidth=50) +
								theme_bw()
	close.steps
	
	far.steps <- ggplot(data.out.far, aes(x=dist)) + 
								geom_histogram(color="black", fill="grey", binwidth=50) +
								theme_bw()
	far.steps
	
	#  Turn angles
	all.turns <- ggplot(data.out, aes(x=turn.angle)) + 
							geom_histogram(color="black", fill="grey", binwidth=0.25) +
							theme_bw()
	all.turns				
	
	close.turns <- ggplot(data.out.close, aes(x=turn.angle)) + 
								geom_histogram(color="black", fill="grey", binwidth=0.25) +
								theme_bw()
	close.turns
	
	far.turns<- ggplot(data.out.far, aes(x=turn.angle)) + 
								geom_histogram(color="black", fill="grey", binwidth=0.25) +
								theme_bw()
	far.turns
	
	#  Create vectors of just turns and step lengths and remove NA's
	#   Step lengths
	all.step_lengths <- as.vector(data.out$dist)
	all.step_lengths <- na.omit(all.step_lengths)
	
	close.step_lengths <- as.vector(data.out.close$dist)
	close.step_lengths <- na.omit(close.step_lengths)
	
	far.step_lengths <- as.vector(data.out.far$dist)
	far.step_lengths <- na.omit(far.step_lengths)
	
	#   Turn angles
	all.turning_angles <- as.vector(data.out$turn.angle)
	all.turning_angles <- na.omit(all.turning_angles)
	
	close.turning_angles <- as.vector(data.out.close$turn.angle)
	close.turning_angles <- na.omit(close.turning_angles)
	
	far.turning_angles <- as.vector(data.out.far$turn.angle)
	far.turning_angles <- na.omit(far.turning_angles)
	
		
		#  Manual trouble shooting of code to pull in data and run function										
		# tmp_d <- loc.data %>% group_by(SwB_Wpt_ID) %>% filter(n() > 2)
		# d_in <- split(as.data.frame(tmp_d), tmp_d$SwB_Wpt_ID, drop = T)
		# out_d <- bind_rows(lapply(1:length(d_in), function(i){
			# print(i)
			# move.fun(d_in[[i]])
		# }))	
			
	