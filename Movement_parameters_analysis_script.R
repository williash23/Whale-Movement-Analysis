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
		#   Arrange by same whale ID and observation order.
		pts2 <- arrange(pts1, SwB_Wpt_ID, ObOrder_Time)
		pts3 <- pts2 %>%
									dplyr::select(SwB_Wpt_ID, whale_easting, whale_northing) %>%
									dplyr::rename( x = whale_easting , y = whale_northing) 
									
		loc.data <- as.data.frame(pts3)
									
		#  Create dummy "dates" vector: this is the other input for the move_param_fun 
		#   Do not need to make dates - can run through move_param_fun without it.
		# dates <- (seq(as.Date("2011/5/22"), as.Date("1995/1/1"), by=-1))
				# If you want to create a column in the locs dataframe for the "dummy" dates column:
				# locs3 <- mutate(locs2, dates = (seq(as.Date("2011/5/22"), as.Date("1995/1/1"), by=-1)))
		
	#  Location data processing: go through data set by "SwB_Wpt_ID" and do movement parameter creation function.		
		data.out <- loc.data %>%
						group_by(SwB_Wpt_ID)%>%
						filter(n()>2) %>%
						do(move.fun(data.frame(x = .$x, y = .$y))) %>%
						as.data.frame(.)
																	
		#  Manual trouble shooting of code to pull in data and run function										
		# tmp_d <- loc.data %>% group_by(SwB_Wpt_ID) %>% filter(n() > 2)
		# d_in <- split(as.data.frame(tmp_d), tmp_d$SwB_Wpt_ID, drop = T)
		# out_d <- bind_rows(lapply(1:length(d_in), function(i){
			# print(i)
			# move.fun(d_in[[i]])
		# }))	
		
				data.out <- loc.data %>%
						group_by(SwB_Wpt_ID)%>%
						filter(n()>2) %>%
						do(move.fun(data.frame(x = .$x, y = .$y))) %>%
						as.data.frame(.)

		
	
	