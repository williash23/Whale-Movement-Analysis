# Sara Williams
# 7/14/2015
# Analysis script to run functions to run movement oarameter generation functions
################################################################################
		#  Load packages
		library(dplyr)

		#  Read in whale observation location data
		locs1<- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/Whale_Pts_SST.csv")
		
		#  Source Movement_parameters_fun
		source(file.path("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/Movement_parameters_fun.R"))
		
		#  Prep locs_data: create "x" and "y" columns
		#   Arrange by same whale ID and observation order.
		locs2 <- arrange(locs1, SwB_Wpt_ID, ObOrder_Time)
		
			#  Create lat/lon XY data from UTM XY data  
			proj4 <- CRS("+proj=utm +zone=8 + datum=WGS84 + ellps=WGS84") 
			tmp <- locs2 %>%
							dplyr::select(x = whale_easting, y = whale_northing) %>%
							as.matrix(.)
			xy_tmp <- SpatialPoints(as.data.frame(tmp), proj4) 
			xy <- as.data.frame(spTransform(xy_tmp, CRS("+proj=longlat + datum=WGS84 + ellps=WGS84")))
			locs3 <- cbind(locs2, xy)
			loc_data <- locs3 %>%
									select(SwB_Wpt_ID, x, y)
		
		#  Prep "dates" dataframe: create "dummy" dates column
		dates <- (seq(as.Date("2011/5/22"), as.Date("1995/1/1"), by=-1))
				# If you want to create a column in the locs dataframe for the "dummy" dates column:
				# locs3 <- mutate(locs2, dates = (seq(as.Date("2011/5/22"), as.Date("1995/1/1"), by=-1)))
				
				
		loc_data %>%
				
		# process_wrap <- function(x){			
							# move_fun(x)
							# p <- proj_fun(x$sst_fn)
							# cr <- crop_fun(x$sst_fn)
						# }			
		
		#  Location data processing: go through data set by "SwB_Wpt_ID" and do movement parameter creation function.
			loc_data %>% 
				group_by(SwB_Wpt_ID) %>%
				do(move_fun(.)) 
				