#  Sara Williams
#  10/7/2015
#  Script to seperate resights of the same whale for behaviors that change within the resighting 
#  burst. Duplicate a resight so that it matches a behavior of the sighting immediately previous.
################################################################################

#  Load packages
library(plyr)
library(dplyr)

#  Read in whale location data (all whale observations)
temp <- dat4
#read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/whale_points_id_errors_rem.csv")
temp2 <-filter(temp, temp$whale_behavior=="BL-Blowing" | temp$whale_behavior=="DF-Dive-fluke-up" | 
							 temp$whale_behavior=="DN-Dive-no-fluke" | temp$whale_behavior=="LF-Lunge-feed" | 
							 temp$whale_behavior=="RE-Resting" | temp$whale_behavior=="SA-Surface-active")
temp2$whale_behavior <- droplevels(temp2$whale_behavior)
temp2 <- temp2 %>%
				   group_by(same_whale_ID)%>%
				   as.data.frame(.)
				   
#  Make a new whale_behavior category that is numeric and attach to dataframe.
val_beh <- revalue(temp2$whale_behavior, c("BL-Blowing"=1, "DF-Dive-fluke-up"=2, "DN-Dive-no-fluke"=3, 
									   "LF-Lunge-feed"=4, "RE-Resting"=5, "SA-Surface-active"=6))
dat <- cbind(temp2, val_beh)

#  From preliminary trajectory from relocations using ADEhabitatLT package, a few step lengths are
#   unreasonable. These are:
# # # tmp <- filter(traj_dat, dist > 10000)
# # # x       				y 					date        dx        				dy     					dist 				 dt       R2n 						 id 			SwB_Wpt_ID
# # # 390794.6 		6531345   	26  			17283.07 			-6212.822 		18365.83 		 1    	581461.2 				10    		2008-05-28-N-032
# # # 412665.7 		6523511 		1382 		-29119.64  		8152.515 			30239.33  	 1         0.0  						502  		2010-05-21-N-039
# # # 383546.1 		6531663 		1383  	28882.04 			-8074.089 		29989.38  	 1 		914417196.6		502   		2010-05-21-N-039

#  Remove erroneous (unreasonable) locations.
dat <- dat %>%
			 filter(whale_dist_to_shore_m > 0) %>%
			 filter(whale_depth_m < 0) %>%
			 filter(SwB_Wpt_ID != "2008-05-28-N-032") %>% 
			 filter(SwB_Wpt_ID != "2010-05-21-N-039") %>%
			 arrange(SwB_Wpt_ID, ObOrder_Time)

dat_beh <- dat %>%
					   dplyr::group_by(SwB_Wpt_ID) %>%
					   mutate(num_beh = n_distinct(val_beh)) %>%
					   mutate(num_obs = max(ObOrder_Time)) %>%
					   mutate(duplicate = "0") %>%
					   as.data.frame(.)

# # # write.csv(dat_beh, "all_locs_not dup.csv")					   
################################################################################

#  Adjust data set: break each group of same whale sightings into sub-groups based on changes
#   to observed behavior.  For example:

# # # EvB_Wpt_Id       			SwB_Wpt_ID    				whale_behavior 			val_beh 
# # # 2008-05-25-N-005 		2008-05-25-N-005        RE-Resting       				2
# # # 2008-05-25-N-006 		2008-05-25-N-005		DF-Dive-fluke-up     	1
# # # BECOMES......
# # # EvB_Wpt_Id       			SwB_Wpt_ID    				whale_behavior 			val_beh 		new_whale_beh			new_val_beh
# # # 2008-05-25-N-005 		2008-05-25-N-005a      RE-Resting       				2						RE-Resting       				2	
# # # 2008-05-25-N-006 		2008-05-25-N-005a		DF-Dive-fluke-up     	1						RE-Resting       				2
# # # 2008-05-25-N-006 		2008-05-25-N-005b		DF-Dive-fluke-up     	1						DF-Dive-fluke-up     	1

# # # EvB_Wpt_Id       			SwB_Wpt_ID    				whale_behavior 			val_beh
# # # 2008-05-26-N-034		2008-05-26-N-034 		SA-Surface-active        2
# # # 2008-05-26-N-035 		2008-05-26-N-034  		DF-Dive-fluke-up         1
# # # 2008-05-26-N-036 		2008-05-26-N-034  		DF-Dive-fluke-up         1
# # # BECOMES......
# # # EvB_Wpt_Id       			SwB_Wpt_ID    				whale_behavior 			val_beh 		new_whale_beh			new_val_beh
# # # 2008-05-26-N-034		2008-05-26-N-034a 		SA-Surface-active        2						SA-Surface-active        2
# # # 2008-05-26-N-035 		2008-05-26-N-034a 		DF-Dive-fluke-up         1						SA-Surface-active        2
# # # 2008-05-26-N-035 		2008-05-26-N-034b 	DF-Dive-fluke-up         1						DF-Dive-fluke-up         1
# # # 2008-05-26-N-036 		2008-05-26-N-034b  	DF-Dive-fluke-up         1						DF-Dive-fluke-up         1
################################################################################

beh_mult <- dat_beh %>%
						  filter(n() > 1) %>%
						  filter(num_beh > 1) %>%
						  group_by(SwB_Wpt_ID) %>%
						  filter(row_number() > 1, val_beh != lag(val_beh)) %>%
						  mutate(duplicate = "1") %>%
						  ungroup(.) %>%
						  as.data.frame(.)

dup_dat <- as.data.frame(bind_rows(dat_beh, beh_mult))
dup_dat$duplicate <- as.numeric(dup_dat$duplicate)
dup_dat <- dup_dat %>%
					   arrange(SwB_Wpt_ID, ObOrder_Time, desc(duplicate)) %>%
					   group_by(SwB_Wpt_ID) %>%
					   mutate(new_beh =  ifelse(duplicate == 1, lag(val_beh), val_beh)) %>%
					   mutate(n_row = n()) %>%
					   as.data.frame(.)	

# # # write.csv(dup_dat, "all_locs_w_dup.csv")		
					   
rem_last_fun <- function(x){
									if(x$duplicate[nrow(x) - 1] == 1 ){
										out <- x %>% slice(1:(n()-1))
										}else{
										out <- x
										}
									}

#  Remove unccesarry duplicated last row and observations of only one row
dup_rem_last <- dup_dat %>%
								  group_by(SwB_Wpt_ID) %>%
								  filter(n() > 1) %>%
								  do(rem_last_fun(.))%>%
								  as.data.frame(.)

# # # write.csv(dup_rem_last, "all_locs_w_dup_last_row_rem.csv")							   

