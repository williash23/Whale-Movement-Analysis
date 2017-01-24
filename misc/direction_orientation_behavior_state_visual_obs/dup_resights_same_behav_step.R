#  Sara Williams
#  10/7/2015
#  Script to seperate resights of the same whale for behaviors that change within the resighting 
#  burst. Duplicate a resight so that it matches a behavior of the sighting immediately previous.
################################################################################

#  Load packages
library(plyr)
library(dplyr)

#  Read in whale location data (all whale observations)
temp <- read.csv("C:/Users/sara.williams/Documents/GitHub/General-Data-Cleaning/Whales_0615_general_clean.csv")

#  Subset to only observations where behavior is recorded.
temp2 <-filter(temp, temp$whale_behavior=="BL-Blowing" | temp$whale_behavior=="DF-Dive-fluke-up" | 
							 temp$whale_behavior=="DN-Dive-no-fluke" | temp$whale_behavior=="LF-Lunge-feed" | 
							 temp$whale_behavior=="RE-Resting" | temp$whale_behavior=="SA-Surface-active")
temp2$whale_behavior <- droplevels(temp2$whale_behavior)

#  Make a new whale_behavior category that is numeric and attach to dataframe.
val_beh <- revalue(temp2$whale_behavior, c("BL-Blowing"=1, "DF-Dive-fluke-up"=2, "DN-Dive-no-fluke"=3, 
									   "LF-Lunge-feed"=4, "RE-Resting"=5, "SA-Surface-active"=6))
dat2 <- cbind(temp2, val_beh)

#  Remove erroneous (unreasonable) locations.
dat3 <- dat2 %>%
			 filter(whale_dist_shore_m > 0) %>%
			 filter(whale_depth_m < 0) %>%
			 arrange(same_whale_ID, ob_order_time)

dat_beh <- dat3 %>%
					   dplyr::group_by(same_whale_ID) %>%
					   mutate(num_beh = n_distinct(val_beh)) %>%
					   mutate(num_obs = n()) %>%
					   mutate(duplicate = "0") %>%
					   as.data.frame(.)

# # # write.csv(dat_beh, "Whales_0615_locations_clean.csv")					   
################################################################################

#  Adjust data set: break each group of same whale sightings into sub-groups based on changes
#   to observed behavior.  For example:

# # # EvB_Wpt_Id       			same_whale_ID  			whale_behavior 			val_beh 
# # # 2008-05-25-N-005 		2008-05-25-N-005        RE-Resting       				2
# # # 2008-05-25-N-006 		2008-05-25-N-005		DF-Dive-fluke-up     	1
# # # BECOMES......
# # # EvB_Wpt_Id       			same_whale_ID 			whale_behavior 			val_beh 		new_whale_beh			new_val_beh
# # # 2008-05-25-N-005 		2008-05-25-N-005a      RE-Resting       				2						RE-Resting       				2	
# # # 2008-05-25-N-006 		2008-05-25-N-005a		DF-Dive-fluke-up     	1						RE-Resting       				2
# # # 2008-05-25-N-006 		2008-05-25-N-005b		DF-Dive-fluke-up     	1						DF-Dive-fluke-up     	1

# # # EvB_Wpt_Id       			same_whale_ID			whale_behavior 			val_beh
# # # 2008-05-26-N-034		2008-05-26-N-034 		SA-Surface-active        2
# # # 2008-05-26-N-035 		2008-05-26-N-034  		DF-Dive-fluke-up         1
# # # 2008-05-26-N-036 		2008-05-26-N-034  		DF-Dive-fluke-up         1
# # # BECOMES......
# # # EvB_Wpt_Id       			same_whale_ID			whale_behavior 			val_beh 		new_whale_beh			new_val_beh
# # # 2008-05-26-N-034		2008-05-26-N-034a 		SA-Surface-active        2						SA-Surface-active        2
# # # 2008-05-26-N-035 		2008-05-26-N-034a 		DF-Dive-fluke-up         1						SA-Surface-active        2
# # # 2008-05-26-N-035 		2008-05-26-N-034b 	DF-Dive-fluke-up         1						DF-Dive-fluke-up         1
# # # 2008-05-26-N-036 		2008-05-26-N-034b  	DF-Dive-fluke-up         1						DF-Dive-fluke-up         1
################################################################################

beh_mult <- dat_beh %>%
						  dplyr::filter(num_beh > 1) %>%
						  group_by(same_whale_ID) %>%
						  dplyr::filter(row_number() > 1, val_beh != lag(val_beh)) %>%
						  mutate(duplicate = "1") %>%
						  ungroup(.) %>%
						  as.data.frame(.)

dup_dat <- as.data.frame(bind_rows(dat_beh, beh_mult))
dup_dat$duplicate <- as.numeric(dup_dat$duplicate)
dup_dat <- dup_dat %>%
					   arrange(same_whale_ID, ob_order_time, desc(duplicate)) %>%
					   group_by(same_whale_ID) %>%
					   mutate(new_beh =  ifelse(duplicate == 1, lag(val_beh), val_beh)) %>%
					   mutate(n_row = n()) %>%
					   as.data.frame(.)	

# # # write.csv(dup_dat, "Whales_0615_w_dup.csv")		
					   
rem_last_fun <- function(x){
									if(x$duplicate[nrow(x) - 1] == 1 ){
										out <- x %>% slice(1:(n()-1))
										}else{
										out <- x
										}
									}

#  Remove unccesarry duplicated last row and observations of only one row
dup_rem_last <- dup_dat %>%
								  group_by(same_whale_ID) %>%
								  filter(n() > 1) %>%
								  do(rem_last_fun(.))%>%
								  as.data.frame(.)

# # # write.csv(dup_rem_last, "Whales_0615_w_dup_last_row_rem.csv")							   
################################################################################

#  Duplicate every sighting except the first and last?

