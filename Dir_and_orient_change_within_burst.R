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
temp <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/all_locs_not_dup.csv")

#  Only use observations where there is more than one observation for comparison of behaviors
temp2 <- temp %>%
				   group_by(SwB_Wpt_ID)%>%
				   filter(n() > 1) %>% 
				   as.data.frame(.)

# levels(dat$whale_direction)
# ""           "A-Away"     "G-Against"  "G/A-Ag/Aw"  "G/T-Ag/to"  "GA-Ag/Away" "GT-Ag/To"   "T-Toward"   "U-Unsure"   "W-With"     "W/A-Wi/Aw"  "W/T-Wi/T0"  "W/T-Wit/To" "WA-With/Aw" "WT-With/To"
# levels(dat$whale_orientation)
# ""                 "0-Parallel"       "4-45-degrees"     "45-45-degrees"    "9-Perpendicular"  "90-Perpendicular" "U-Unsure"         "U - Unsure" 

direct <- revalue(temp2$whale_direction, c("U-Unsure"=0, "A-Away"=1, "T-Toward"=2, "G-Against"=3, "W-With"=4, "G/A-Ag/Aw"=5,  "GA-Ag/Away"=5, "G/T-Ag/to"=6, "GT-Ag/To"=6,  
			 "W/A-Wi/Aw"=7, "WA-With/Aw"=7, "W/T-Wi/T0"=8, "W/T-Wit/To"=8, "WT-With/To"=8))
orient <-revalue(temp2$whale_orientation, c("U-Unsure"=0,"U - Unsure"=0, "0-Parallel"=1, "4-45-degrees"=3, "45-45-degrees"=3, "9-Perpendicular"=2, "90-Perpendicular"=2))
			 
dat <- cbind(temp2, direct, orient)			 

#  Create variable for CPA, before, after
		split.cpa <- function(x){
			
			nms <- c("before", "after", "closest")
			
			x[x > 0] <- 2
			x[x < 0] <- 1
			x[x == 0] <- 3
			
			out <- nms[x]
		out
		}
		
dat_new <- dat %>%
						group_by(SwB_Wpt_ID) %>%
						arrange(ObOrder_CPA) %>%
						mutate(pos = ObOrder_Time - ObOrder_Time[first(ObOrder_CPA)], ba = split.cpa(pos)) %>%
						as.data.frame(.)
						
just_before_CPA <- dat_new %>%
										 filter(pos == -1 | pos == 0) %>%
										 group_by(SwB_Wpt_ID) %>%
										 arrange(pos) %>%
										 mutate(match_dir_lag = (direct == lag(direct)))%>%
										 mutate(match_ort_lag = (orient == lag(orient)))%>%
										 as.data.frame(.)
just_after_CPA <- dat_new %>%
										 filter(pos == 1 | pos == 0) %>%
										 group_by(SwB_Wpt_ID) %>%
										 arrange(pos) %>%
										 mutate(match_dir_lag = (direct == lag(direct)))%>%
										 mutate(match_ort_lag = (orient == lag(orient)))%>%
										 as.data.frame(.)
										 
#  Same direction within burst?
dat_change <- dat_new %>%
							  group_by(SwB_Wpt_ID)%>%
							  mutate(match_dir_lag = (direct == lag(direct)))%>%
							  mutate(match_ort_lag = (orient == lag(orient)))%>%
							  mutate(match_dir_first = (direct == first(direct)))%>%
							  mutate(match_ort_first = (orient == first(orient)))%>%
							  as.data.frame(.)


# # # c("U-Unsure"="Unsure", "A-Away"="Away", "T-Toward"="Toward", "G-Against"="Against", "W-With"="With", "G/A-Ag/Aw"="Against/Away",  "GA-Ag/Away"="Against/Away",,
# # # "G/T-Ag/to"="Against/Toward", "GT-Ag/To"="Against/Toward", "W/A-Wi/Aw"="With/Away", "WA-With/Aw"="With/Away", "W/T-Wi/T0"="With/Toward", 
# # # "W/T-Wit/To"="With/Toward", "WT-With/To"="With/Toward"))