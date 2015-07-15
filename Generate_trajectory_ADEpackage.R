#  Sara Williams
#  7/14/2015
#  Script to get trajectories from repeated observations of same whale and to summarize
#  step lengths and turning angles 
################################################################################

#  Load packages
require(stats)
require(plyr)
require(dplyr)
library(adehabitatLT)						
library(adehabitatMA)

#  Read in whale location data (all whale observations)
locs1 <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/Whale_Pts.csv")
head(locs1)

#  Arrange locations in order of the same whale (so in a group of multiple observations)
#   and then arrange also by observation order
locs2 <- arrange(locs1, SwB_Wpt_ID, ObOrder_Time)
locs3 <- mutate(locs2, dummy_date = (seq(as.Date("2011/5/22"), as.Date("1995/1/1"), by=-1)))

					proj4 <- CRS("+proj=utm +zone=8 + datum=WGS84 + ellps=WGS84") 
					
					tmp <- locs3 %>%
								dplyr::select(x = whale_easting, y = whale_northing) %>%
								as.matrix(.)
					xy_tmp <- SpatialPoints(as.data.frame(tmp), proj4) 
					xy <- as.data.frame(spTransform(xy_tmp, CRS("+proj=longlat + datum=WGS84 + ellps=WGS84")))
locs4 <- cbind(locs3, xy)
					
whales_traj <- as.ltraj(xy =locs4[,c("x","y")], date =as.POSIXct(locs4$dummy_date), id = locs4$SwB_Wpt_ID)


whales_traj
head(whales_traj[[1]])
plot(whales_traj[[1]])

#################################################################################################
final2 <- x %>%
	group_by(SwB_Wpt_ID) %>%
	summarise(sum(ObOrder2 > 3))

colnames(final2) <- c("id", "SumObOrder2")
	
locs2 <- final2[ which(final2$SumObOrder2 >0),] 	
locs2

##################################################################################################
######## Example code from adehabitatLT vignette #######################################################
data(puechabonsp)
head(puechabonsp)
locs <- puechabonsp$relocs
locs <- as.data.frame(locs)
head(locs)

da <- as.character(locs$Date)
head(da)
da <- as.POSIXct(strptime(as.character(locs$Date),"%y%m%d"))

puech <- as.ltraj(xy = locs[,c("X","Y")], date = da, id = locs$Name)
puech
head(puech[[1]])
plot(puech)
##################################################################################################