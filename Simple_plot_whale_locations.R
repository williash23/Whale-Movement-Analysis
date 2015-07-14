library(ggmap)
library(sp)
library(rgdal)

########################
#  Get coordinates of GLBA if needed
geocode("Glacier Bay National Park")

#  Read in whale observation data. This data frame currently has whale locations in UTM (zone 8N) format only.
data <- read.csv("C:/Users/sara.williams/Desktop/R_Spatial_Analysis/Whale_Pts.csv")

#  Create bounding box for study area map.
bounds <- c(-137.5, 57.0, -134.5, 59.0)
#  Create background map to hold whale location points.
SA <- get_map(location=bounds, source="stamen", maptype="watercolor", crop=FALSE) 

#  Convert UTM coordinates in data file to XY data for displaying points.
coordinates(data) <- c("whale_easting", "whale_northing")
		proj4string(data) <- CRS("+proj=utm +zone=8 + datum=WGS84 + ellps=WGS84") 
		xy <- as.data.frame(spTransform(data, CRS("+proj=longlat + datum=WGS84 + ellps=WGS84")))

#  Display map and whale locations.
ggmap(SA)+
geom_point(aes(x = x, y = y), data = xy, alpha = .5, color="darkred", size = 1)


#  Display map and whale locations - locations color-coded by behavior type.
ggmap(SA)+
geom_point(aes(x = x, y = y, colour=Behavior, size=Behavior), data = xy, alpha = .5, size = 1)

#  Bin by behavior and show counts of each behavior in each grid cell
ggmap(SA)+
stat_bin2d(aes(x = x, y = y, colour=Behavior, size=Behavior), 
				data = xy, alpha = .5, size = 1, bins=30)

#  Connect points from same whale and add labels and title.
ggmap(SA)+
geom_line(aes(x = x, y = y, group=SwB_Wpt_ID), data = xy, alpha = .5, color="darkred", size = 1)+
scale_size(range=c(3,20))+
labs(x = 'Longitude', y = 'Latitude')+
ggtitle('Whale locations in GLBA and adjacent waters')


#  Subset to only feeding
xy_feed <- subset(xy, Behavior=="LF-Lunge-feed")

#  Display map and whale locations - density plot by year.
ggmap(SA)+
stat_density2d(aes(x = x, y = y, fill= ..level..), size = 2, bins=5, data = xy, geom='polygon')+
scale_fill_gradient('Observation Density')+
facet_wrap(~Behavior)

+
scale_alpha(range = c(.4, .75), guide = FALSE)+
guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))

+
facet_wrap(~ Year)




#  Add in shapefile of ship tracks. 
ships <- readOGR(dsn="C:/Users/sara.williams/Documents/UM/Alaska/AK_2013/2013_Sara/3_FA", layer="Trk_2013S") 
proj4string(ships)
shipsDF <- as.data.frame(ships)

##  Display map and ship locations.
ggmap(SA)+
geom_path(aes(x=Longitude, y=Latitude, group=Event_ID), data=shipsDF, color ="orangered4", alpha = .4, size = .5)




proj4string(shpData) # describes data’s current coordinate reference system 
#  To change to correct projection: 
shpData <- spTransform(shpData, CRS("+proj=longlat +datum=WGS84"))




ggmap(SA)+
geom_line(aes(x = X_PROJ, y = Y_PROJ, group=Event_ID), data = xy, alpha = .5, color="white", size = 0.5)


