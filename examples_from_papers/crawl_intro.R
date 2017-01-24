# Sara Williams
# 12/8/2015; updated 2/1/2016, 3/9/2016
# Data prep for model running script.
################################################################################

#  Load packages
library(dplyr)
library(adehabitatLT)    
library(tidyr)
library(stringr)
library(aspace)
library(sp)
library(rgdal)

################################################################################

#  Load data
dat <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/Whales_0615_locations_clean.csv")

#  Generate step lengths and turning angles using ADEpackage
locs_a <- arrange(dat, same_whale_ID, ob_order_time)
locs_b <- locs_a %>%
                dplyr::select(same_whale_ID, unique_event_ID, X_whale_UTM, Y_whale_UTM) %>%
                dplyr::rename(X = X_whale_UTM, Y = Y_whale_UTM)
locs_b$same_whale_ID <- droplevels(locs_b$same_whale_ID)

time_1 <- dat %>%
                 dplyr::select(TimeTxt, unique_event_ID) %>%
                 as.data.frame()
time_2 <- as.data.frame(str_split_fixed(time_1$TimeTxt, ":", 3))
time_3 <- cbind(time_2, time_1$unique_event_ID)
names(time_3)[4] <- "ID"
time_3$V1 <- as.character(time_3$V1)
time_3$V2 <- as.character(time_3$V2)
time_3$V3 <- as.character(time_3$V3)
time_3$V1 <- as.integer(time_3$V1)
time_3$V2 <- as.integer(time_3$V2)
time_3$V3 <- as.integer(time_3$V3)
time_3$V1 <- as.numeric(time_3$V1)
time_3$V2 <- as.numeric(time_3$V2)
time_3$V3 <- as.numeric(time_3$V3)

time_4 <- time_3 %>%
                 rename (hr = V1, min = V2, sec = V3) %>%
                 mutate(sec_frac = sec/60, min_sec = sec_frac+min) %>%
                 mutate(min_frac = min_sec/60, time = hr+min_frac)

time_locs_all <- cbind(locs_b, time_4$time)
names(time_locs_all)[5] <- "Time"

locs_1 <- filter(time_locs_all, same_whale_ID == "2008-05-06-N-005")
final_dat <- dplyr::select(Timelocs_1, Time, X, Y)

coordinates(final_dat) <- ~X+Y
UTM <- CRS("+proj=utm +zone=8 +datum=WGS84")
proj4string(final_dat) <- UTM


initial = list(a=c(coordinates(final_dat)[1,1],0,
                   coordinates(final_dat)[1,2],0),
                   P=diag(c(10000^2,54000^2,10000^2,5400^2)))