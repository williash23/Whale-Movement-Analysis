
#  Load packages
library(dplyr)
library(adehabitatLT)    
library(tidyr)
library(stringr)
library(aspace)
################################################################################

#  Load data
dat_raw <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/Whales_0615_locations_clean.csv")

#  Data clean up and manipulation
tmp1 <- dat_raw %>%
            filter(year > 2009) %>% ## per info from Karin, PTB (bearing) was kind of weird before 2010
            filter(ship_whale_dist < 8000) %>%
            filter(count == 1) %>%            
            group_by(same_whale_ID) %>%
            filter(n() > 1) %>%
            ungroup() %>%
            as.data.frame()
tmp2 <- arrange(tmp1, same_whale_ID, ob_order_time)
tmp3 <- tmp2 %>%
              dplyr::select(same_whale_ID, X_whale_UTM, Y_whale_UTM, ship_whale_dist, 
                                     ship_whale_bearing, whale_behavior, ob_order_time, TimeTxt, same_whale_ID) %>%
              dplyr::rename(X = X_whale_UTM, Y = Y_whale_UTM)
num_cues <- tmp3 %>%
                    group_by(same_whale_ID) %>%
                    mutate(num_cues = n()) %>%
                    ungroup() %>%
                    as.data.frame()

#  Add time difference between successive observations to observations.
t1 <- tmp3
t2 <- as.data.frame(str_split_fixed(t1$TimeTxt, ":", 3))
t3 <- cbind(t2, t1$X, t1$Y, t1$ob_order_time, t1$ship_whale_dist, 
                    t1$ship_whale_bearing, t1$whale_behavior, t1$TimeTxt, t1$same_whale_ID)
names(t3)[4] <- "X"
names(t3)[5] <- "Y"
names(t3)[6] <- "ob_order_time"
names(t3)[7] <- "ship_whale_dist"
names(t3)[8] <- "ship_whale_bearing"
names(t3)[9] <- "whale_behavior"
names(t3)[10] <- "time_txt"
names(t3)[11] <- "same_whale_ID"
t3$V1 <- as.numeric(as.integer(as.character(t3$V1)))
t3$V2 <- as.numeric(as.integer(as.character(t3$V2)))
t3$V3 <- as.numeric(as.integer(as.character(t3$V3)))
#   time_diff is the difference in time between a sighting and the sighting before it
t4 <- t3 %>%
         dplyr::rename (hr = V1, min = V2, sec = V3) %>%
         dplyr::mutate(sec_frac = sec/60, min_sec = sec_frac+min) %>%
         dplyr::mutate(min_frac = min_sec/60, time_dec = hr+min_frac) %>%
         group_by(same_whale_ID) %>%
         mutate(time_diff = lead(time_dec) - time_dec) %>%
         mutate(time_diff_sec = time_diff*3600) %>%
         ungroup() %>%
         as.data.frame()

#  Create trajectory (ltraj) object for steps and turns
tmp4 <- t4 %>%
               dplyr::select(X, Y, ob_order_time, time_txt, time_dec, ship_whale_dist,
                                     ship_whale_bearing, whale_behavior, same_whale_ID, time_diff_sec)






#  Whale activity
#   Deep dive > 120s, surface interval dive < 50
surf <- filter(tmp4, time_diff_sec < 50)
surf$same_whale_ID <- droplevels(surf$same_whale_ID)

traj_surf <- as.ltraj(xy = surf[,c("X","Y")], id = surf$same_whale_ID, typeII = FALSE)
#  Convert into dataframe
traj_surf_df_tmp <- ld(traj_surf)
traj_surf_df <- traj_surf_df_tmp %>%
                 dplyr::select(id, x, y, dist, rel.angle, abs.angle, dx, dy) 
names(traj_surf_df) <- c("same_whale_ID", "X", "Y", "step", "turn", "abs_angle", "dx", "dy")
surf_dat <- filter(traj_surf_df, step < 1943)

deep <- filter(tmp4, time_diff_sec > 120)
deep$same_whale_ID <- droplevels(deep$same_whale_ID)

traj_deep <- as.ltraj(xy = deep[,c("X","Y")], id = deep$same_whale_ID, typeII = FALSE)
#  Convert into dataframe
traj_deep_df_tmp <- ld(traj_deep)
traj_deep_df <- traj_deep_df_tmp %>%
                 dplyr::select(id, x, y, dist, rel.angle, abs.angle, dx, dy) 
names(traj_deep_df) <- c("same_whale_ID", "X", "Y", "step", "turn", "abs_angle", "dx", "dy")
deep_dat <- filter(traj_deep_df, step < 1943)





#  Distance to ship at first sighting
#   Near (<1000m), far away (> 3000m)
near_ids <- tmp4 %>%
                    group_by(same_whale_ID) %>%
                    slice(1) %>%
                    filter(ship_whale_dist < 1000) %>%
                    ungroup() %>%
                    as.data.frame(.)
near <- semi_join(tmp4, near_ids, by = c("same_whale_ID"))
near$same_whale_ID <- droplevels(near$same_whale_ID)
#near <- semi_join(tmp4, near_ids, by = c("same_whale_ID", "X", "Y")  ?? correct join?

traj_near <- as.ltraj(xy = near[,c("X","Y")], id = near$same_whale_ID, typeII = FALSE)
#  Convert into dataframe
traj_near_df_tmp <- ld(traj_near)
traj_near_df <- traj_near_df_tmp %>%
                 dplyr::select(id, x, y, dist, rel.angle, abs.angle, dx, dy) 
names(traj_near_df) <- c("same_whale_ID", "X", "Y", "step", "turn", "abs_angle", "dx", "dy")
near_dat <- filter(traj_near_df, step < 1943)


far_ids <- tmp4 %>%
                  group_by(same_whale_ID) %>%
                  slice(1) %>%
                  filter(ship_whale_dist > 3000) %>%
                  ungroup() %>%
                  as.data.frame(.)
far <- semi_join(tmp4, far_ids, by = "same_whale_ID") 
far$same_whale_ID <- droplevels(far$same_whale_ID)

traj_far <- as.ltraj(xy = far[,c("X","Y")], id = far$same_whale_ID, typeII = FALSE)
#  Convert into dataframe
traj_far_df_tmp <- ld(traj_far)
traj_far_df <- traj_far_df_tmp %>%
                 dplyr::select(id, x, y, dist, rel.angle, abs.angle, dx, dy) 
names(traj_far_df) <- c("same_whale_ID", "X", "Y", "step", "turn", "abs_angle", "dx", "dy")
far_dat <- filter(traj_far_df, step < 1943)





#  Relative bearing (approach) at first sighting
#   Side approach( greater than +/- 40 deg relative to bow), front approach (less than +/- 20)
side_ids <- tmp4 %>%
                    group_by(same_whale_ID) %>%
                    slice(1) %>%
                    filter(ship_whale_bearing < -40 | ship_whale_bearing > 40) %>%
                    ungroup() %>%
                    as.data.frame(.)
side <- semi_join(tmp4, side_ids, by = c("same_whale_ID"))

traj_side <- as.ltraj(xy = side[,c("X","Y")], id = side$same_whale_ID, typeII = FALSE)
#  Convert into dataframe
traj_side_df_tmp <- ld(traj_side)
traj_side_df <- traj_side_df_tmp %>%
                 dplyr::select(id, x, y, dist, rel.angle, abs.angle, dx, dy) 
names(traj_side_df) <- c("same_whale_ID", "X", "Y", "step", "turn", "abs_angle", "dx", "dy")
side_dat <- filter(traj_side_df, step < 1943)


front_ids <- tmp4 %>%
                    group_by(same_whale_ID) %>%
                    slice(1) %>%
                    filter(ship_whale_bearing < 20 | ship_whale_bearing < -20) %>%
                    ungroup() %>%
                    as.data.frame(.)
front <- semi_join(tmp4, front_ids, by = c("same_whale_ID"))

traj_front <- as.ltraj(xy = front[,c("X","Y")], id = front$same_whale_ID, typeII = FALSE)
#  Convert into dataframe
traj_front_df_tmp <- ld(traj_front)
traj_front_df <- traj_front_df_tmp %>%
                 dplyr::select(id, x, y, dist, rel.angle, abs.angle, dx, dy) 
names(traj_front_df) <- c("same_whale_ID", "X", "Y", "step", "turn", "abs_angle", "dx", "dy")
front_dat <- filter(traj_front_df, step < 1943)


#  Check for differences
tst_far_1 <- dplyr::select(far_dat, X, Y, step, turn)
tst_far_2 <- dplyr::select(far, X, Y, step, turn) ## from "far" is from original analysis prep script
which_diff_far1 <-  setdiff(tst_far_1, tst_far_2) #  which_diff are rows that appear in tst_far_1, but not tst_far_2
which_diff_far2 <-  setdiff(tst_far_2, tst_far_1) #  which_diff are rows that appear in tst_far_2, but not tst_far_1


tst_near_1 <- dplyr::select(near_dat, X, Y, step, turn)
tst_near_2 <- dplyr::select(near, X, Y, step, turn) ## from "near" is from original analysis prep script
which_diff_near1 <- setdiff(tst_near_1, tst_near_2) #  which_diff are rows that appear in tst_near_1, but not tst_near_2
which_diff_near2 <- setdiff(tst_near_2, tst_near_1) #  which_diff are rows that appear in tst_near_2, but not tst_near_1
