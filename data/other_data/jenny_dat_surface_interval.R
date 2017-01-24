j_dat_tmp1 <-read.csv("C:/Users/sara.williams/Desktop/j_dat_all.csv")
#j_dat_tmp1 <-read.csv("C:/Users/sara.williams/Downloads/Transcript_08-21-2016_ip.csv")

# j_dat_tmp1$event_id_new <- paste(j_dat_tmp1$vr_source, j_dat_tmp1$event_id, sep = "_")
# j_dat_tmp2 <- j_dat_tmp1 %>%
                       # group_by(event_id_new) %>%
                       # filter(n() > 1) %>%
                       # as.data.frame()
# j_dat <- j_dat_tmp2 %>%
             # group_by(event_id_new, order_cue) %>%
             # slice(1) %>%
             # as.data.frame()


j_dat_tmp2 <- j_dat_tmp1 %>%
                       group_by(event_id) %>%
                       filter(n() > 1) %>%
                       as.data.frame()
j_dat <- j_dat_tmp2 %>%
             group_by(event_id, order_cue) %>%
             slice(1) %>%
             as.data.frame()

fluke_up <- filter(j_dat_tmp2, cue == "fu")
fluke_up_ids <- as.data.frame(distinct(fluke_up))
last_cue_fu <- semi_join(j_dat_tmp1, fluke_up_ids, by = "event_id") %>%
                         group_by(event_id) %>%
                         arrange(order_all, order_cue) %>%
                         filter(row_number()==n()) %>%
                         as.data.frame()


#  Add time difference between successive observations to observations.
time_1 <- j_dat %>%
                 dplyr::select(event_id, cue, order_all, order_cue, sighted) %>%
                 as.data.frame()
time_2 <- as.data.frame(str_split_fixed(time_1$sighted, ":", 3))
time_3 <- cbind(time_2, time_1$event_id, time_1$cue, time_1$order_all, time_1$order_cue, time_1$sighted)

names(time_3)[4] <- "event_id"
names(time_3)[5] <- "cue"
names(time_3)[6] <- "order_all"
names(time_3)[7] <- "order_cue"
names(time_3)[8] <- "sighted"

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
                 dplyr::rename (hr = V1, min = V2, sec = V3) %>%
                 dplyr::mutate(sec_frac = sec/60, min_sec = sec_frac+min) %>%
                 dplyr::mutate(min_frac = min_sec/60, time = hr+min_frac)
                         
j_dat_time_diff <- time_4 %>%
                              filter(event_id != "11_g") %>%
                              group_by(event_id) %>%
                              mutate(time_diff = lead(time) - time) %>%
                              mutate(time_diff_sec = time_diff*3600) %>%
                              mutate(time_diff_sec_shift = lag(time_diff_sec)) %>%
                              dplyr::select(event_id, cue, order_all, order_cue, time_diff_sec, 
                              sighted, time_diff_sec_shift) %>%
                              as.data.frame()

#summary(j_dat_time_diff)
# mean   : 27.80
# n = 146
# sd = 34.6

final_dat <- full_join(j_dat_time_diff, j_dat_tmp1, by = c("event_id", "order_all", "order_cue")


long_dives <- j_dat_time_diff %>%
                       filter(time_diff_sec_shift > 60)