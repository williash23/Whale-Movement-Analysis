# Sara Williams
# 2/22/2016; 
# Generate figures from parameter estimates from movement state models
#   and run in JAGS. For data seperated into surface interval and deep dive by time between
#   observations.
################################################################################

#  Load packages
library(CircStats)
library(reshape)
library(ggplot2)
library(dplyr)
library(rjags)
library(emojifont)

################################################################################
# Load output from MCMC posterior distribtuion iterations (JAGS object)
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/single_fit.RData")

load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/single_fit_transit.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/single_fit_station.RData")

load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/new_categories_output/single_fit_close.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/new_categories_output/single_fit_far.RData")

load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/new_categories_output/single_fit_side.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/new_categories_output/single_fit_front.RData")

load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/new_categories_output/single_fit_dive.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/new_categories_output/single_fit_surf.RData")
################################################################################

# NOTE: due to differences in formulation of Weibull distribution between JAGS and R:
###### scale = (1/lambda)^(1/v) ######

#  Deep dive vs surface interval
    #   Select rows from posterior distribution for each parameter estimate.
    #   model_fit_object[rows, columns, sheets]
    niter <- 10000
    x_d <- single_fit_dive
    x_su <- single_fit_surf
    keep_1 <- sample(1:niter, 5000, replace = F)
    keep_2 <- sample(1:niter, 5000, replace = F)
    keep_3 <- sample(1:niter, 5000, replace = F)

          chain_1_d <- x_d[[1]]
          sims_1_d <- chain_1_d[keep_1, c(2, 3, 4, 1)]
          chain_1_su<- x_su[[1]]
          sims_1_su <- chain_1_su[keep_1, c(2, 3, 4, 1)]
          
          chain_2_d <- x_d[[2]]
          sims_2_d <- chain_2_d[keep_2, c(2, 3, 4, 1)]
          chain_2_su <- x_su[[2]]
          sims_2_su <- chain_2_su[keep_2, c(2, 3, 4, 1)]
          
          chain_3_d <- x_d[[3]]
          sims_3_d <- chain_3_d[keep_3, c(2, 3, 4, 1)]
          chain_3_su <- x_su[[3]]
          sims_3_su <- chain_3_su[keep_3, c(2, 3, 4, 1)]

          sims_d <- rbind(sims_1_d, sims_2_d, sims_3_d)
          sims_su <- rbind(sims_1_su, sims_2_su, sims_3_su)

          steps_d <- numeric(length = nrow(sims_d))
          turns_d <- numeric(length = nrow(sims_d))
          steps_su <- numeric(length = nrow(sims_su))
          turns_su <- numeric(length = nrow(sims_su))
          
           rwcauchy <- function(n, mu = 0, rho = 0) {
              u = runif(n)
              V = cos(2 * pi * u)
              c = 2 * rho/(1 + rho^2)
              t <- (sign(runif(n) - 0.5) * acos((V + c)/(1 + c * V)) + mu)%%(2 * pi)
              return(t)
          }
      
          for(i in 1:nrow(sims_d)){
            steps_d[i] <- rweibull(1, sims_d[i,3], (1/sims_d[i,4])^(1/sims_d[i,3]))
            turns_d[i] <- rwcauchy(1, sims_d[i,1], sims_d[i,2])
            #rwrappedcauchy(1, sims_d[i,1],  sims_d[i,2])
            }
          for(i in 1:nrow(sims_su)){
            steps_su[i] <- rweibull(1, sims_su[i,3], (1/sims_su[i,4])^(1/sims_su[i,3]))
            turns_su[i] <- rwcauchy(1, sims_su[i,1], sims_su[i,2])
            #rwrappedcauchy(1, sims_su[i,1],  sims_su[i,2])
            }

        post_sims_d <- as.data.frame(cbind(steps_d, turns_d))
        post_sims_d$turns_d[post_sims_d$turns_d>pi]=post_sims_d$turns_d[post_sims_d$turns_d>pi]-2*pi
        post_sims_su <- as.data.frame(cbind(steps_su, turns_su))
        post_sims_su$turns_su[post_sims_su$turns_su>pi]=post_sims_su$turns_su[post_sims_su$turns_su>pi]-2*pi
        
        post_sims_d_plot <- post_sims_d %>% 
                                          mutate(iter_num = 1:nrow(post_sims_d)) %>%
                                          mutate(type = "deep dive")  %>%
                                          dplyr::rename(steps = steps_d, turns = turns_d)
        post_sims_su_plot <- post_sims_su %>% 
                                          mutate(iter_num = 1:nrow(post_sims_su)) %>%
                                          mutate(type = "surface interval") %>%
                                          dplyr::rename(steps = steps_su, turns = turns_su)
        both <- rbind(post_sims_d_plot, post_sims_su_plot)
        both$type <- as.factor(both$type)
################################################################################      

#  Density plot
### Colors: #78B7C5 #EBCC2A
steps_both_plot <- ggplot(both, aes(steps, colour = type, fill = type)) + 
                               geom_density(alpha= 0.7) +
                               xlab("Step length (km)") +
                               ylab("Frequency") +
                               #xlim(c(0, 5)) +
                               #ylim(c(0, 2.5)) +
                               theme_bw() +
                               theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                               axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
steps_both_col <- steps_both_plot + 
                              scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                              scale_fill_manual(values = c("grey47","grey87"))
steps_both_col

# steps_pt_both <- ggplot(both, aes(steps, iter_num, colour = type)) +
                             # geom_point(alpha = 0.3) +
                             # xlab("Step length (km)") +
                             # ylab("Iteration number") +
                             # theme_bw() +
                             # theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                             # axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
# steps_pt_both_col <- steps_pt_both +
                                   # scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                                   # scale_fill_manual(values = c("grey47","grey87"), guide = FALSE)

#  Density plot
turns_both_plot <- ggplot(both, aes(turns, colour = type, fill = type)) + 
                               geom_density(alpha= 0.7) +
                               xlab("Turn angle (rad)") +
                               ylab("Frequency") +
                               xlim(c(-4, 4)) +
                               theme_bw() +
                               theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                               axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
turns_both_col <- turns_both_plot + 
                             scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                              scale_fill_manual(values = c("grey47","grey87"))
turns_both_col

# #  Point plot
# turns_pt_both <- ggplot(both, aes(turns, iter_num, colour = type)) +
                             # geom_point(alpha = 0.3) +
                             # xlab("Turn angle (rad)") +
                             # ylab("Iteration number") +
                             # theme_bw() +
                             # theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                             # axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
# turns_pt_both_col <- turns_pt_both +
                                   # scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                                   # scale_fill_manual(values = c("grey47","grey87"), guide = FALSE)
################################################################################



#  First observation close to ship or far from ship
    #   Select rows from posterior distribution for each parameter estimate.
    #   model_fit_object[rows, columns, sheets]
    niter <- 10000
    x_c <- single_fit_close
    x_f <- single_fit_far
    keep_1 <- sample(1:niter, 5000, replace = F)
    keep_2 <- sample(1:niter, 5000, replace = F)
    keep_3 <- sample(1:niter, 5000, replace = F)

          chain_1_c <- x_c[[1]]
          sims_1_c <- chain_1_c[keep_1, c(2, 3, 4, 1)]
          chain_1_f<- x_f[[1]]
          sims_1_f <- chain_1_f[keep_1, c(2, 3, 4, 1)]
          
          chain_2_c <- x_c[[2]]
          sims_2_c <- chain_2_c[keep_2, c(2, 3, 4, 1)]
          chain_2_f <- x_f[[2]]
          sims_2_f <- chain_2_f[keep_2, c(2, 3, 4, 1)]
          
          chain_3_c <- x_c[[3]]
          sims_3_c <- chain_3_c[keep_3, c(2, 3, 4, 1)]
          chain_3_f <- x_f[[3]]
          sims_3_f <- chain_3_f[keep_3, c(2, 3, 4, 1)]

          sims_c <- rbind(sims_1_c, sims_2_c, sims_3_c)
          sims_f <- rbind(sims_1_f, sims_2_f, sims_3_f)

          steps_c <- numeric(length = nrow(sims_c))
          turns_c <- numeric(length = nrow(sims_c))
          steps_f <- numeric(length = nrow(sims_f))
          turns_f <- numeric(length = nrow(sims_f))
          
           rwcauchy <- function(n, mu = 0, rho = 0) {
              u = runif(n)
              V = cos(2 * pi * u)
              c = 2 * rho/(1 + rho^2)
              t <- (sign(runif(n) - 0.5) * acos((V + c)/(1 + c * V)) + mu)%%(2 * pi)
              return(t)
          }
      
          for(i in 1:nrow(sims_c)){
            steps_c[i] <- rweibull(1, sims_c[i,3], (1/sims_c[i,4])^(1/sims_c[i,3]))
            turns_c[i] <- rwcauchy(1, sims_c[i,1], sims_c[i,2])
            #rwrappedcauchy(1, sims_c[i,1],  sims_c[i,2])
            }
          for(i in 1:nrow(sims_f)){
            steps_f[i] <- rweibull(1, sims_f[i,3], (1/sims_f[i,4])^(1/sims_f[i,3]))
            turns_f[i] <- rwcauchy(1, sims_f[i,1], sims_f[i,2])
            #rwrappedcauchy(1, sims_f[i,1],  sims_f[i,2])
            }

        post_sims_c <- as.data.frame(cbind(steps_c, turns_c))
        post_sims_c$turns_c[post_sims_c$turns_c>pi]=post_sims_c$turns_c[post_sims_c$turns_c>pi]-2*pi
        post_sims_f <- as.data.frame(cbind(steps_f, turns_f))
        post_sims_f$turns_f[post_sims_f$turns_f>pi]=post_sims_f$turns_f[post_sims_f$turns_f>pi]-2*pi
        
        post_sims_c_plot <- post_sims_c %>% 
                                          mutate(iter_num = 1:nrow(post_sims_c)) %>%
                                          mutate(type = "first observation close to ship")  %>%
                                          dplyr::rename(steps = steps_c, turns = turns_c)
        post_sims_f_plot <- post_sims_f %>% 
                                          mutate(iter_num = 1:nrow(post_sims_f)) %>%
                                          mutate(type = "first observation far from ship") %>%
                                          dplyr::rename(steps = steps_f, turns = turns_f)
        both <- rbind(post_sims_c_plot, post_sims_f_plot)
        both$type <- as.factor(both$type)
################################################################################      

#  Density plot
### Colors: #78B7C5 #EBCC2A
steps_both_plot <- ggplot(both, aes(steps, colour = type, fill = type)) + 
                               geom_density(alpha= 0.7) +
                               xlab("Step length (km)") +
                               ylab("Frequency") +
                               #xlim(c(0, 5)) +
                               #ylim(c(0, 2.5)) +
                               theme_bw() +
                               theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                               axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
steps_both_col <- steps_both_plot + 
                              scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                              scale_fill_manual(values = c("grey47","grey87"))
steps_both_col

# steps_pt_both <- ggplot(both, aes(steps, iter_num, colour = type)) +
                             # geom_point(alpha = 0.3) +
                             # xlab("Step length (km)") +
                             # ylab("Iteration number") +
                             # theme_bw() +
                             # theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                             # axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
# steps_pt_both_col <- steps_pt_both +
                                   # scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                                   # scale_fill_manual(values = c("grey47","grey87"), guide = FALSE)

#  Density plot
turns_both_plot <- ggplot(both, aes(turns, colour = type, fill = type)) + 
                               geom_density(alpha= 0.7) +
                               xlab("Turn angle (rad)") +
                               ylab("Frequency") +
                               xlim(c(-4, 4)) +
                               theme_bw() +
                               theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                               axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
turns_both_col <- turns_both_plot + 
                             scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                              scale_fill_manual(values = c("grey47","grey87"))
turns_both_col

# #  Point plot
# turns_pt_both <- ggplot(both, aes(turns, iter_num, colour = type)) +
                             # geom_point(alpha = 0.3) +
                             # xlab("Turn angle (rad)") +
                             # ylab("Iteration number") +
                             # theme_bw() +
                             # theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                             # axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
# turns_pt_both_col <- turns_pt_both +
                                   # scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                                   # scale_fill_manual(values = c("grey47","grey87"), guide = FALSE)
################################################################################



#  First observation at side of ship or at front of ship
    #   Select rows from posterior distribution for each parameter estimate.
    #   model_fit_object[rows, columns, sheets]
    niter <- 10000
    x_si <- single_fit_side
    x_fr <- single_fit_front
    keep_1 <- sample(1:niter, 5000, replace = F)
    keep_2 <- sample(1:niter, 5000, replace = F)
    keep_3 <- sample(1:niter, 5000, replace = F)

          chain_1_si <- x_si[[1]]
          sims_1_si <- chain_1_si[keep_1, c(2, 3, 4, 1)]
          chain_1_fr<- x_fr[[1]]
          sims_1_fr <- chain_1_fr[keep_1, c(2, 3, 4, 1)]
          
          chain_2_si <- x_si[[2]]
          sims_2_si <- chain_2_si[keep_2, c(2, 3, 4, 1)]
          chain_2_fr <- x_fr[[2]]
          sims_2_fr <- chain_2_fr[keep_2, c(2, 3, 4, 1)]
          
          chain_3_si <- x_si[[3]]
          sims_3_si <- chain_3_si[keep_3, c(2, 3, 4, 1)]
          chain_3_fr <- x_fr[[3]]
          sims_3_fr <- chain_3_fr[keep_3, c(2, 3, 4, 1)]

          sims_si <- rbind(sims_1_si, sims_2_si, sims_3_si)
          sims_fr <- rbind(sims_1_fr, sims_2_fr, sims_3_fr)

          steps_si <- numeric(length = nrow(sims_si))
          turns_si <- numeric(length = nrow(sims_si))
          steps_fr <- numeric(length = nrow(sims_fr))
          turns_fr <- numeric(length = nrow(sims_fr))
          
           rwcauchy <- function(n, mu = 0, rho = 0) {
              u = runif(n)
              V = cos(2 * pi * u)
              c = 2 * rho/(1 + rho^2)
              t <- (sign(runif(n) - 0.5) * acos((V + c)/(1 + c * V)) + mu)%%(2 * pi)
              return(t)
          }
      
          for(i in 1:nrow(sims_si)){
            steps_si[i] <- rweibull(1, sims_si[i,3], (1/sims_si[i,4])^(1/sims_si[i,3]))
            turns_si[i] <- rwcauchy(1, sims_si[i,1], sims_si[i,2])
            #rwrappedcauchy(1, sims_si[i,1],  sims_si[i,2])
            }
          for(i in 1:nrow(sims_fr)){
            steps_fr[i] <- rweibull(1, sims_fr[i,3], (1/sims_fr[i,4])^(1/sims_fr[i,3]))
            turns_fr[i] <- rwcauchy(1, sims_fr[i,1], sims_fr[i,2])
            #rwrappedcauchy(1, sims_fr[i,1],  sims_fr[i,2])
            }

        post_sims_si <- as.data.frame(cbind(steps_si, turns_si))
        post_sims_si$turns_si[post_sims_si$turns_si>pi]=post_sims_si$turns_si[post_sims_si$turns_si>pi]-2*pi
        post_sims_fr <- as.data.frame(cbind(steps_fr, turns_fr))
        post_sims_fr$turns_fr[post_sims_fr$turns_fr>pi]=post_sims_fr$turns_fr[post_sims_fr$turns_fr>pi]-2*pi
        
        post_sims_si_plot <- post_sims_si %>% 
                                          mutate(iter_num = 1:nrow(post_sims_si)) %>%
                                          mutate(type = "first observation at side of ship")  %>%
                                          dplyr::rename(steps = steps_si, turns = turns_si)
        post_sims_fr_plot <- post_sims_fr %>% 
                                          mutate(iter_num = 1:nrow(post_sims_fr)) %>%
                                          mutate(type = "first observation at front of ship") %>%
                                          dplyr::rename(steps = steps_fr, turns = turns_fr)
        both <- rbind(post_sims_si_plot, post_sims_fr_plot)
        both$type <- as.factor(both$type)
################################################################################      

#  Density plot
### Colors: #78B7C5 #EBCC2A
steps_both_plot <- ggplot(both, aes(steps, colour = type, fill = type)) + 
                               geom_density(alpha= 0.7) +
                               xlab("Step length (km)") +
                               ylab("Frequency") +
                               #xlim(c(0, 5)) +
                               #ylim(c(0, 2.5)) +
                               theme_bw() +
                               theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                               axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
steps_both_col <- steps_both_plot + 
                              scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                              scale_fill_manual(values = c("grey47","grey87"))
steps_both_col

# steps_pt_both <- ggplot(both, aes(steps, iter_num, colour = type)) +
                             # geom_point(alpha = 0.3) +
                             # xlab("Step length (km)") +
                             # ylab("Iteration number") +
                             # theme_bw() +
                             # theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                             # axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
# steps_pt_both_col <- steps_pt_both +
                                   # scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                                   # scale_fill_manual(values = c("grey47","grey87"), guide = FALSE)

#  Density plot
turns_both_plot <- ggplot(both, aes(turns, colour = type, fill = type)) + 
                               geom_density(alpha= 0.7) +
                               xlab("Turn angle (rad)") +
                               ylab("Frequency") +
                               xlim(c(-4, 4)) +
                               theme_bw() +
                               theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                               axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
turns_both_col <- turns_both_plot + 
                             scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                             scale_fill_manual(values = c("grey47","grey87"))
turns_both_col

# #  Point plot
# turns_pt_both <- ggplot(both, aes(turns, iter_num, colour = type)) +
                             # geom_point(alpha = 0.3) +
                             # xlab("Turn angle (rad)") +
                             # ylab("Iteration number") +
                             # theme_bw() +
                             # theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                             # axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
# turns_pt_both_col <- turns_pt_both +
                                   # scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                                   # scale_fill_manual(values = c("grey47","grey87"), guide = FALSE)
################################################################################


#  All data in single movement mode model
    #  Single movement mode model
    # Insert model output object
    x <- single_fit
    niter <- 10000
    keep_1 <- sample(1:niter, 5000, replace = F)
    keep_2 <- sample(1:niter, 5000, replace = F)
    keep_3 <- sample(1:niter, 5000, replace = F)

    #  One movement mode:
    #   Select rows from posterior distribution for each parameter estimate.
    #   model_fit_object[rows, columns, sheets]
          chain_1 <- x[[1]]
          sims_1 <- chain_1[keep_1, c(2, 3, 4, 1)]
          
          chain_2 <- x[[2]]
          sims_2 <- chain_2[keep_2, c(2, 3, 4, 1)]
          
          chain_3 <- x[[3]]
          sims_3 <- chain_3[keep_3, c(2, 3, 4, 1)]
        
          sims <- rbind(sims_1, sims_2, sims_3)

          steps <- numeric(length = nrow(sims))
          turns <- numeric(length = nrow(sims))
          
          rwcauchy <- function(n, mu = 0, rho = 0) {
              u = runif(n)
              V = cos(2 * pi * u)
              c = 2 * rho/(1 + rho^2)
              t <- (sign(runif(n) - 0.5) * acos((V + c)/(1 + c * V)) + mu)%%(2 * pi)
              return(t)
          }

          for(i in 1:nrow(sims)){
            steps[i] <- rweibull(1, sims[i,3], (1/sims[i,4])^(1/sims[i,3]))
            turns[i] <- rwcauchy(1, sims[i,1], sims[i,2])
            #rwrappedcauchy(1, sims[i,1],  sims[i,2])
            
            }

        post_sims <- as.data.frame(cbind(steps, turns))
        post_sims$turns[post_sims$turns>pi]=post_sims$turns[post_sims$turns>pi]-2*pi

        post_sims_plot <- post_sims %>% 
                                          mutate(iter_num = 1:nrow(post_sims))


#  All data - two single movement modes model
    #   Select rows from posterior distribution for each parameter estimate.
    #   model_fit_object[rows, columns, sheets]
    niter <- 10000
    x_t <- single_fit_transit
    x_s <- single_fit_station
    keep_1 <- sample(1:niter, 5000, replace = F)
    keep_2 <- sample(1:niter, 5000, replace = F)
    keep_3 <- sample(1:niter, 5000, replace = F)

          chain_1_t <- x_t[[1]]
          sims_1_t <- chain_1_t[keep_1, c(2, 3, 4, 1)]
          chain_1_s<- x_s[[1]]
          sims_1_s <- chain_1_s[keep_1, c(2, 3, 4, 1)]
          
          chain_2_t <- x_t[[2]]
          sims_2_t <- chain_2_t[keep_2, c(2, 3, 4, 1)]
          chain_2_s <- x_s[[2]]
          sims_2_s <- chain_2_s[keep_2, c(2, 3, 4, 1)]
          
          chain_3_t <- x_t[[3]]
          sims_3_t <- chain_3_t[keep_3, c(2, 3, 4, 1)]
          chain_3_s <- x_s[[3]]
          sims_3_s <- chain_3_s[keep_3, c(2, 3, 4, 1)]

          sims_t <- rbind(sims_1_t, sims_2_t, sims_3_t)
          sims_s <- rbind(sims_1_s, sims_2_s, sims_3_s)

          steps_t <- numeric(length = nrow(sims_t))
          turns_t <- numeric(length = nrow(sims_t))
          steps_s <- numeric(length = nrow(sims_s))
          turns_s <- numeric(length = nrow(sims_s))
          
           rwcauchy <- function(n, mu = 0, rho = 0) {
              u = runif(n)
              V = cos(2 * pi * u)
              c = 2 * rho/(1 + rho^2)
              t <- (sign(runif(n) - 0.5) * acos((V + c)/(1 + c * V)) + mu)%%(2 * pi)
              return(t)
          }
      
          for(i in 1:nrow(sims_t)){
            steps_t[i] <- rweibull(1, sims_t[i,3], (1/sims_t[i,4])^(1/sims_t[i,3]))
            turns_t[i] <- rwcauchy(1, sims_t[i,1], sims_t[i,2])
            #rwrappedcauchy(1, sims_t[i,1],  sims_t[i,2])
            }
          for(i in 1:nrow(sims_s)){
            steps_s[i] <- rweibull(1, sims_s[i,3], (1/sims_s[i,4])^(1/sims_s[i,3]))
            turns_s[i] <- rwcauchy(1, sims_s[i,1], sims_s[i,2])
            #rwrappedcauchy(1, sims_s[i,1],  sims_s[i,2])
            }

        post_sims_t <- as.data.frame(cbind(steps_t, turns_t))
        post_sims_t$turns_t[post_sims_t$turns_t>pi]=post_sims_t$turns_t[post_sims_t$turns_t>pi]-2*pi
        post_sims_s <- as.data.frame(cbind(steps_s, turns_s))
        post_sims_s$turns_s[post_sims_s$turns_s>pi]=post_sims_s$turns_s[post_sims_s$turns_s>pi]-2*pi
        
        post_sims_t_plot <- post_sims_t %>% 
                                          mutate(iter_num = 1:nrow(post_sims_t)) %>%
                                          mutate(type = "transit")  %>%
                                          dplyr::rename(steps = steps_t, turns = turns_t)
        post_sims_s_plot <- post_sims_s %>% 
                                          mutate(iter_num = 1:nrow(post_sims_s)) %>%
                                          mutate(type = "station") %>%
                                          dplyr::rename(steps = steps_s, turns = turns_s)
        both <- rbind(post_sims_t_plot, post_sims_s_plot)
################################################################################      

#  Density plot
# steps_plot <- ggplot(post_sims_plot) + 
                          # geom_density(aes(steps), fill = "#EBCC2A", colour = "#EBCC2A") +
                          # xlab("Step length (km)") +
                          # ylab("Frequency") +
                          # xlim(c(0, 5)) +
                          # ylim(c(0, 2.5)) +
                          # theme_bw() +
                          # theme(panel.grid.minor = element_blank())
# steps_plot

steps_plot <- ggplot(post_sims_plot) + 
                          geom_density(aes(steps), fill = "grey47", colour = "grey47") +
                          xlab("Step length (km)") +
                          ylab("Frequency") +
                          #xlim(c(0, 5)) +
                          #ylim(c(0, 2.5)) +
                          theme_bw() +
                          theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                          axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
steps_plot

### Colors: #78B7C5 #EBCC2A

steps_both_plot <- ggplot(both, aes(steps, colour = type, fill = type)) + 
                               geom_density(alpha= 0.7) +
                               xlab("Step length (km)") +
                               ylab("Frequency") +
                               #xlim(c(0, 5)) +
                               #ylim(c(0, 2.5)) +
                               theme_bw() +
                               theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                               axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
steps_both_col <- steps_both_plot + 
                              scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                              scale_fill_manual(values = c("grey47","grey87"), guide = FALSE)
steps_both_col

# #  Point plot
# steps_pt <- ggplot(post_sims_plot) +
                   # geom_point(aes(steps, iter_num), colour = "grey47", alpha = 0.3) +
                   # xlab("Step length (km)") +
                   # ylab("Iteration number") +
                   # theme_bw() +
                   # theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                   # axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
# steps_pt

# steps_pt_both <- ggplot(both) +
                             # geom_point(aes(steps, iter_num), colour = type, alpha = 0.3) +
                             # xlab("Step length (km)") +
                             # ylab("Iteration number") +
                             # theme_bw() +
                             # theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                             # axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
# steps_pt_both_col <- steps_pt_both +
                                   # scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                                   # scale_fill_manual(values = c("grey47","grey87"), guide = FALSE)

#  Density plot
turns_plot <- ggplot(post_sims_plot) + 
                      geom_density(aes(turns), fill = "grey47", colour = "grey47") +
                      xlab("Turn angle (rad)") +
                      ylab("Frequency") +
                      xlim(c(-4, 4)) +
                      theme_bw() +
                      theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                      axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
turns_plot

turns_both_plot <- ggplot(both, aes(turns, colour = type, fill = type)) + 
                               geom_density(alpha= 0.7) +
                               xlab("Turn angle (rad)") +
                               ylab("Frequency") +
                               xlim(c(-4, 4)) +
                               theme_bw() +
                               theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                               axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
turns_both_col <- turns_both_plot + 
                             scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                              scale_fill_manual(values = c("grey47","grey87"))
turns_both_col

# #  Point plot
# turns_pt <- ggplot(post_sims_plot, aes(turns, iter_num)) +
                   # geom_point(colour = "black", alpha = 0.3) +
                   # xlab("Turn angle (rad)") +
                   # ylab("Iteration number") +
                   # theme_bw()
# turns_pt
################################################################################
