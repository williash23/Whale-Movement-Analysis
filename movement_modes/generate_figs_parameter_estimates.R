# Sara Williams
# 2/22/2016; 
# Generate figures from parameter estimates from movement state models
#   and run in JAGS.
################################################################################

#  Load packages
library(CircStats)
library(reshape)
library(ggplot2)
library(dplyr)
library(rjags)

################################################################################
# Load output from MCMC posterior distribtuion iterations (JAGS object)

load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/single_fit.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/single_fit_cone.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/single_fit_transit.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/single_fit_station.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/double_fit.RData")
################################################################################

# NOTE: due to differences in formulation of Weibull distribution between JAGS and R:
###### scale = (1/lambda)^(1/v) ######

#  Single movement mode model
    # Insert model output object
    x <- single_fit
    niter <- 15000
    keep_1 <- sample(1:niter, 15000, replace = F)
    keep_2 <- sample(1:niter, 15000, replace = F)
    keep_3 <- sample(1:niter, 15000, replace = F)

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

#  Save simulated step length and turning angle values generated from MCMC posterior distribtions for later use.
save(post_sims, file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/post_sims.RData")

#  Two movement modes model
    #   Select rows from posterior distribution for each parameter estimate.
    #   model_fit_object[rows, columns, sheets]
    niter <- 10000
    x_t <- single_fit_transit
    x_s <- single_fit_station
    keep_1 <- sample(1:niter, 2000, replace = F)
    keep_2 <- sample(1:niter, 2000, replace = F)
    keep_3 <- sample(1:niter, 2000, replace = F)

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
      
          for(i in 1:nrow(sims_t)){
            steps_t[i] <- rweibull(1, sims_t[i,3], (1/sims_t[i,4])^(1/sims_t[i,3]))
            turns_t[i] <- turns[i] <- rwcauchy(1, sims_t[i,1], sims_t[i,2])
            #rwrappedcauchy(1, sims_t[i,1],  sims_t[i,2])
            }
          for(i in 1:nrow(sims_s)){
            steps_s[i] <- rweibull(1, sims_s[i,3], (1/sims_s[i,4])^(1/sims_s[i,3]))
            turns_s[i] <- turns[i] <- rwcauchy(1, sims_s[i,1], sims_s[i,2])
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
steps_plot <- ggplot(post_sims_plot) + 
                          geom_density(aes(steps), fill = "#EBCC2A", colour = "#EBCC2A") +
                          xlab("Step length (km)") +
                          ylab("Frequency") +
                          xlim(c(0, 10)) +
                          ylim(c(0, 1.5)) +
                          theme_bw() +
                          theme(panel.grid.minor = element_blank())
steps_plot

steps_plot <- ggplot(post_sims_plot) + 
                          geom_density(aes(steps), fill = "grey47", colour = "grey47") +
                          xlab("Step length (km)") +
                          ylab("Frequency") +
                          xlim(c(0, 5)) +
                          ylim(c(0, 1.5)) +
                          theme_bw() +
                          theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                          axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
steps_plot

### Colors: #78B7C5 #EBCC2A

steps_both_plot <- ggplot(both, aes(steps, colour = type, fill = type)) + 
                               geom_density(alpha= 0.7) +
                               xlab("Step length (km)") +
                               ylab("Frequency") +
                               xlim(c(0, 5)) +
                               ylim(c(0, 2.5)) +
                               theme_bw() +
                               theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                               axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
steps_both_col <- steps_both_plot + 
                              scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                              scale_fill_manual(values = c("grey47","grey87"), guide = FALSE)
steps_both_col

#  Point plot
steps_pt <- ggplot(post_sims_plot) +
                   geom_point(aes(steps, iter_num), colour = "grey47", alpha = 0.3) +
                   xlab("Step length (km)") +
                   ylab("Iteration number") +
                   theme_bw() +
                   theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                   axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
steps_pt

steps_pt_both <- ggplot(both) +
                             geom_point(aes(steps, iter_num), colour = type, alpha = 0.3) +
                             xlab("Step length (km)") +
                             ylab("Iteration number") +
                             theme_bw() +
                             theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                             axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
steps_pt_both_col <- steps_pt_both +
                                   scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                                   scale_fill_manual(values = c("grey47","grey87"), guide = FALSE)

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

#  Point plot
turns_pt <- ggplot(post_sims_plot, aes(turns, iter_num)) +
                   geom_point(colour = "black", alpha = 0.3) +
                   xlab("Turn angle (rad)") +
                   ylab("Iteration number") +
                   theme_bw()
turns_pt
################################################################################

#  Generate many simulations of N steps
#  Identify which models posterior iterations to use
mod <- post_sims
#   Number of movement paths
npath <- 10
#   Length of walk
nobs <- 3
#   Matrices to hold simulations
mat_X <- matrix(nrow = nobs, ncol = npath)
mat_Y <- matrix(nrow = nobs, ncol = npath)

s <- 0.65
e <- 0.85
dirs <- seq(s, e, ((e-s)/(npath-1)))

for(j in 1:npath){
keep <- sample(1:nrow(mod), nobs, replace = T)
#keep <- keep_start:(keep_start+nsteps-1)

   # make distributed steps
   steps_sim <- mod[keep, 1]
   # make clustered turning angles
   theta_sim <- mod[keep, 2]
   theta_sim[1] <- 0 ## control for first 2 turn angles so that each path starts in same direction
   theta_sim[2] <- 0.79 #dirs[j]
   # cumulative angle (absolute orientation)
   phi_sim <- cumsum(theta_sim)
   # step length components -- 
   dX_sim <- steps_sim*cos(phi_sim)
   dX_sim[1] <- 0 ##MAKE EACH START AT [0,0]
   dY_sim <- steps_sim*sin(phi_sim)
   
   # actual X-Y values
   X_sim <- as.matrix(cumsum(dX_sim))
   Y_sim <- as.matrix(cumsum(dY_sim))
   mat_X[,j] <- X_sim
   mat_Y[,j] <- Y_sim
   }

#  Data from for locations.
df_X_tmp <- as.data.frame(mat_X)
df_Y_tmp <- as.data.frame(mat_Y)
df_X <- melt(df_X_tmp)
df_Y <- melt(df_Y_tmp)
df_XY_tmp <- cbind(df_X , df_Y, rep(1:nobs))
names(df_XY_tmp) <- c("walk_num", "X", "walk_num_rep", "Y", "location_num")
df_XY <- df_XY_tmp %>%
               dplyr::select(walk_num, location_num, X, Y) %>%
               mutate(model = "Single State")

#  Figure of all simulated locations/paths
sim_crw <- ggplot() + 
                   geom_path(data = df_XY, alpha= 0.3, aes(x = X, y = Y, group = walk_num), colour = "grey47", size = 0.75) +
                   #geom_path(data = df_XY2, aes(x = X, y = Y, group = walk_num)) +
                   xlab("X (km)") +
                   xlim(c(-5, 5)) +
                   ylab("Y (km)") +
                   ylim(c(-5, 5)) +
                   theme_bw() +
                   theme(panel.grid.minor = element_blank(), axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
sim_crw
################################################################################

library(emojifont)
load.emojifont('OpenSansEmoji.ttf')

 # post_sims_95 <- post_sims %>%
                           # filter(steps < (quantile(post_sims$steps, probs = 0.95))
                                    # & steps > (quantile(post_sims$steps, probs = 0.05))) %>%
                           # filter(turns < (quantile(post_sims$turns, probs = 0.90))
                                    # & turns > (quantile(post_sims$turns, probs = 0.10))
                                    
### Show all possible turn angles but show variation in step length
circ <- seq(2.80, 4.71, 0.01) #### 180 degree arc of interest from ship perspective
#circ <- seq(0.54, 5.74, 0.01) #### 90% credible interval of turn angles
#  Identify which models posterior iterations to use
mod <- post_sims
#   Number of movement paths
npath <- length(circ)
#   Length of walk
nobs <- 2
#   Matrices to hold simulations
mat_X <- matrix(nrow = nobs, ncol = npath)
mat_Y <- matrix(nrow = nobs, ncol = npath)

for(i in 1:npath){
     keep <- sample(1:nrow(mod), nobs, replace = FALSE)
    #keep <- keep_start:(keep_start+nsteps-1)

       # make distributed steps
       steps_sim <- mod[keep, 1]
       # make clustered turning angles
       theta_sim <- mod[keep, 2]
       theta_sim[1] <- 0
       theta_sim[2] <- 0.785
       #theta_sim[2] <- circ[i]
       # cumulative angle (absolute orientation)
       phi_sim <- cumsum(theta_sim)
       # step length components -- 
       dX_sim <- steps_sim*cos(phi_sim)
       dX_sim[1] <- 0 ##MAKE EACH START AT [0,0]
       dY_sim <- steps_sim*sin(phi_sim)
       
       # actual X-Y values
       X_sim <- as.matrix(cumsum(dX_sim))
       Y_sim <- as.matrix(cumsum(dY_sim))
       mat_X[,i] <- X_sim
       mat_Y[,i] <- Y_sim
   }

#  Data from for locations.
df_X_tmp <- as.data.frame(mat_X)
df_Y_tmp <- as.data.frame(mat_Y)
df_X <- melt(df_X_tmp)
df_Y <- melt(df_Y_tmp)
df_XY_tmp <- cbind(df_X , df_Y, rep(1:nobs))
names(df_XY_tmp) <- c("walk_num", "X", "walk_num_rep", "Y", "location_num")
df_XY <- df_XY_tmp %>%
               dplyr::select(walk_num, location_num, X, Y) %>%
               mutate(model = "Single State")


ship_line <- data.frame(x1 = -2, x2 = -2, y1 = 3.5, y2 = 0.5)
ship_loc1 <- data.frame(x = -2, y = 3.8, label = (emoji('ship')))
ship_loc2 <- data.frame(x = -2, y = 0.2, label = (emoji('ship')))
whale_loc <- data.frame(x = 0.4, y = 0, label = (emoji('whale')))

#  Figure of all simulated locations/paths and ship moving - 7.5 min at 13 kts == 3 km of ship movement.
sim_step_crw <- ggplot() + 
                            geom_polygon(data = df_XY, aes(x = X, y = Y, group = walk_num), colour = "gray47" fill = "gray47") +
                            #geom_path(data = df_XY2, aes(x = X, y = Y, group = walk_num)) +
                            xlab("X (km)") +
                            xlim(c(-8, 8)) +
                            ylab("Y (km)") +
                            ylim(c(6, -6)) +
                            theme_bw() +
                            theme(panel.grid.minor = element_blank(), panel.border = element_blank(), axis.text = element_text(size = rel(2)), 
                            axis.title = element_text(size = rel(2))) +
                            geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), colour = "black", size = 1.5, 
                            data = ship_line, arrow = arrow(length = unit(0.15, "inches"), type = "closed")) +
                            geom_text(data = ship_loc1, aes(x, y, label=label), family="OpenSansEmoji", size=15, colour = "darkred") +
                            geom_text(data = ship_loc2, aes(x, y, label=label), family="OpenSansEmoji", size=15, colour = "darkred", alpha = 0.2) +
                            geom_text(data = whale_loc, aes(x, y, label=label), family="OpenSansEmoji", size=14, colour="darkblue")
sim_step_crw 



ship_line <- data.frame(x1 = -2, x2 = -2, y1 = 4.7, y2 = 1)
ship_loc <- data.frame(x=-2, y=5, label = (emoji('ship')))
whale_loc <- data.frame(x = 0.4, y = 0, label = (emoji('whale')))

#  Figure of all simulated locations/paths and ship
sim_step_crw <- ggplot() + 
                            geom_path(data = df_XY, alpha= 0.05, aes(x = X, y = Y, group = walk_num), colour = "gray47", size = 1) +
                            #geom_path(data = df_XY2, aes(x = X, y = Y, group = walk_num)) +
                            xlab("X (km)") +
                            xlim(c(-7, 3)) +
                            ylab("Y (km)") +
                            ylim(c(5, -5)) +
                            theme_bw() +
                            theme(panel.grid.minor = element_blank(), axis.text = element_text(size = rel(2)), 
                            axis.title = element_text(size = rel(2))) +
                            geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), colour = "black", size = 1.5, 
                            data = ship_line, arrow = arrow(length = unit(0.15, "inches"), type = "closed")) +
                            geom_text(data = ship_loc, aes(x, y, label=label), family="OpenSansEmoji", size=12, colour = "darkred") +
                            geom_text(data = whale_loc, aes(x, y, label=label), family="OpenSansEmoji", size=15, colour="darkblue")
sim_step_crw 

################################################################################


#  Function from J. Nowak: summarise mcmc runs
#   Takes a mcmc.list as called by PopR with year attached, rnd is a 
#   logical input when T values < 1 rounded to 2 digits while larger 
#   values are rounded to zero digits and static_covars is a  
#   character vector of covariate names that do not vary in time
#   Returns a data frame summarising each parameter by mean, sd and 95%
#   credible interval boundaries

 jags_reduce <- function(x, rnd = T, static_covars = NULL){
        tmp <- as.matrix(x)

        summ <- data.frame(Parameter = colnames(tmp)) %>%
          mutate(
            Mean = apply(tmp, 2, mean, na.rm = T),
            SD = apply(tmp, 2, sd, na.rm = T),
            lo = apply(tmp, 2, quantile, 0.025, na.rm = T),
            hi = apply(tmp, 2, quantile, 0.975, na.rm = T),
            Parameter = gsub("\\[.*\\]", "", Parameter)) 

        if(rnd){
          summ <- summ %>%
           group_by(Parameter) %>%
           mutate(Mean = round(Mean, ifelse(mean(Mean) < 3, 2, 0)),
              lo = round(lo, ifelse(mean(Mean) < 3, 2, 0)),
              hi = round(hi, ifelse(mean(Mean) < 3, 2, 0)),
              SD = round(SD, ifelse(mean(Mean) < 3, 2, 0))) %>%
            ungroup()
        }

       rownames(summ) <- NULL

    return(summ)
    }

jags_reduce(single_fit)
# jags_reduce(single_fit_transit)
# jags_reduce(single_fit_station)
# jags_reduce(single_fit_cone)
# ################################################################################
