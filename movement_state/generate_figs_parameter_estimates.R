# Sara Williams
# 2/22/2016; 
# Generate figures from parameter estimates from movement state models
#   and run in JAGS.
################################################################################

#  Load packages
library(circular)
library(reshape)
library(ggplot2)
library(dplyr)
library(CircStats)

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



# Single model output
    x <- single_fit
    niter <- 10000
    keep_1 <- sample(1:niter, replace = F)
    keep_2 <- sample(1:niter, replace = F)
    keep_3 <- sample(1:niter, replace = F)

    #  Select rows from posterior distribution for each parameter estimate.
    #   single_fit[rows, columns, sheets]
          chain_1 <- x[[1]]
          sims_1 <- chain_1[keep_1, c(5, 4, 2, 3)]

          chain_2 <- x[[2]]
          sims_2 <- chain_2[keep_2, c(5, 4, 2, 3)]

          chain_3 <- x[[3]]
          sims_3 <- chain_3[keep_3, c(5, 4, 2, 3)]

          sims <- rbind(sims_1, sims_2, sims_3)

          steps <- numeric(length = nrow(sims))
          turns <- numeric(length = nrow(sims))
      
          for(i in 1:nrow(sims)){
            steps[i] <- rweibull(1, sims[i,1], sims[i,2])
            turns[i] <- rwrappedcauchy(1, sims[i,3],  sims[i,4])
            }

        post_iters <- as.data.frame(cbind(steps, turns))
        post_iters$turns[post_iters$turns>pi]=post_iters$turns[post_iters$turns>pi]-2*pi

        post_iters_plot <- post_iters %>% 
                                      mutate(iter_num = 1:nrow(post_iters))

#  Density plot
steps_plot <- ggplot(post_iters) + 
                      geom_density(aes(steps), fill = "light grey", colour = "grey") +
                      xlab("Step length (km)") +
                      ylab("Frequency") +
                      xlim(c(0, 4)) +
                      ylim(c(0, 2.5)) +
                      theme_bw() +
                      theme(panel.grid.minor = element_blank())
steps_plot
#  Point plot
steps_pt <- ggplot(post_iters_plot, aes(steps, iter_num)) +
                   geom_point(colour = "black", alpha = 0.3) +
                   xlab("Step length (km)") +
                   ylab("Iteration number") +
                   theme_bw()
steps_pt

#  Density plot
turns_plot <- ggplot(post_iters) + 
                      geom_density(aes(turns), fill = "light grey", colour = "grey") +
                      xlab("Turn angle (rad)") +
                      ylab("Frequency") +
                      xlim(c(-4, 4)) +
                      theme_bw() +
                      theme(panel.grid.minor = element_blank())
turns_plot
#  Point plot
turns_pt <- ggplot(post_iters_plot, aes(turns, iter_num)) +
                   geom_point(colour = "black", alpha = 0.3) +
                   xlab("Turn angle (rad)") +
                   ylab("Iteration number") +
                   theme_bw()
turns_pt

#  Generate many simulations of N steps of movement from state 1 and state 2
#   Number of simulations
nsims <- 10000
#   Length of walk
nsteps <- 3

#   Matrix to hold state 1 simulations
mat_X <- matrix(nrow = nsteps, ncol = nsims)
mat_Y <- matrix(nrow = nsteps, ncol = nsims)

for(j in 1:nsims){
keep <- sample(1:nrow(post_iters), nsteps, replace = F)

   # make distributed steps
   steps_sim <-post_iters[keep, 1]
   # make clustered turning angles
   theta_sim <-post_iters[keep, 2]
   # cumulative angle (absolute orientation)
   phi_sim <- cumsum(theta_sim)
   # step length components -- MAKE EACH START AT [0,0]
   dX_sim <- steps_sim*cos(phi_sim)
   dX_sim[1] <- 0
   dY_sim <- steps_sim*sin(phi_sim)
   dY_sim[1] <- 0
   # actual X-Y values
   X_sim <- as.matrix(cumsum(dX_sim))
   Y_sim <- as.matrix(cumsum(dY_sim))
   mat_X[,j] <- X_sim
   mat_Y[,j] <- Y_sim
}

#   Data from for locations in state 1 movement mode.
df_X_tmp <- as.data.frame(mat_X)
df_Y_tmp <- as.data.frame(mat_Y)
df_X <- melt(df_X_tmp)
df_Y <- melt(df_Y_tmp)
df_XY_tmp <- cbind(df_X , df_Y, rep(1:nsteps))
names(df_XY_tmp) <- c("walk_num", "X", "walk_num_rep", "Y", "location_num")
df_XY <- df_XY_tmp %>%
               dplyr::select(walk_num, location_num, X, Y) %>%
               mutate(model = "Single State")

#   Figure of all simulated locations/paths from state 1 and state 2 movement modes.
sim_crw <- ggplot(df_XY) + 
                    geom_line(aes(x = X, y = Y, group = walk_num), alpha = 0.1)+
                    xlab("X (km)") +
                    ylab("Y (km)") +
                    theme_bw()
sim_crw 

df_XY$X[df_XY$location_num == 2] <- 0.079597698
df_XY$Y[df_XY$location_num == 2] <- 0.08254940

























# Double model output
    niter <- 1000
    keep_1 <- sample(1:niter, replace = F)
    keep_2 <- sample(1:niter, replace = F)
    keep_3 <- sample(1:niter, replace = F)

    #  Select rows from posterior distribution for each parameter estimate.
    #   single_fit[rows, columns, sheets]
          chain_1 <- double_fit[[1]]
          sims_1_1 <- chain_1[keep_1, c(11, 9, 5, 7)]
          sims_1_2 <- chain_1[keep_1, c(12, 10, 6, 8)]

          chain_2 <- double_fit[[2]]
          sims_2_1 <- chain_2[keep_2, c(11, 9, 5, 7)]
          sims_2_2 <- chain_2[keep_2, c(12, 10, 6, 8)]
          
          chain_3 <- double_fit[[3]]
          sims_3_1 <- chain_3[keep_3, c(11, 9, 5, 7)]
          sims_3_2 <- chain_3[keep_3, c(12, 10, 6, 8)]

          sims_1 <- rbind(sims_1_1, sims_2_1, sims_3_1)
          sims_2 <- rbind(sims_1_2, sims_2_2, sims_3_2)

           steps_1 <- numeric(length = nrow(sims_1))
           turns_1 <- numeric(length = nrow(sims_1))
           steps_2 <- numeric(length = nrow(sims_2))
           turns_2 <- numeric(length = nrow(sims_2))
      
          for(i in 1:nrow(sims_1)){
            steps_1[i] <- rweibull(1, sims_1[i,1], sims_1[i,2])
            turns_1[i] <- rwrappedcauchy(1, sims_1[i,3],  sims_1[i,4])
            }
            
          for(i in 1:nrow(sims_2)){
            steps_2[i] <- rweibull(1, sims_2[i,1], sims_2[i,2])
            turns_2[i] <- rwrappedcauchy(1, sims_2[i,3],  sims_1[i,4])
            }

        post_iters_1 <- as.data.frame(cbind(steps_1, turns_1))                     
        post_iters_2 <- as.data.frame(cbind(steps_2, turns_2))                
                      
steps_plot_1 <- ggplot(post_iters_1) + 
                          geom_histogram(aes(steps_1), binwidth = 0.1, fill = "light grey", colour = "grey") +
                          xlab("Step length (km)") +
                          ylab("Frequency") +
                          ylim(c(0, 625)) +
                          theme_bw()
steps_plot_2 <- ggplot(post_iters_2) + 
                          geom_histogram(aes(steps_2), binwidth = 0.1, fill = "light grey", colour = "grey") +
                          xlab("Step length (km)") +
                          ylab("Frequency") +
                          ylim(c(0, 625)) +
                          theme_bw()

turns_plot_1 <- ggplot(post_iters_1) + 
                          geom_histogram(aes(turns_1),binwidth = 0.1, fill = "light grey", colour = "grey") +
                          xlab("Turn angle (rad)") +
                          ylab("Frequency") +
                          xlim(c(-3.25, 3.25)) +
                          theme_bw()                     

                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
# #  Parameter estimates from "XXX" model ouput
# #  Check out just the steps and turning angles.
# #   rweibull(n, shape, scale = 1)
# #   rwrappedcauchy(n, mu = circular(0), rho = exp(-1), control.circular=list())

# library(RColorBrewer)
# myColors <- gray.colors(2, start = 0.6, end = 0.8)
# names(myColors) <- levels(steps_df$variable)
# colScale <- scale_colour_manual(name = "variable",values = myColors)

# #  Randomly select a values for each paramter given estimate mean and 95% CI 
# #   and to generate a density plot from 1000 values of a step length and turning angle from these 
# #   selected values.

 
# #   Make steps from state 1
# steps_1 <-rweibull(100, 0.93, 0.425)
# #   Make steps from state 2
# steps_2 <-rweibull(100, 1.65,  0.4034)
# #   Dataframe for steps
# steps_df_tmp <- as.data.frame(cbind(steps_1, steps_2))
# names(steps_df_tmp)[1] <- "Stationary"
# names(steps_df_tmp)[2] <- "Transit"
# steps_df <- melt(steps_df_tmp)
# #   Figure of step lengths distribution from state 1 and state 2
# compare_steps <- ggplot(steps_df) + 
                                # geom_density(aes(x = value, colour = variable, linetype = variable), size = 1, alpha = .25) +
                                # xlab("Step length (km)") +
                                # theme_bw()
# compare_steps + colScale


# #   Make turning angles from state 1
# theta_1 <- rwrappedcauchy(100, 0.1737,  335)
# theta_1[theta_1>pi]=theta_1[theta_1>pi]-2*pi
# #   Make turning angles from state 2
# theta_2 <- rwrappedcauchy(100,  -0.1691, 0.2766)
# theta_2[theta_2>pi]=theta_2[theta_2>pi]-2*pi
# #   Dataframe for turning anlges
# angles_df_tmp <- as.data.frame(cbind(theta_1, theta_2))
# names(angles_df_tmp)[1] <- "Stationary"
# names(angles_df_tmp)[2] <- "Transit"
# angles_df <- melt(angles_df_tmp)
# #   Figure of turning angle distribution from state 1 and state 2
# compare_angles <- ggplot(angles_df) + 
                                 # geom_density(aes(x = value, colour = variable, linetype = variable), size = 1, alpha = .25) +
                                 # xlab("Turn angle (radian)") +
                                 # theme_bw()
# compare_angles + colScale








#  Generate many simulations of N steps of movement from state 1 and state 2
#   Number of simulations
sims <- 10
#   Length of walk
N <- 3
#   Matrix to hold state 1 simulations
mat_X_1 <- matrix(nrow=N, ncol = sims)
mat_Y_1 <- matrix(nrow=N, ncol = sims)
for(j in 1: sims){
   # make distributed steps
   steps_sim_1 <-rweibull(N, 0.912288, 394.524400)
   # make clustered turning angles
   theta_sim_1 <- rwrappedcauchy(N, -0.011679,  0.268873)
   theta_sim_1[theta_sim_1>pi]=theta_sim_1[theta_sim_1>pi]-2*pi
   # cumulative angle (absolute orientation)
   phi_sim_1 <- cumsum(theta_sim_1)
   # step length components
   dX_sim_1 <- steps_sim_1*cos(phi_sim_1)
   dY_sim_1 <- steps_sim_1*sin(phi_sim_1)
   # actual X-Y values
   X_sim_1 <- as.matrix(cumsum(dX_sim_1))
   Y_sim_1 <- as.matrix(cumsum(dY_sim_1))
   mat_X_1[,j] <- X_sim_1
   mat_Y_1[,j] <- Y_sim_1
}
#   Matrix to hold state 2 simulations
mat_X_2 <- matrix(nrow=N, ncol = sims)
mat_Y_2 <- matrix(nrow=N, ncol = sims)
for(j in 1: sims){
   # make distributed steps
   steps_sim_2 <-rweibull(N, 1.060, 652.5)
   # make clustered turning angles
   theta_sim_2 <- rwrappedcauchy(N,   0.1217,  0.2485)
   theta_sim_2[theta_sim_2>pi]=theta_sim_2[theta_sim_2>pi]-2*pi
   # cumulative angle (absolute orientation)
   phi_sim_2<- cumsum(theta_sim_2)
   # step length components
   dX_sim_2 <- steps_sim_2*cos(phi_sim_2)
   dY_sim_2 <- steps_sim_2*sin(phi_sim_2)
   # actual X-Y values
   X_sim_2 <- as.matrix(cumsum(dX_sim_2))
   Y_sim_2 <- as.matrix(cumsum(dY_sim_2))
   mat_X_2[,j] <- X_sim_2
   mat_Y_2[,j] <- Y_sim_2
}
#   Data from for locations in state 1 movement mode.
df_X_1_tmp <- as.data.frame(mat_X_1)
df_Y_1_tmp <- as.data.frame(mat_Y_1)
df_X_1 <- melt(df_X_1_tmp)
df_Y_1 <- melt(df_Y_1_tmp)
df_XY_1_tmp <- cbind(df_X_1 , df_Y_1)
names(df_XY_1_tmp) <- c("walk_num", "X", "walk_num_rep", "Y")
df_XY_1 <- df_XY_1_tmp %>%
               dplyr::select(walk_num, X, Y) %>%
               mutate(mov_mod = "State 1")
#   Data from for locations in state 2 movement mode.
df_X_2_tmp <- as.data.frame(mat_X_2)
df_Y_2_tmp <- as.data.frame(mat_Y_2)
df_X_2 <- melt(df_X_2_tmp)
df_Y_2 <- melt(df_Y_2_tmp)
df_XY_2_tmp <- cbind(df_X_2 , df_Y_2)
names(df_XY_2_tmp) <- c("walk_num", "X", "walk_num_rep", "Y")
df_XY_2 <- df_XY_2_tmp %>%
                   dplyr::select(walk_num, X, Y) %>%
                   mutate(mov_mod = "State 2")
#  Datafarme of all locations together.
sim_crw_df <- as.data.frame(bind_rows(df_XY_1  , df_XY_2))
#   Figure of all simulated locations/paths from state 1 and state 2 movement modes.
compare_sim_crw <- ggplot(sim_crw_df) + 
                                     geom_line(aes(x = X, y = Y, 
                                                               group = walk_num, 
                                                               colour = mov_mod), alpha = 0.4) +
                                     xlab("X (km)") +
                                     ylab("Y (km)") +
                                     theme_bw()
compare_sim_crw 



steps_1_sim <-as.data.frame(rweibull(100, 0.949, 0.43))
names(steps_1_sim)[1] <- "step_length"

steps_1_obs <- as.data.frame((l_sing[l_sing < 1500])/1000)
names(steps_1_obs)[1] <- "step_length"

steps_1_sim_plot <- ggplot(steps_1_sim) + 
                          geom_density(aes(step_length)) +
                          xlab("Step length (m)") +
                          theme_minimal() +
                          theme(axis.title.x = element_blank())
                          
steps_1_obs_plot <- ggplot(steps_1_obs) +
                                  geom_histogram(aes(step_length)) +
                                  theme_minimal() +
                                  theme(axis.title.x = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(steps_1_sim_plot), ggplotGrob(steps_1_obs_plot), size = "last"))


steps_2_sim <-as.data.frame(rweibull(100, 1.673, 333))
names(steps_2_sim)[1] <- "step_length"

steps_2_obs <- as.data.frame((l_sing[l_sing > 1500])/1000)
names(steps_2_obs)[1] <- "step_length"

steps_2_sim_plot <- ggplot(steps_2_sim) + 
                          geom_density(aes(step_length)) +
                          xlab("Step length (m)") +
                          theme_minimal() +
                          theme(axis.title.x = element_blank())
                          
steps_2_obs_plot <- ggplot(steps_2_obs) +
                                  geom_histogram(aes(step_length)) +
                                  theme_minimal() +
                                  theme(axis.title.x = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(steps_2_sim_plot), ggplotGrob(steps_2_obs_plot), size = "last"))



theta_1 <- rwrappedcauchy(1000, -0.011679,  0.268873)
theta_1[theta_1>pi]=theta_1[theta_1>pi]-2*pi


