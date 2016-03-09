# Sara Williams
# 2/22/2016; 
# Generate figures from parameter estimates from models formulated after code from Morales et al. 2004
#   and run in JAGS.
################################################################################

#  Load packages
library(circular)
library(reshape)
library(ggplot2)
library(dplyr)

################################################################################
#  Parameter estimates from "double switch" model ouput
#  Check out just the steps and turning angles.
#   rweibull(n, shape, scale = 1)
#   rwrappedcauchy(n, mu = circular(0), rho = exp(-1), control.circular=list())

#   Make steps from state 1
steps_1 <-rweibull(5000, 1.287852, 0.154374)
#   Make turning angles from state 1
theta_1 <- rwrappedcauchy(5000, -0.72212,  0.07696)
theta_1[theta_1>pi]=theta_1[theta_1>pi]-2*pi
#   Make steps from state 2
steps_2 <-rweibull(5000,  1.117634, 0.53264)
#   Make turning angles from state 2
theta_2 <- rwrappedcauchy(5000, -0.19032,  0.22263)
theta_2[theta_2>pi]=theta_2[theta_2>pi]-2*pi
#   Dataframe for steps
steps_df_tmp <- as.data.frame(cbind(steps_1, steps_2))
names(steps_df_tmp)[1] <- "Stationary"
names(steps_df_tmp)[2] <- "Transit"
steps_df <- melt(steps_df_tmp)
#   Figure of step lengths distribution from state 1 and state 2
compare_steps <- ggplot(steps_df) + 
                   geom_density(aes(x = value, colour = variable, fill = variable), alpha = .25) +
                     xlab("Step length (km)") +
                     theme_bw()
compare_steps
#   Dataframe for turning anlges
angles_df_tmp <- as.data.frame(cbind(theta_1, theta_2))
names(angles_df_tmp)[1] <- "Stationary"
names(angles_df_tmp)[2] <- "Transit"
angles_df <- melt(angles_df_tmp)
#   Figure of turning angle distribution from state 1 and state 2
compare_angles <- ggplot(angles_df) + 
                      geom_density(aes(x = value, colour = variable, fill = variable), alpha = .25) +
                      xlab("Turn angle (radian)") +
                      theme_bw()
compare_angles

#  Generate many simulations of N steps of movement from state 1 and state 2
#   Number of simulations
sims <- 1000
#   Length of walk
N <- 3
#   Matrix to hold state 1 simulations
mat_X_1 <- matrix(nrow=N, ncol = sims)
mat_Y_1 <- matrix(nrow=N, ncol = sims)
for(j in 1: sims){
   # make distributed steps
   steps_sim_1 <-rweibull(N, 1.287852, 0.154374)
   # make clustered turning angles
   theta_sim_1 <- rwrappedcauchy(N, -0.72212,  0.07696)
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
   steps_sim_2 <-rweibull(N, 1.117634, 0.53264)
   # make clustered turning angles
   theta_sim_2 <- rwrappedcauchy(N,  -0.19032,  0.22263)
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