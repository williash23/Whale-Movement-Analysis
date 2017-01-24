# Sara Williams
# 7/26/2016; 
# Generate plots from MCMC posterior distribution object (iterations)
#   and run in JAGS.
################################################################################

#  Load packages
library(plyr)
library(packcircles)
library(ggmcmc)
library(gridExtra)
library(ggthemes)
################################################################################



#   Number of movement paths
npath <- 10000
#   Length of walk
nobs <- 3

mod <- post_sims_su_plot
#   Matrices to hold simulations
mat_X <- matrix(nrow = nobs, ncol = npath)
mat_Y <- matrix(nrow = nobs, ncol = npath)

for(i in 1:npath){
     keep <- sample(1:nrow(mod), nobs, replace = FALSE)
 
       # make distributed steps
       steps_sim <- mod[keep, 1]
       # make clustered turning angles
       theta_sim <- mod[keep, 2]
       theta_sim[1] <- 0
       #theta_sim[2] <- 0
       #theta_sim[2] <- circ[i]
       # cumulative angle (absolute orientation)
       phi_sim <- cumsum(theta_sim)
       # step length components -- 
       dX_sim <- steps_sim*cos(phi_sim)
       dX_sim[1] <- 0 ##MAKE EACH START AT [0,0]
       #dX_sim[2] <- 0
       dY_sim <- steps_sim*sin(phi_sim)
       dY_sim[1] <- 0
       #dY_sim[2] <- 0.5
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

#  Make data suitable for plotting as points
df_XY_pt <- df_XY %>%
                    filter(location_num == 2) 





mod_2 <- post_sims_d_plot
#   Matrices to hold simulations
mat_X_2 <- matrix(nrow = nobs, ncol = npath)
mat_Y_2 <- matrix(nrow = nobs, ncol = npath)

for(i in 1:npath){
     keep <- sample(1:nrow(mod_2), nobs, replace = FALSE)
 
       # make distributed steps
       steps_sim_2 <- mod_2[keep, 1]
       # make clustered turning angles
       theta_sim_2 <- mod_2[keep, 2]
       theta_sim_2[1] <- 0
       #theta_sim[2] <- 0
       #theta_sim[2] <- circ[i]
       # cumulative angle (absolute orientation)
       phi_sim_2 <- cumsum(theta_sim_2)
       # step length components -- 
       dX_sim_2 <- steps_sim_2*cos(phi_sim_2)
       dX_sim_2[1] <- 0 ##MAKE EACH START AT [0,0]
       #dX_sim[2] <- 0
       dY_sim_2 <- steps_sim_2*sin(phi_sim_2)
       dY_sim_2[1] <- 0
       #dY_sim[2] <- 0.5
       # actual X-Y values
       X_sim_2 <- as.matrix(cumsum(dX_sim_2))
       Y_sim_2 <- as.matrix(cumsum(dY_sim_2))
       mat_X_2[,i] <- X_sim_2
       mat_Y_2[,i] <- Y_sim_2
   }

#  Data from for locations.
df_X_tmp_2 <- as.data.frame(mat_X_2)
df_Y_tmp_2 <- as.data.frame(mat_Y_2)
df_X_2 <- melt(df_X_tmp_2)
df_Y_2 <- melt(df_Y_tmp_2)
df_XY_tmp_2 <- cbind(df_X_2 , df_Y_2, rep(1:nobs))
names(df_XY_tmp_2) <- c("walk_num", "X", "walk_num_rep", "Y", "location_num")
df_XY_2 <- df_XY_tmp_2 %>%
                    dplyr::select(walk_num, location_num, X, Y) %>%
                    mutate(model = "Single State")

#  Make data suitable for plotting as points
df_XY_pt_2 <- df_XY_2 %>%
                         filter(location_num == 2) 




whale_loc <- data.frame(x = 0, y = 0)
one_move_sim_pt_z <- ggplot() + 
                                        geom_point(data = df_XY_pt_2, alpha= 0.1, aes(x = X, y = Y, group = walk_num), colour = "dark blue", fill = "dark blue", size = 2) +
                                        xlab("X (km)") +
                                        xlim(c(-3, 3)) +
                                        ylab("Y (km)") +
                                        ylim(c(3, -3)) +
                                        theme_bw() +
                                        theme(panel.grid.minor = element_blank(), axis.text = element_text(size = rel(2)), 
                                        axis.title = element_text(size = rel(2))) +
                                        geom_point(data = df_XY_pt, alpha= 0.1, aes(x = X, y = Y, group = walk_num), colour = "dark red", fill = "dark red", size = 2) +
                                        geom_point(data = whale_loc, aes(x, y), size=3, colour="black") +
                                        #geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), colour = "black", size = 1, 
                                                                          #data = whale_seg, arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
                                        annotate("point", x = 0, y = 0, colour ="black", size = 280, shape=1) 
                                        #geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), colour = "black", size = 1.5, 
                                                                          #data = ship_line, arrow = arrow(length = unit(0.15, "inches"), type = "closed")) +
                                        #geom_text(data = ship_loc1, aes(x, y, label=label), family="OpenSansEmoji", size=12, colour = "#dc143c") +
                                        #geom_text(data = ship_loc2, aes(x, y, label=label), family="OpenSansEmoji", size=12, colour = "#dc143c", alpha = 0.2) 
                                        
one_move_sim_pt_z
































ship_line <- data.frame(x1 = -2, x2 = -2, y1 = 1.8, y2 = 0.5)
ship_loc1 <- data.frame(x = -2, y = 1.8, label = (emoji('ship')))
ship_loc2 <- data.frame(x = -2, y = 0.2, label = (emoji('ship')))
#whale_loc <- data.frame(x = 0.2, y = 0, label = (emoji('whale')))

whale_seg <- data.frame(x1= 6.25, y1=0, x2 = 7, y2 = 0)


#  Figure of all simulated location (points) and ship
one_move_sim_pt <- ggplot() + 
                                    geom_point(data = df_XY_pt, alpha= 0.1, aes(x = X, y = Y, group = walk_num), colour = "grey50", fill = "grey50", size = 2) +
                                    xlab("X (km)") +
                                    xlim(c(-10, 10)) +
                                    ylab("Y (km)") +
                                    ylim(c(10, -10)) +
                                    theme_bw() +
                                    theme(panel.grid.minor = element_blank(), axis.text = element_text(size = rel(2)), 
                                                axis.title = element_text(size = rel(2))) +
                                    geom_point(data = whale_loc, aes(x, y), size=3, colour="black") +
                                    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), colour = "black", size = 1, 
                                                                      data = whale_seg, arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
                                    annotate("point", x = 0, y = 0, colour ="black", size = 280, shape=1) 
                                    #annotate("text", x = -7.75, y = 0, label = "pi", parse = T, size = 10) +
                                    #annotate("text", x = 7.75, y = 0, label = "0", parse = T, size = 10) +
                                    #annotate("text", x = 0, y = 7.75, label = "pi/2", parse = T, size = 10) +
                                    #annotate("text", x = 0, y = -7.75, label = "3 pi/2", parse = T, size = 10) 
                                    #geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), colour = "black", size = 1.5, 
                                                                      #data = ship_line, arrow = arrow(length = unit(0.15, "inches"), type = "closed")) +
                                    #geom_text(data = ship_loc1, aes(x, y, label=label), family="OpenSansEmoji", size=12, colour = "#dc143c") +
                                    #geom_text(data = ship_loc2, aes(x, y, label=label), family="OpenSansEmoji", size=12, colour = "#dc143c", alpha = 0.2) 
                                    
one_move_sim_pt

#  Same as figure above but zoomed in
whale_seg <- data.frame(x1=1.22, y1=0, x2 = 1.35, y2 = 0)
















#  Make circle plot - includes all turn angles but uses pulls from posterior for steps
#   radiu = poster pull step length
#   Instead of using entire circles using pie pieces???
#  Generate vectors to hold data with package 'packcircles'
r <- post_sims_plot$steps
id <- seq(1, length(r), 1)
x <- rep(0, length(r)) 
y <- rep(0, length(r)) 
points <- data.frame(x, y, r, id)
points <- dplyr::arrange(points, by = r)

#  Generate circle vertices for each iterations step length
#   Here, each iteration's step length is the radius of the circle
dat <- circlePlotData(points, npoints = 100, xyr.cols = 1:3, id.col = 4)

#  Plot all of these circles (45000 circles from the 45000 simulated values from 45000 MCMC iterations)
p <- ggplot() +
        geom_polygon(data=dat, aes(x, y, group=id), fill = "grey47", alpha = 0.1, linetype = 0) +
        coord_equal() +
        theme_bw() +
        xlab("X (km)") +
        ylab("Y (km)") +
        theme(panel.grid.minor = element_blank(), axis.text = element_text(size = rel(2)), 
        axis.title = element_text(size = rel(2)))
p





color <-graycol(n = 50)
dr <- 0.05
emptyplot(xlim = c(-2, 2), col = color[length(color)],
main = "filledcircle")
filledcircle(r1 = 1, mid = c(1, 1), dr = dr,
col = shadepalette(endcol = "darkblue"))

filledcircle(r1 = post_sims$steps[1], r2 = 0, mid = c(0,0), dr = 0.01, from = 0, to = post_sims$turns[1], col = shadepalette(endcol = "darkblue"))




 draw.arc(x=0,y=0,radius=post_sims$steps[1],angle1=0,angle2=post_sims$turns[1],n=0.05)