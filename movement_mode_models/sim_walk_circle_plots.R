# Sara Williams
# 7/26/2016; 
# Generate plots from MCMC posterior distribution object (iterations)
#   and run in JAGS.
################################################################################

#  Load packages
library(dplyr)
library(packcircles)
library(ggmcmc)
library(gridExtra)
library(ggthemes)
################################################################################



#   Number of movement paths

#   Length of walk
nobs <- 1

mod <- post_sims_si_plot
npath <- 1000

loc_3_df <- matrix(nrow = npath, ncol =2)

for(i in 1:npath){

     keep <- sample(1:nrow(mod), nobs, replace = FALSE)

	start_x <- 433525.87
	start_y <- 6507454.23
	sec_x <-  433525.87
	sec_y <- start_y + 160

	step_l <- mod[keep,1] * 1000
	turn_a <- mod[keep, 2] 
	new_x <- sec_x + (step_l * cos(turn_a))
	new_y <- sec_y + (step_l * sin(turn_a))

	loc_3_df[i,1] <- new_x
	loc_3_df[i,2] <- new_y
	}
	
loc_1 <- as.data.frame(cbind(rep(start_x, npath), rep(start_y, npath), 
	rep(1, npath), as.data.frame(seq(1:npath))))
loc_2 <- as.data.frame(cbind(rep(sec_x, npath), rep(sec_y, npath), 
	rep(2, npath), as.data.frame(seq(1:npath))))
loc_3 <- cbind(loc_3_df, as.data.frame(rep(3, npath)), as.data.frame(seq(1:npath)))
	
names(loc_1)[1:4] <- c("X", "Y", "loc_num", "path_num")
names(loc_2)[1:4] <- c("X", "Y", "loc_num", "path_num")
names(loc_3)[1:4] <- c("X", "Y", "loc_num", "path_num")	
	
	
loc_df <- rbind(loc_1, loc_2, loc_3) %>%
	dplyr::arrange(path_num, loc_num)
loc_df2 <- loc_df %>%dplyr::filter(loc_num == 3)


		# p <- ggplot() + 
			# #geom_line(data = loc_df2, aes(x = X, Y = y, group = path_num), 
					# #colour = "grey60", size = 0.5, alpha= 0.1) +
			# geom_point(data = loc_df2, aes(x = X, y = Y), colour = "grey60", size = 1, alpha= 0.1) +
			# #ylim(6507000, 6508400) +
			# theme_bw() +
			# ylab("Longitude (UTM)/n") +
			# xlab("Latitude (UTM)")
		# p

pts <- loc_df2
coordinates(pts) <- ~X+Y
proj4string(pts) <- CRS("+proj=utm +zone=8 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

r_mat <- matrix(0, 100, 100)
r_template <- raster(r_mat)
extent(r_template) <- extent(pts)
projection(r_template) <- CRS("+proj=utm +zone=8 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

r <- r_template
tab <- table(cellFromXY(r, pts))
r[as.numeric(names(tab))] <- tab

	range01 <- function(x){(x-x_min)/(x_max - x_min)}
	x_min <- r@data@min
	x_max <- r@data@max
	r2 <- raster::calc(r, range01)
	#r2[r2 == 0] <- NA
	
	
r_mat <- matrix(0, 50, 50)
r_template <- raster(r_mat)
extent(r_template) <- extent(pts)
projection(r_template) <- CRS("+proj=utm +zone=8 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

	stack <- stack(tmp)
	for(p in 2:1000){
		loc_2_df <- matrix(nrow = npath, ncol =2)
		loc_3_df <- matrix(nrow = npath, ncol =2)

		for(i in 1:npath){

			 keep <- sample(1:nrow(mod), nobs, replace = FALSE)

			start_x <- 433525.87
			start_y <- 6507454.23
			#sec_x <-  433525.87
			#sec_y <- start_y + 160

			step_l <- mod[keep,1] * 1000
			turn_a <- mod[keep, 2]
			new_x <- start_x + (step_l * cos(turn_a))
			new_y <- start_y + (step_l * sin(turn_a))
			
			step_l <- mod[keep,1] * 1000
			turn_a <- mod[keep, 2]
			new2_x <- new_x + (step_l * cos(turn_a))
			new2_y <- new_y + (step_l * sin(turn_a))
			
			loc_2_df[i,1] <- new_x
			loc_2_df[i,2] <- new_y
			
			loc_3_df[i,1] <- new2_x
			loc_3_df[i,2] <- new2_y
			}
			
		loc_1 <- as.data.frame(cbind(rep(start_x, npath), rep(start_y, npath), 
			rep(1, npath), as.data.frame(seq(1:npath))))
		#loc_2 <- as.data.frame(cbind(rep(sec_x, npath), rep(sec_y, npath), 
			#rep(2, npath), as.data.frame(seq(1:npath))))
		
		loc_2 <- cbind(loc_2_df, as.data.frame(rep(2, npath)), as.data.frame(seq(1:npath)))
		loc_3 <- cbind(loc_3_df, as.data.frame(rep(3, npath)), as.data.frame(seq(1:npath)))
		
		names(loc_1)[1:4] <- c("X", "Y", "loc_num", "path_num")
		names(loc_2)[1:4] <- c("X", "Y", "loc_num", "path_num")
		names(loc_3)[1:4] <- c("X", "Y", "loc_num", "path_num")	
			
			
		loc_df <- rbind(loc_1, loc_2, loc_3) %>%
			dplyr::arrange(path_num, loc_num)
		loc_df2 <- loc_df %>%dplyr::filter(loc_num == 3)

		pts <- loc_df2
		coordinates(pts) <- ~X+Y
		proj4string(pts) <- CRS("+proj=utm +zone=8 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

	r <- r_template
	tab <- table(cellFromXY(r, pts))
	r[as.numeric(names(tab))] <- tab

	x_min <- r@data@min
	x_max <- r@data@max
	tmp <- raster::calc(r, range01)

	stack[[p]] <- tmp
	}
	
	s <- sum(stack)
	x_min <- s@data@min
	x_max <- s@data@max
	s_p <- raster::calc(s, range01)
	
cuts <- seq(0,1,0.005)
pal <- colorRampPalette(c("white","red"))

par(yaxs="i",las=0)
plot(s_p, ylab = "Latitude (UTM)", xlab = "Longitude (UTM)", breaks=cuts, col = pal(length(cuts)),
	xlim = c(433000, 434000), ylim = c(6507000, 6508000))
box(bty="l")
#segments(start_x, start_y, sec_x, sec_y, cex = 50)
#points(start_x, start_y, col = "black", pch=19, cex = 0.9)
points(sec_x, sec_y, col = "black", pch=19, cex = 0.9)
grid(nx=NULL,ny=NULL,lty=1,lwd=1,col="gray")	
	
	
	
	
	
	
	start_x <- 433525.87
	start_y <- 6507454.23
	#sec_x <-  433525.87
	#sec_y <- start_y + 160


	z <- simm.crw(date=1:2, 
			h = 2412, 
			r = 0.203,
			x0 = c(start_x, start_y), 
			burst= "1",
			typeII=FALSE,
			proj4string=CRS("+proj=utm +zone=8 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))	

	for(i in 2:10000){
		
		ids <- as.character(i +1)
		tmp <- simm.crw(date=1:2, 
			h = 2412, 
			r = 0.203,
			x0 = c(start_x, start_y), 
			burst= ids, 
			typeII=FALSE, 
			proj4string=CRS("+proj=utm +zone=8 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))	

		z <- c(z, tmp)
		}
		
		
		
	#plot(z, pch = 19)
		
		
		
		paths <- paths_tmp %>%
			mutate(step_num = rep(seq(1, 2, 1), length(z))) %>%
			mutate(path_num = rep(1:2, each = length(z))) %>%
			dplyr::select(x, y, step_num, path_num) %>%
			filter(step_num == 2)
		
		pts <- paths
		coordinates(pts) <- ~x+y
		proj4string(pts) <- CRS("+proj=utm +zone=8 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

		
	r <- r_template
	tab <- table(cellFromXY(r, pts))
	r[as.numeric(names(tab))] <- tab

	x_min <- r@data@min
	x_max <- r@data@max
	r2<- raster::calc(r, range01)

cuts <- seq(0,1,0.05)
pal <- colorRampPalette(c("white","red"))

par(yaxs="i",las=0)
plot(r2, ylab = "Latitude (UTM)", xlab = "Longitude (UTM)", breaks=cuts, col = pal(length(cuts))) 
	#xlim = c(433000, 434000), ylim = c(6507000, 6508000))
box(bty="l")
#segments(start_x, start_y, sec_x, sec_y, cex = 50)
#points(start_x, start_y, col = "black", pch=19, cex = 0.9)
points(sec_x, sec_y, col = "black", pch=19, cex = 0.9)
grid(nx=NULL,ny=NULL,lty=1,lwd=1,col="gray")	
	
	
		paths <- paths_tmp %>%
			mutate(step_num = rep(seq(1, 2, 1), length(z))) %>%
			mutate(path_num = rep(1:2, each = length(z))) %>%
			dplyr::select(x, y, step_num, path_num) 
	
		p <- ggplot() + 
			geom_point(data = paths, aes(x = x, y = y), 
					colour = "darkred", size = 2, alpha= 0.15) +
			#geom_point(data = loc_df2, aes(x = X, y = Y), colour = "grey60", size = 1, alpha= 0.1) +
			#ylim(6507000, 6508400) +
			theme_bw() +
			ylab("Longitude (UTM)") +
			xlab("Latitude (UTM)")
		p
		

	
	
	
	
	
	
	
	
for(i in 1:npath){
     keep <- sample(1:nrow(mod), nobs, replace = FALSE)
 
       # make distributed steps
       steps_sim <- mod[keep, 1]
       # make clustered turning angles
       theta_sim <- mod[keep, 4]
       theta_sim[1] <- 0
       theta_sim[2] <- 0
       theta_sim[3] <- theta_sim[3]
       # cumulative angle (absolute orientation)
       phi_sim <- cumsum(theta_sim)
       # step length components -- 
       dX_sim <- steps_sim*cos(phi_sim)
       dX_sim[1] <- 0 ##MAKE EACH START AT [0,0]
       dX_sim[2] <- 0.4
       dY_sim <- steps_sim*sin(phi_sim)
       dY_sim[1] <- 0
       dY_sim[2] <- 0
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
                    filter(location_num > 1 ) 
					#mutate(new_X = ifelse(location_num == 2, 0.4, X))
p <- ggplot() + 
	geom_line(data = df_XY_pt, alpha= 0.1, aes(x = X, y = Y, group = walk_num), 
	#geom_point(data = df_XY_pt, alpha= 0.1, aes(x = X, y = Y), 
		colour = "grey60", size = 1) +
	theme_bw()
p



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
       # theta_sim_2[1] <- 0
       #theta_sim[2] <- 0
       #theta_sim[2] <- circ[i]
       # cumulative angle (absolute orientation)
       phi_sim_2 <- theta_sim_2
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
	geom_point(data = df_XY_pt, alpha= 0.1, aes(x = X, y = Y, group = walk_num), colour = "dark blue", fill = "dark blue", size = 2) +
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