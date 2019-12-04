# Sara Williams
# 3/18/2017
# Whale ship strike risk simulation - function to run simulaiton for 13 knot ship speed
#   Have to load script: sim_once_to_start.R before running this one
################################################################################

library(raster)
library(rgeos)
library(sp)
library(sf)
library(dplyr)
library(tidyr)
library(rasterVis)

UTM <- CRS("+proj=utm +zone=8 +datum=WGS84")


#  Load data
# Load output from MCMC posterior distribtuion iterations (JAGS object)
load("C:/Users/saraw/OneDrive/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit_surf.RData")

rwcauchy <- function(n, mu = 0, rho = 0) {
  u = runif(n)
  V = cos(2 * pi * u)
  c = 2 * rho/(1 + rho^2)
  t <- (sign(runif(n) - 0.5) * acos((V + c)/(1 + c * V)) + mu)%%(2 * pi)
  return(t)
}


#  Generate simulated step lengths and turn angles for whale movements from parameter 
#   estimate posterior distritbutions from movement analysis
#   Currently using single movement mode with all data model

tmp <- read.csv("C:/Users/saraw/Documents/Whales/data/Whales_0615_general_clean.csv") %>%
	group_by(same_whale_ID) %>%
	dplyr::filter(count < 2) %>%
	dplyr::filter(n() > 2) %>%
	dplyr::select(unique_event_ID, same_whale_ID, TimeTxt, X_whale_UTM, Y_whale_UTM) %>%
	as.data.frame() %>%
	separate(TimeTxt, c("hr", "mn", "sc")) 
tmp$hr <- as.numeric(tmp$hr)
tmp$mn <- as.numeric(tmp$mn)
tmp$sc <- as.numeric(tmp$sc)

mu_x <- mean(tmp$X_whale_UTM, na.rm = TRUE)
mu_y <- mean(tmp$Y_whale_UTM, na.rm = TRUE)
 
tmp2 <- tmp %>%
	dplyr::mutate(tim_sec = (hr*60*60) + (mn*60) + sc) %>%
	group_by(same_whale_ID) %>%
	mutate(tim_dif_sec = (lag(tim_sec) - tim_sec) *-1) %>%
	as.data.frame() 
tmp2$tim_dif_sec[tmp2$tim_dif_sec > 1000] <- NA
tmp2$tim_dif_sec[tmp2$tim_dif_sec < -1000] <- NA
tmp3 <- tmp2 %>% dplyr::filter(tim_dif_sec < 50 & tim_dif_sec > 0)
# mean swim speed  = 2.5182
# sd = 10.2

 

	#  Mean parameter estimates for "short" dive length
	mu <- 0.07 
	rho <-0.21
	v <- 1.47
	lambda <- 34.03
 
  niter <- 10000
  nsamp <- 2
  
  x <- single_fit_surf

	
	#  Background set up
	#  13 knots = 6.68778 m/s
	#  19 knots = 9.77444 m/s
	ship_spd <- 9.77444
	whale_spd <- 2.518
	btw_surf <- 30
	mov_time <- 180
	ship_vert_dist <- mov_time * ship_spd
	
	perp_dist_ship_whale <- whale_spd*mov_time
	vert_dist_ship_whale <- ship_vert_dist
	
	# Whale starting location
	init_x <- 0  #mu_x
	init_y <- 0 #mu_y
	
	# Number of whale surfacings
	n_obs <- 6 #roundup(mov_time/btw_surf)
	
	# Ship starting location
	init_ship_x <- init_x - perp_dist_ship_whale
	init_ship_y <- init_y - vert_dist_ship_whale
	
	# Ship movement
	num_ship_pts <- round(vert_dist_ship_whale/(ship_spd*btw_surf))
	
	ship_init_loc_df <- as.data.frame(cbind(init_ship_x, init_ship_y))
	colnames(ship_init_loc_df) <- c("x", "y")
	
	ship_loc1 <- c(init_ship_x, init_ship_y)
	ship_loc2 <- c(init_ship_x, ship_loc1[2] + (ship_spd * 30))
	ship_loc3 <- c(init_ship_x, ship_loc2[2] + (ship_spd * 30))
	ship_loc4 <- c(init_ship_x, ship_loc3[2] + (ship_spd * 30))
	ship_loc5 <- c(init_ship_x, ship_loc4[2] + (ship_spd * 30))
	ship_loc6 <- c(init_ship_x, ship_loc5[2] + (ship_spd * 30))
	ship_loc7 <- c(init_ship_x, ship_loc6[2] + (ship_spd * 30))
	ship_loc_df <- as.data.frame(rbind(ship_loc1, ship_loc2, ship_loc3, ship_loc4, 
		ship_loc5, ship_loc6, ship_loc7))
	colnames(ship_loc_df) <- c("x", "y")
	coordinates(ship_loc_df) <- ~x+y
	ship_loc_sf <- st_as_sf(ship_loc_df)  %>%
		st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
		mutate(obs_num = 1:n())
	ship_loc_sp <- as(ship_loc_sf, "Spatial")
	
	ship_ln <- st_as_sf(ship_loc_df) %>%
		dplyr::summarise(do_union = FALSE) %>%
		st_cast("LINESTRING") %>%
		st_set_crs("+proj=utm +zone=8 +datum=WGS84")
	
	ship_strk_area <- st_as_sf(ship_loc_sp[7,])
    #st_as_sf(ship_loc_sp[6:7,]) %>%
		#dplyr::summarise(do_union = FALSE) %>%
		#st_cast("LINESTRING") 
		
	ship_strk_area_b <- st_buffer(ship_strk_area, 50)
	ship_ln_b <- st_buffer(ship_ln, 10)

	
n_strk_sim <- 10000
strk_sim <- as.vector(rep(NA, length = n_strk_sim))
path_sim <- list()
#
#tst  <- adehabitatLT::simm.crw(date=1:n_obs, h = lambda, r = rho,
#  x0=c(init_x, init_y), id="A1", burst=1,
#  typeII=FALSE, proj4string=CRS("+proj=utm +zone=8 +datum=WGS84"))

for(i in 1:n_strk_sim){
	
  keep_1 <- sample(1:niter, nsamp, replace = F)
  keep_2 <- sample(1:niter, nsamp, replace = F)
  keep_3 <- sample(1:niter, nsamp, replace = F)
  
  chain_1 <- x[[1]]
  sims_1 <- chain_1[keep_1, c(2, 3, 4, 1)]
  chain_2 <- x[[2]]
  sims_2 <- chain_2[keep_2, c(2, 3, 4, 1)]
  chain_3 <- x[[3]]
  sims_3 <- chain_3[keep_3, c(2, 3, 4, 1)]
  
  sims <- rbind(sims_1, sims_2, sims_3)
  
  steps <- numeric(length = nrow(sims))
  turns <- numeric(length = nrow(sims))
  
  for(q in 1:nrow(sims)){
    steps[q] <- rweibull(1, sims[q,3], (1/sims[q,4])^(1/sims[q,3]))
    turns[q] <- rwcauchy(1, sims[q,1], sims[q,2])
    #rwrappedcauchy(1, sims[i,1],  sims[i,2])
    }
  
  post_sims <- as.data.frame(cbind(steps, turns))
  post_sims$turns[post_sims$turns>pi]=post_sims$turns[post_sims$turns>pi]-2*pi

	mod <- post_sims
	mod$steps <- mod$steps*1000

	obs_mat <- matrix(nrow = n_obs, ncol = 2)
	obs_mat[1,1] <- init_x
	obs_mat[1,2] <- init_y
	obs_mat[2,1] <- (abs(rweibull(1, v, (1/lambda)^(1/v)))*-1) * 1000
	obs_mat[2,2] <- init_y
	
	for(m in 3:(n_obs)){
		
		obs_mat[m,1] <- obs_mat[m-1, 1] + (sin(mod$turns[m-2]) * mod$steps[m-2])
		obs_mat[m,2] <- obs_mat[m-1, 2] + (cos(mod$turns[m-2]) * mod$steps[m-2])
		}
	
	mov_df <- as.data.frame(cbind(obs_mat, c(seq(1:n_obs))))
	names(mov_df) <- c("x", "y", "obs")
	coordinates(mov_df) <- ~x+y
	mov_pts <- st_as_sf(mov_df) %>%
		st_set_crs("+proj=utm +zone=8 +datum=WGS84") 
	mov_ln <- mov_pts %>%
		dplyr::summarise(do_union = FALSE) %>%
		st_cast("LINESTRING") 
			
	whale_b <- st_buffer(mov_ln, 10)
	whale_end <- mov_pts %>%
    dplyr::filter(obs == 6) 
			
	whale_end_b <- st_buffer(whale_end, 10)
	
	strk_tmp <- st_intersects(whale_end_b, ship_strk_area_b, sparse = FALSE)
	strk <- ifelse(strk_tmp[1,1] == FALSE, 0, 1)
	
	mov_ln_num <- mov_pts[3:n_obs,] %>%
		dplyr::summarise(do_union = FALSE) %>%
		st_cast("LINESTRING") %>%
		mutate(path_num = i)

	path_sim[[i]] <- mov_ln_num
	strk_sim[i] <- strk
	}
	
	sum(strk_sim)
	
	
  
  paths <- do.call(rbind, path_sim)
		
	frst_step <- mov_pts[1:2,] %>%
		dplyr::summarise(do_union = FALSE) %>%
		st_cast("LINESTRING")
	
	frst_obs <- mov_pts[1,] 
	sec_obs <- mov_pts[2,]
	obs1_sp <- as(frst_obs, "Spatial")
	obs2_sp <- as(sec_obs, "Spatial")
	paths_sp <- as(paths, "Spatial")
	ship_init_sp <- as(ship_loc_sf[1,], "Spatial")
	

	
	area_rst <- raster(extent(paths_sp), crs = projection(paths_sp), ncols = 100, nrow = 100)
	#x <- rasterize(xy, r, fun='count')
  
  
	lengths <- sapply(1:ncell(area_rst), function(i){
	  tmp_rst <- area_rst
	  tmp_rst[i] <- 1
	  tmp_shp <- rasterToPolygons(tmp_rst)

	  if (gIntersects(paths_sp, tmp_shp)) {
		paths_sp_crp <- crop(paths_sp, tmp_shp)
		paths_sp_crp_length <- gLength(paths_sp_crp)
		return(paths_sp_crp_length)
	  } else {
		return(0)
	  }
	})
	
	
	area_rst[] <- lengths
	area_rst2 <- area_rst/(cellStats(area_rst, sum))
	ship_sp <- as(ship_ln, "Spatial")
	ship_strk_area_sp <- as(ship_strk_area_b, "Spatial")
	
	colfunc <- colorRampPalette(c("#CCCCCC", "black"))
	colfunc2 <- colorRampPalette(c("#F7EEEE", "darkred"))
	myTheme <- rasterTheme(region = c("white", colfunc2(15)))
	
	p <- levelplot(area_rst2, margin = FALSE, xlab = "Movement in X direction (m)",
            ylab = "Movement in Y direction (m)", par.settings = myTheme,
			xlim=c(init_x-2000, init_x+2000), ylim=c(init_y-2000, init_y+2000))
	p + layer(sp.lines(ship_strk_area_sp, lwd = 1, col='black')) + 
	layer(sp.points(ship_loc_sp, pch = 17, cex = 1, col = colfunc(7))) +
	layer(sp.points(obs1_sp, pch = 19, col='dodgerblue2')) +
	layer(sp.points(ship_init_sp, pch = 17, cex = 1, col='goldenrod2')) 	
	#layer(sp.lines(obs2_sp, col='gold'))
	
	
  
  
  
	r <- area_rst2
	r.min = cellStats(r, "min")
	r.max = cellStats(r, "max")
  
  rescale <- function(x, x.min = NULL, x.max = NULL, new.min = 0, new.max = 1) {
    if(is.null(x.min)) x.min = min(x)
    if(is.null(x.max)) x.max = max(x)
    new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
    }

	r.scale <- rescale(r, r.min, r.max, 0, 1)
	r.scale[is.na(r.scale)] <- 0

	myTheme <- rasterTheme(region = c("white", colfunc2(15)))
	
	p <- levelplot(r.scale, margin = FALSE, xlab = "Movement in X direction (m)",
            ylab = "Movement in Y direction (m)", par.settings = myTheme,
			xlim=c(-2000, 2000), ylim=c(-2000, 2000))
	p + layer(sp.lines(ship_strk_area_sp, lwd = 1, col='black')) + 
	layer(sp.points(ship_loc_sp, pch = 17, cex = 1, col = colfunc(7))) +
	layer(sp.points(obs1_sp, pch = 19, col='dodgerblue2')) +
	layer(sp.points(ship_init_sp, pch = 17, cex = 1, col='goldenrod2')) 	
	#layer(sp.lines(obs2_sp, col='gold'))
	
	
	
	
	
	
	
	
	
	ptsx <- c(-1000, 0, -1000)
	ptsy <- c(-354, 0, 354)
	df <- as.data.frame(cbind(ptsx, ptsy))

	pts <- df
	coordinates(pts) <- ~ptsx+ptsy

	pts_sf <- st_as_sf(mov_df) %>%
		st_set_crs("+proj=utm +zone=8 +datum=WGS84") 
		pts_sf <- st_as_sf(pts) %>%
		st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
		dplyr::summarise(do_union = FALSE) %>%
		st_cast("LINESTRING")

	ln1 <- as(pts_sf, 'Spatial')
	
		
	p <- levelplot(area_rst2, margin = FALSE, xlab = "Movement in X direction (m)",
            ylab = "Movement in Y direction (m)", par.settings = myTheme,
			xlim=c(-2000, 2000), ylim=c(-2000, 2000))
	p + layer(sp.lines(ship_strk_area_sp, lwd = 1, col='black')) + 
	layer(sp.points(ship_loc_sp, pch = 17, cex = 1, col = colfunc(7))) +
	layer(sp.points(obs1_sp, pch = 19, col='dodgerblue2')) +
	#layer(sp.lines(ln1, col='dodgerblue2')) 	
	#layer(sp.lines(obs2_sp, col='gold'))
	
  
  
  
  
  
  #  Background set up
  #  13 knots = 6.68778 m/s
  #  19 knots = 9.77444 m/s
  ship_spd <- 9.77444
  whale_spd <-  2.84939 #6.68778
  btw_surf <- 30
  mov_time <- 180
  ship_vert_dist <- mov_time * ship_spd
  
  perp_dist_ship_whale <- 230
  vert_dist_ship_whale <- ship_vert_dist
  
  # Whale starting location
  init_x <- 0 #mu_x
  init_y <- 0 #mu_y
  
  # Number of whale surfacings
  n_obs <- 7 #roundup(mov_time/btw_surf)
  
  # Ship starting location
  init_ship_x <- init_x - perp_dist_ship_whale
  init_ship_y <- init_y - vert_dist_ship_whale
  
  # Ship movement
  num_ship_pts <- n_obs
  
  ship_init_loc_df <- as.data.frame(cbind(init_ship_x, init_ship_y))
  colnames(ship_init_loc_df) <- c("x", "y")
  
    
    
  ship_turn_seq <- seq(0, -0.35, -0.035)
  ship_turn_strk_sim <- numeric(length(ship_turn_seq))
  all_paths_ls <- list()
  
  for(t in 1:10){
    
    ship_turn <- ship_turn_seq[t]
    ship_dat <- as.data.frame(cbind(c(0, rep(ship_turn, 6)),
      rep((ship_spd * 30), 7)))
    colnames(ship_dat) <- c("turns", "steps")
    
    ship_mat <- matrix(nrow = 7, ncol = 2)
    ship_mat[1,1] <- init_ship_x
    ship_mat[1,2] <- init_ship_y
  
    for(m in 2:7){
      
      ship_mat[m,1] <- ship_mat[m-1, 1] + (sin(ship_dat$turns[m]) * ship_dat$steps[m])
      ship_mat[m,2] <- ship_mat[m-1, 2] + (cos(ship_dat$turns[m]) * ship_dat$steps[m])
      }
      
    ship_loc_df <- as.data.frame(ship_mat)
    colnames(ship_loc_df) <- c("x", "y")
    coordinates(ship_loc_df) <- ~x+y
    ship_loc_sf <- st_as_sf(ship_loc_df)  %>%
      st_set_crs("+proj=utm +zone=8 +datum=WGS84") %>%
      mutate(obs_num = 1:n())
    ship_loc_sp <- as(ship_loc_sf, "Spatial")
    
    ship_ln <- st_as_sf(ship_loc_df) %>%
      dplyr::summarise(do_union = FALSE) %>%
      st_cast("LINESTRING") %>%
      st_set_crs("+proj=utm +zone=8 +datum=WGS84")
    
    ship_strk_area <- st_as_sf(ship_loc_sp[7,])
      #st_as_sf(ship_loc_sp[6:7,]) %>%
      #dplyr::summarise(do_union = FALSE) %>%
      #st_cast("LINESTRING") 
      
    ship_strk_area_b <- st_buffer(ship_strk_area, 50)
    ship_ln_b <- st_buffer(ship_ln, 10)
    
      
    n_strk_sim <- 10000
    strk_sim <- as.vector(rep(NA, length = n_strk_sim))
    path_sim <- list()
    
    for(i in 1:n_strk_sim){
      
        x <- single_fit_surf
        keep_1 <- sample(2500:7500, 1000, replace = F)
        keep_2 <- sample(2500:7500, 1000, replace = F)
        keep_3 <- sample(2500:7500, 1000, replace = F)
        
        chain_1 <- x[[1]]
        sims_1 <- chain_1[keep_1, c(2, 3, 4, 1)]
        chain_2 <- x[[2]]
        sims_2 <- chain_2[keep_2, c(2, 3, 4, 1)]
        chain_3 <- x[[3]]
        sims_3 <- chain_3[keep_3, c(2, 3, 4, 1)]
        
        sims_tmp <- rbind(sims_1, sims_2, sims_3)
        keep_tmp <- sample(1:nrow(sims_tmp), nobs-2, replace = FALSE)
        sims <- sims_tmp[keep_tmp, ]
        
        steps <- numeric(length = nrow(sims))
        turns <- numeric(length = nrow(sims))
        
        for(s in 1:nrow(sims)){
          steps[s] <- rweibull(1, sims[s,3], (1/sims[s,4])^(1/sims[s,3]))
          turns[s] <- rwcauchy(1, sims[s,1], sims[s,2])
          #rwrappedcauchy(1, sims[i,1],  sims[i,2])
          }
    
      mod <- as.data.frame(cbind(steps, turns))
      mod$turns[mod$turns>pi]=mod$turns[mod$turns>pi]-2*pi
      mod$steps <- mod$steps*1000
    
      obs_mat <- matrix(nrow = nobs, ncol = 2)
      obs_mat[1,1] <- init_x
      obs_mat[1,2] <- init_y
      obs_mat[2,1] <- (abs(rweibull(1, v, (1/lambda)^(1/v)))*-1) * 1000
      obs_mat[2,2] <- init_y
      
      for(m in 3:(nobs)){
        
        obs_mat[m,1] <- obs_mat[m-1, 1] + (sin(mod$turns[m-2]) * mod$steps[m-2])
        obs_mat[m,2] <- obs_mat[m-1, 2] + (cos(mod$turns[m-2]) * mod$steps[m-2])
        }
      
      mov_df <- as.data.frame(cbind(obs_mat, c(seq(1:nobs))))
      names(mov_df) <- c("x", "y", "obs")
      coordinates(mov_df) <- ~x+y
      mov_pts <- st_as_sf(mov_df) %>%
        st_set_crs("+proj=utm +zone=8 +datum=WGS84") 
      mov_ln <- mov_pts %>%
        dplyr::summarise(do_union = FALSE) %>%
        st_cast("LINESTRING") 
          
    
      whale_b <- st_buffer(mov_ln, 10)
      whale_end <- mov_pts %>%
        dplyr::filter(obs > 5) 
        #dplyr::filter(obs == 7) #%>%
        #dplyr::summarise(do_union = FALSE) %>%
        #st_cast("LINESTRING") 
          
      whale_end_b <- st_buffer(whale_end, 10)
      
      strk_tmp <- st_intersects(whale_end_b, ship_strk_area_b, sparse = FALSE)
      #strk <- ifelse(strk_tmp[1,1] == FALSE, 0, 1)
      strk <- ifelse(strk_tmp[1,1] == FALSE & strk_tmp[2,1] == FALSE, 0,
        ifelse(strk_tmp[1,1] == FALSE & strk_tmp[2,1] == TRUE, 1,
        ifelse(strk_tmp[1,1] == TRUE & strk_tmp[2,1] == FALSE, 1, 1)))
      
      mov_ln_num <- mov_pts[3:nobs,] %>%
        dplyr::summarise(do_union = FALSE) %>%
        st_cast("LINESTRING") %>%
        mutate(path_num = i)
    
      path_sim[[i]] <- mov_ln_num
      strk_sim[i] <- strk
      }
      
      tot_strk <- sum(strk_sim)
      
      ship_turn_strk_sim[t] <- tot_strk
      paths <- do.call(rbind, path_sim)
      all_paths_ls[[t]] <- paths
      }
    
  
 
 
 
 
 
		
	frst_step <- mov_pts[1:2,] %>%
		dplyr::summarise(do_union = FALSE) %>%
		st_cast("LINESTRING")
	
	frst_obs <- mov_pts[1,] 
	sec_obs <- mov_pts[2,]
	obs1_sp <- as(frst_obs, "Spatial")
	obs2_sp <- as(sec_obs, "Spatial")
	paths_sp <- as(paths, "Spatial")
	ship_init_sp <- as(ship_loc_sf[1,], "Spatial")
	

	
	area_rst <- raster(extent(paths_sp), crs = projection(paths_sp), ncols = 100, nrow = 100)
	
	lengths <- sapply(1:ncell(area_rst), function(i){
	  tmp_rst <- area_rst
	  tmp_rst[i] <- 1
	  tmp_shp <- rasterToPolygons(tmp_rst)

	  if (gIntersects(paths_sp, tmp_shp)) {
		paths_sp_crp <- crop(paths_sp, tmp_shp)
		paths_sp_crp_length <- gLength(paths_sp_crp)
		return(paths_sp_crp_length)
	  } else {
		return(0)
	  }
	})
	
	
	area_rst[] <- lengths
	area_rst2 <- area_rst/(cellStats(area_rst, sum))
	ship_sp <- as(ship_ln, "Spatial")
	ship_strk_area_sp <- as(ship_strk_area_b, "Spatial")
	
	colfunc <- colorRampPalette(c("#CCCCCC", "black"))
	colfunc2 <- colorRampPalette(c("#F7EEEE", "darkred"))
	myTheme <- rasterTheme(region = c("white", colfunc2(15)))
	
	p <- levelplot(area_rst2, margin = FALSE, xlab = "Movement in X direction (m)",
            ylab = "Movement in Y direction (m)", par.settings = myTheme,
			xlim=c(-2000, 2000), ylim=c(-2000, 2000))
	p + layer(sp.lines(ship_strk_area_sp, lwd = 1, col='black')) + 
	layer(sp.points(ship_loc_sp, pch = 17, cex = 1, col = colfunc(7))) +
	layer(sp.points(obs1_sp, pch = 19, col='dodgerblue2')) +
	layer(sp.points(ship_init_sp, pch = 17, cex = 1, col='goldenrod2')) 	
	#layer(sp.lines(obs2_sp, col='gold'))
	
	
	r <- area_rst2
	r.min = cellStats(r, "min")
	r.max = cellStats(r, "max")
  
  rescale <- function(x, x.min = NULL, x.max = NULL, new.min = 0, new.max = 1) {
    if(is.null(x.min)) x.min = min(x)
    if(is.null(x.max)) x.max = max(x)
    new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
    }

	r.scale <- rescale(r, r.min, r.max, 0, 1)
	r.scale[is.na(r.scale)] <- 0

	myTheme <- rasterTheme(region = c("white", colfunc2(15)))
	
	p <- levelplot(r.scale, margin = FALSE, xlab = "Movement in X direction (m)",
            ylab = "Movement in Y direction (m)", par.settings = myTheme,
			xlim=c(-2000, 2000), ylim=c(-2000, 2000))
	p + layer(sp.lines(ship_strk_area_sp, lwd = 1, col='black')) + 
	layer(sp.points(ship_loc_sp, pch = 17, cex = 1, col = colfunc(7))) +
	layer(sp.points(obs1_sp, pch = 19, col='dodgerblue2')) +
	layer(sp.points(ship_init_sp, pch = 17, cex = 1, col='goldenrod2')) 	
	#layer(sp.lines(obs2_sp, col='gold'))
	
	
