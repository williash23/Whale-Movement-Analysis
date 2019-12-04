mround <- function(x,base){ 
        base*round(x/base) 
} 

acumulated.distrib= function(sample,x){
    minors= 0
    for(n in sample){
        if(n<=x){
            minors= minors+1
        }
    }
    return (minors/length(sample))
}



ecdf_plot_all <- post_sims_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) 
ecdf_plot_su <- post_sims_su_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) 
ecdf_plot_d <- post_sims_d_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) 
ecdf_plot_n <- post_sims_n_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) 
ecdf_plot_f <- post_sims_f_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) 
ecdf_plot_fr <- post_sims_fr_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) 
ecdf_plot_si <- post_sims_si_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) 


ecdf_turns_75_su <- as.numeric(quantile(ecdf(ecdf_plot_su$abs_turns_deg), probs = 0.75)) 
ecdf_turns_75_d <- as.numeric(quantile(ecdf(ecdf_plot_d$abs_turns_deg), probs = 0.75)) 


ecdf_steps_75_n <- as.numeric(quantile(ecdf(ecdf_plot_n$steps), probs = 0.75)) 
ecdf_steps_75_f <- as.numeric(quantile(ecdf(ecdf_plot_f$steps), probs = 0.75)) 


acumulated.distrib(ecdf_plot_all$abs_turns_deg, 45)

par(mfrow = c(2,2))
hist(ecdf_plot_n$abs_turns_deg, breaks = 10)
hist(ecdf_plot_n$turns_deg, breaks = 10)
plot(ecdf(ecdf_plot_n$abs_turns_deg))
plot(ecdf(ecdf_plot_n$turns_deg))

ecdf_step_50 <- as.numeric(quantile(ecdf(ecdf_plot_all$steps), probs = 0.5)) * 1000
ecdf_step_75 <- as.numeric(quantile(ecdf(ecdf_plot_all$steps), probs = 0.75)) * 1000
ecdf_step_95 <- as.numeric(quantile(ecdf(ecdf_plot_all$steps), probs = 0.95)) *1000

ecdf_turns_50 <- as.numeric(quantile(ecdf(ecdf_plot_n$turns_deg), probs = 0.5)) 
ecdf_turns_75 <- as.numeric(quantile(ecdf(ecdf_plot_all$abs_turns_deg), probs = 0.75)) 


ship_dist <- 2500

avoid_turn_95 <- NISTradianTOdeg(atan(ecdf_step_95/ship_dist))
avoid_turn_75 <- NISTradianTOdeg(atan(ecdf_step_75/ship_dist))
avoid_turn_50 <- NISTradianTOdeg(atan(ecdf_step_50/ship_dist))

avoid_dist_95 <- sqrt(ecdf_step_95^2 + ship_dist^2)

start_x <- 433525.87
start_y <- 6507454.23
start_x <- start_x + 200
start_y <- start_y + 1250

#par(mfrow = c(1,2))

plot(start_x, start_y, type="n",
	ylim = c(6504000,6514000),
	xlim = c(431000, 435000),
	xlab = "Longitude (UTM)", ylab = "Latitude (UTM)")
	
points(start_x, start_y, pch = 19, col = "black", cex = 1, lwd = 2)

draw.arc(3.66955475, 2.48709418, x = start_x,  y =start_y, radius = ecdf_step_50)	
draw.arc(3.66955475, 2.48709418, x = start_x,  y =start_y, radius = ecdf_step_75)
draw.arc(3.66955475, 2.48709418, x = start_x,  y =start_y, radius = ecdf_step_95)
segments(start_x, start_y, start_x - 2000 , start_y, lty = 3, lwd = 2)

segments(start_x, start_y, 432650, 6511100, lty = 2, lwd = 2)
segments(start_x, start_y, 432620, 6506700, lty = 2, lwd = 2)

#draw.sector(315, 225, x = start_x, y = start_y), rou1 = ecdf_step_75)

circ_dat <- draw.circle(start_x, start_y, radius = c(ecdf_step_95, ecdf_step_75, ecdf_step_50), 
	col = c("grey90", "grey70", "grey50"), border = "black", lty = 1, lwd = 1)
points(start_x, start_y, pch = 19, col = "black", cex = 1.5, lwd = 2)

segments(start_x, start_y, 433000, 6511600, lty = 3, lwd = 2)
segments(start_x, start_y - ship_dist, start_x - ecdf_step_75, start_y, lty = 3, lwd = 2)
segments(start_x, start_y - ship_dist, start_x - ecdf_step_50, start_y, lty = 3, lwd = 2)
#segments(start_x, start_y + ship_dist, start_x + ecdf_step_95, start_y, lty = 3, lwd = 2)
#segments(start_x, start_y + ship_dist, start_x + ecdf_step_75, start_y, lty = 3, lwd = 2)
#segments(start_x, start_y + ship_dist, start_x + ecdf_step_50, start_y, lty = 3, lwd = 2)
segments(start_x, start_y - ship_dist, start_x, start_y, lwd = 2)
segments(start_x - ecdf_step_95, start_y, start_x + ecdf_step_95, start_y, lwd = 2)

points(start_x, start_y - ship_dist, pch = 19, col = "red", cex = 1.5, lwd = 2)

text(start_x, start_y + ecdf_step_95 + 460, "95%", cex = 0.6)
text(start_x, start_y + ecdf_step_75 + 180, "75%", cex = 0.6)
text(start_x, start_y + ecdf_step_50 + 65, "50%", cex = 0.6)
#text(start_x - 50, 6508550, "2500 m", cex = 0.6, srt = 90)
text(start_x - 150, start_y + 75, "297", cex = 0.6)
text(start_x - 415, start_y + 75, "564", cex = 0.6)
text(start_x - 850, start_y + 75, "1163 m", cex = 0.6)
text(start_x - 50, start_y  - 1200, "2500 m", cex = 0.6, srt = 90)



plot(loc_df2$X, loc_df2$Y, type="n",
	#ylim = c(6503000,6514000),
	#xlim = c(431000, 435000),
	xlab = "Longtidue (UTM)", ylab = "Latitude (UTM)")
circ_dat <- draw.circle(start_x, start_y, radius = c(ecdf_step_95, ecdf_step_75, ecdf_step_50), 
	col = c("grey90", "grey70", "grey50"), border = "black", lty = 1, lwd = 1)
points(loc_df2$X, loc_df2$Y, pch = 19, col = rgb(0, 0, 0, 0.2), cex = 0.5)
points(start_x, start_y, pch = 19, col = "black", cex = 1.5, lwd = 2)
points(start_x, start_y - ship_dist, pch = 19, col = "red", cex = 1.5, lwd = 2)


turns.ordered <-  sort(ecdf_plot_su $abs_turns_deg)

plot(turns.ordered, (1:n)/n, type = 's', ylim = c(0, 1)) #xlab = 'Sample Quantiles of Ozone', ylab = '', main = 'Empirical Cumluative Distribution\nOzone Pollution in New York')


ecdf_plot_all <- post_sims_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) 
ecdf_plot_su <- post_sims_su_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) 
ecdf_plot_d <- post_sims_d_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) 
ecdf_plot_n <- post_sims_n_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) 
ecdf_plot_f <- post_sims_f_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) 
ecdf_plot_fr <- post_sims_fr_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) 
ecdf_plot_si <- post_sims_si_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) 


	
par(mfrow = c(1,2))
#box(bty="l")
par(yaxs="i",las=1)
plot(ecdf(ecdf_plot_all$steps),col="black",lwd=3, ylab = "Empirical cumulative probability", 
	xlab="Step length (km)", xlim = c(0, 2.5), ylim = c(0, 1), main=NULL, xaxs = "i")
segments(quantile(ecdf_plot_all$steps, probs = 0.5), 0.5, quantile(ecdf_plot_all$steps, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_all$steps, probs = 0.5), 0.5, 0, 0.5, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_all$steps, probs = 0.75), 0.75, quantile(ecdf_plot_all$steps, probs = 0.75), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_all$steps, probs = 0.75), 0.75, 0, 0.75, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_all$steps, probs = 0.95), 0.95, quantile(ecdf_plot_all$steps, probs = 0.95), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_all$steps, probs = 0.95), 0.95, 0, 0.95, 
	lty = 2, lwd = 1.5)
axis(side = c(1,2), lty = "solid")
abline(h = 1, lty = "solid")

#box(bty="l")
par(yaxs="i",las=1)
plot(ecdf(ecdf_plot_all$abs_turns_deg),col="black",lwd=3, ylab = "", 
	xlab="Turn angle (deg)", xlim = c(0, 180), ylim = c(0, 1), main=NULL, xaxs = "i")
segments(quantile(ecdf_plot_all$abs_turns_deg, probs = 0.5), 0.5, quantile(ecdf_plot_all$abs_turns_deg, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_all$abs_turns_deg, probs = 0.5), 0.5, 0, 0.5, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_all$abs_turns_deg, probs = 0.75), 0.75, quantile(ecdf_plot_all$abs_turns_deg, probs = 0.75), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_all$abs_turns_deg, probs = 0.75), 0.75, 0, 0.75, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_all$abs_turns_deg, probs = 0.95), 0.95, quantile(ecdf_plot_all$abs_turns_deg, probs = 0.95), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_all$abs_turns_deg, probs = 0.95), 0.95, 0, 0.95, 
	lty = 2, lwd = 1.5)
axis(side = c(1,2), lty = "solid")
abline(h = 1, lty = "solid")






par(mfrow = c(1,2))
#box(bty="l")
par(yaxs="i",las=1)
plot(ecdf(ecdf_plot_d$steps),col="black",lwd=3, 	ylab = "Empirical cumulative probability", 
	xlab="Step length (km)", xlim = c(0, 2.5), ylim = c(0, 1), main=NULL, xaxs = "i")
segments(quantile(ecdf_plot_d$steps, probs = 0.5), 0.5, quantile(ecdf_plot_d$steps, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_d$steps, probs = 0.5), 0.5, 0, 0.5, 
	lty = 2, lwd = 1.5)
lines(ecdf(ecdf_plot_su$steps),col="black",lwd=1)
segments(quantile(ecdf_plot_su$steps, probs = 0.5), 0.5, quantile(ecdf_plot_su$steps, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_su$steps, probs = 0.5), 0.5, 0, 0.5, 
	lty = 2, lwd = 1.5)
axis(side = c(1,2), lty = "solid")
abline(h = 1, lty = "solid")


#box(bty="l")
par(yaxs="i",las=1)
plot(ecdf(ecdf_plot_d$abs_turns_deg),col="black",lwd=3, 	ylab = "", 
	xlab="Turn angle (deg)", xlim = c(0, 180), ylim = c(0, 1), main=NULL, xaxs = "i")
segments(quantile(ecdf_plot_d$abs_turns_deg, probs = 0.5), 0.5,
    quantile(ecdf_plot_d$abs_turns_deg, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_d$abs_turns_deg, probs = 0.5), 0.5, 0, 0.5, 
	lty = 2, lwd = 1.5)
lines(ecdf(ecdf_plot_su$abs_turns_deg),col="black",lwd=1)
segments(quantile(ecdf_plot_su$abs_turns_deg, probs = 0.5), 0.5, quantile(ecdf_plot_su$abs_turns_deg, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_su$abs_turns_deg, probs = 0.5), 0.5, 0, 0.5, 
	lty = 2, lwd = 1.5)
axis(side = c(1,2), lty = "solid")
abline(h = 1, lty = "solid")







par(mfrow = c(1,2))

#box(bty="l")
par(yaxs="i",las=1)
plot(ecdf(ecdf_plot_n$steps),col="black",lwd=3, 	ylab = "Empirical cumulative probability", 
	xlab="Step length (km)", xlim = c(0, 2.5), ylim = c(0, 1), main=NULL, xaxs = "i")
segments(quantile(ecdf_plot_n$steps, probs = 0.5), 0.5, quantile(ecdf_plot_n$steps, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_n$steps, probs = 0.5), 0.5, 0, 0.5, 
	lty = 2, lwd = 1.5)
lines(ecdf(ecdf_plot_f$steps),col="black",lwd=1)
segments(quantile(ecdf_plot_f$steps, probs = 0.5), 0.5, quantile(ecdf_plot_f$steps, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_f$steps, probs = 0.5), 0.5, 0, 0.5, 
	lty = 2, lwd = 1.5)
axis(side = c(1,2), lty = "solid")
abline(h = 1, lty = "solid")


#box(bty="l")
par(yaxs="i",las=1)
plot(ecdf(ecdf_plot_n$abs_turns_deg),col="black",lwd=3, 	ylab = "", 
	xlab="Turn angle (deg)", xlim = c(0, 180), ylim = c(0, 1), main=NULL, xaxs = "i")
segments(quantile(ecdf_plot_n$abs_turns_deg, probs = 0.5), 0.5, quantile(ecdf_plot_n$abs_turns_deg, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
#segments(quantile(ecdf_plot_n$turns_deg, probs = 0.5), 0.5, 0, 0.5, 
	#lty = 2, lwd = 1.5)
lines(ecdf(ecdf_plot_f$abs_turns_deg),col="black",lwd=1)
segments(quantile(ecdf_plot_f$abs_turns_deg, probs = 0.5), 0.5, quantile(ecdf_plot_f$abs_turns_deg, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_f$abs_turns_deg, probs = 0.5), 0.5, 0, 0.5, 
	lty = 2, lwd = 1.5)
axis(side = c(1,2), lty = "solid")
abline(h = 1, lty = "solid")







par(mfrow = c(1,2))

par(yaxs="i",las=1)
plot(ecdf(ecdf_plot_fr$steps),col="black",lwd=3, 	ylab = "Empirical cumulative probability", 
	xlab="Step length (km)", xlim = c(0, 2.5), ylim = c(0, 1), main=NULL, xaxs = "i")
segments(quantile(ecdf_plot_fr$steps, probs = 0.5), 0.5, quantile(ecdf_plot_fr$steps, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_fr$steps, probs = 0.5), 0.5, 0, 0.5, 
	lty = 2, lwd = 1.5)
lines(ecdf(ecdf_plot_si$steps),col="black",lwd=1)
axis(side = c(1,2), lty = "solid")
abline(h = 1, lty = "solid")


par(yaxs="i",las=1)
plot(ecdf(ecdf_plot_fr$abs_turns_deg),col="black",lwd=3, 	ylab = "", 
	xlab="Turn angle (deg)", xlim = c(0, 180), ylim = c(0, 1), main=NULL, xaxs = "i")
lines(ecdf(ecdf_plot_si$abs_turns_deg),col="black",lwd=1)
segments(quantile(ecdf_plot_si$abs_turns_deg, probs = 0.5), 0.5, quantile(ecdf_plot_si$abs_turns_deg, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_si$abs_turns_deg, probs = 0.5), 0.5, 0, 0.5, 
	lty = 2, lwd = 1.5)
axis(side = c(1,2), lty = "solid")
abline(h = 1, lty = "solid")













plot(ecdf(ecdf_plot_fr$turns_deg),col="black",lwd=3, 	ylab = "Empirical cumulative probability", 
	xlab="Turn angle (deg)", xlim = c(0, 360), ylim = c(0, 1), main=NULL, xaxs = "i")
segments(quantile(ecdf_plot_fr$turns_deg, probs = 0.5), 0.5, quantile(ecdf_plot_fr$turns_deg, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
#segments(quantile(ecdf_plot_fr$turns_deg, probs = 0.5), 0.5, 0, 0.5, 
	#lty = 2, lwd = 1.5)
plot(ecdf(ecdf_plot_si$abs_turns,col="black",lwd=1)
segments(quantile(ecdf_plot_si$abs_turns, probs = 0.5), 0.5, quantile(ecdf_plot_si$turns_deg, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(ecdf_plot_si$abs_turns, probs = 0.5), 0.5, 0, 0.5, 
	lty = 2, lwd = 1.5)
axis(side = c(1,2), lty = "solid")
abline(h = 1, lty = "solid")

















par(mfrow = c(1,2), yaxs="i",las=1, omi = c(0.1, 0.1, 0.1, 0.1) )
hist(post_sims_plot$steps, prob=TRUE,col="grey80",border="white", xlab="Step length (km)",
	ylab = "Density", xlim = c(0, 2.5), ylim = c(0, 2.5), breaks =  seq(0, 4.5, .125) , main=NULL)
box(bty="l")
lines(density(post_sims_plot$steps),col="black",lwd=2)
abline(v = mean(post_sims_plot$steps), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")	
	
	
#par(yaxs="i",las=1)
hist(post_sims_plot$turns_deg, prob=TRUE,col="grey80",border="white", xlab="Turn angle (deg)",
	ylab = NULL, xlim = c(-180, 180), ylim = c(0, 0.005), breaks = 25, main=NULL)
#main="Distribution of Respirable Particle Concentrations")
box(bty="l")
lines(density(post_sims_plot$turns_deg),col="black",lwd=2)
abline(v = mean(post_sims_plot$turns_deg), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")


par(mfrow = c(2,2), yaxs="i",las=1, omi = c(0.1, 0.1, 0.1, 0.1) )
box(bty="l")
#par(yaxs="i",las=1)
plot(ecdf(post_sims_plot$steps),col="black",lwd=2, 	ylab = "Empirical cumulative probability", 
	xlab="Step length (km)", xlim = c(0, 2.5), ylim = c(0, 1), main=NULL, xaxs = "i")
abline(v = quantile(post_sims_plot$steps, probs = c(seq(0, 1, 0.25), 0.95)), lty = 2)
#lines(ecdf(post_sims_plot$steps),col="black",lwd=2, lty = 2)

#box(bty="l")
#par(yaxs="i",las=1)
plot(ecdf(post_sims_plot$turns_deg),col="black",lwd=2, 	ylab = "Empirical cumulative probability", 
	xlab="Turn angle (deg)", ylab = NULLxlim = c(-180, 180), ylim = c(0, 1), main=NULL, xaxs = "i")
abline(v = quantile(post_sims_plot$turns_deg, probs = c(seq(0, 1, 0.25), 0.95)), lty = 2)
#lines(ecdf(post_sims_plot$turns_deg),col="black",lwd=2, lty = 2)


box(bty="l")
par(mfrow = c(2,2), yaxs="i",las=1, omi = c(0.1, 0.1, 0.1, 0.1) )
#par(yaxs="i",las=1)
hist(post_sims_d_plot$steps, prob=TRUE,col="grey80",border="white", xlab="Step length (km)",
	ylab = "Density", xlim = c(0, 2.5), ylim = c(0, 1), breaks =  seq(0, 4.5, .125) , main=NULL)
lines(density(post_sims_d_plot$steps),col="black",lwd=2)
abline(v = mean(post_sims_d_plot$steps), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")	
	
#par(yaxs="i",las=1)
hist(post_sims_d_plot$turns_deg, prob=TRUE,col="grey80",border="white", xlab="Turn angle (deg)",
	ylab =NULL, xlim = c(-180, 180), ylim = c(0, 0.005), breaks = 25, main=NULL)
#box(bty="l")
lines(density(post_sims_d_plot$turns_deg),col="black",lwd=2)
abline(v = mean(post_sims_d_plot$turns_deg), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")
	
#par(yaxs="i",las=1)
hist(post_sims_su_plot$steps, prob=TRUE,col="grey80",border="white", xlab="Step length (km)",
	ylab = "Density", xlim = c(0, 2.5), ylim = c(0, 10), breaks =  seq(0, 4.5, .125), main=NULL)
#box(bty="l")
lines(density(post_sims_su_plot$steps),col="black",lwd=2)
abline(v = mean(post_sims_su_plot$steps), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")	

#par(yaxs="i",las=1)
hist(post_sims_su_plot$turns_deg, prob=TRUE,col="grey80",border="white", xlab="Turn angle (deg)",
	ylab = NULL, xlim = c(-180, 180), ylim = c(0, 0.005), breaks = 25, main=NULL)
#box(bty="l")
lines(density(post_sims_su_plot$turns_deg),col="black",lwd=2)
abline(v = mean(post_sims_su_plot$turns_deg), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")



box(bty="l")
par(yaxs="i",las=1)
plot(ecdf(post_sims_d_plot$steps),col="black",lwd=3, 	ylab = "Empirical cumulative probability", 
	xlab="Step length (km)", xlim = c(0, 2.5), ylim = c(0, 1), main=NULL, xaxs = "i")
segments(quantile(post_sims_d_plot$steps, probs = 0.5), 0.5, quantile(post_sims_d_plot$steps, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(post_sims_d_plot$steps, probs = 0.5), 0.5, 0, 0.5, 
	lty = 2, lwd = 1.5)
lines(ecdf(post_sims_su_plot$steps),col="black",lwd=3, lty = 5)
segments(quantile(post_sims_su_plot$steps, probs = 0.5), 0.5, quantile(post_sims_su_plot$steps, probs = 0.5), 0, 
	lty = 2, lwd = 1.5)
segments(quantile(post_sims_su_plot$steps, probs = 0.5), 0.5, 0, 0.5, 
	lty = 2, lwd = 1.5)


box(bty="l")
par(yaxs="i",las=1)
plot(ecdf(post_sims_plot$turns_deg),col="black",lwd=2, 	ylab = "Empirical cumulative probability", 
	xlab="Turn angle (deg)", xlim = c(-180, 180), ylim = c(0, 1), main=NULL, xaxs = "i")
abline(v = quantile(post_sims_plot$turns_deg, probs = c(seq(0, 1, 0.25), 0.95)), lty = 2)
#lines(ecdf(post_sims_plot$turns_deg),col="black",lwd=2, lty = 2)








box(bty="l")
par(mfrow = c(2,2), yaxs="i",las=1, omi = c(0.1, 0.1, 0.1, 0.1) )
par(yaxs="i",las=1)
hist(post_sims_n_plot$steps, prob=TRUE,col="grey80",border="white", xlab="Step length (km)",
	ylab = "Density", xlim = c(0, 2.5), ylim = c(0, 6), breaks = seq(0, 4.5, .125), main=NULL)
box(bty="l")
lines(density(post_sims_n_plot$steps),col="black",lwd=2)
abline(v = mean(post_sims_n_plot$steps), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")	
	
par(yaxs="i",las=1)
hist(post_sims_n_plot$turns_deg, prob=TRUE,col="grey80",border="white", xlab="Turn angle (deg)",
	ylab = NULL, xlim = c(-180, 180), ylim = c(0, 0.005), breaks = 25, main=NULL)
#main="Distribution of Respirable Particle Concentrations")
box(bty="l")
lines(density(post_sims_n_plot$turns_deg),col="black",lwd=2)
abline(v = mean(post_sims_n_plot$turns_deg), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")

#box(bty="l")
#par(mfrow = c(2,2), yaxs="i",las=1, omi = c(0.1, 0.1, 0.1, 0.1) )
hist(post_sims_f_plot$steps, prob=TRUE,col="grey80",border="white", xlab="Step length (km)",
	ylab = "Density", xlim = c(0, 2.5), ylim = c(0, 1.2), breaks =  seq(0, 4.5, .125), main=NULL)
#box(bty="l")
lines(density(post_sims_f_plot$steps),col="black",lwd=2)
abline(v = mean(post_sims_f_plot$steps), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")	
	
#par(yaxs="i",las=1)
hist(post_sims_f_plot$turns_deg, prob=TRUE,col="grey80",border="white", xlab="Turn angle (deg)",
	ylab = NULL, xlim = c(-180, 180), ylim = c(0, 0.005), breaks = 25, main=NULL)
#box(bty="l")
lines(density(post_sims_f_plot$turns_deg),col="black",lwd=2)
abline(v = mean(post_sims_f_plot$turns_deg), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")





#turns_side_1k <- cbind(NISTradianTOdeg(abs(theta_side_1k)), rep(0.00008, length(theta_side_1k)))
#turns_front_1k <- cbind(NISTradianTOdeg(abs(theta_front_1k)), rep(0.0001, length(theta_front_1k)))
#steps_side_1k <- cbind(l_side_1k, rep(0.01, length(l_side_1k)))
#steps_front_1k <- cbind(l_front_1k, rep(0.01, length(l_front_1k)))


    
box(bty="l")
par(mfrow = c(2,2), yaxs="i",las=1, omi = c(0.1, 0.1, 0.1, 0.1) )
par(yaxs="i",las=1)
hist(post_sims_si_plot$steps[post_sims_si_plot$steps < 2.5], prob=TRUE,col="grey80",border="white", xlab="Step length (km)",
	ylab = "Density", xlim = c(0, 2.5), ylim = c(0, 2.5), breaks =  seq(0, 4.5, .125), main=NULL)
box(bty="l")
lines(density(post_sims_si_plot$steps),col="black",lwd=2)
#points(steps_side_1k, col = "red", pch="|")
abline(v = mean(post_sims_si_plot$steps), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")	
	
#par(yaxs="i",las=1)
hist(post_sims_si_plot$turns_deg, prob=TRUE,col="grey80",border="white", xlab="Turn angle (deg)",
	ylab = NULL, xlim = c(-180, 180), ylim = c(0, 0.005), breaks = 25, main=NULL)
#box(bty="l")
lines(density(post_sims_si_plot$turns_deg),col="black",lwd=2)
#points(turns_side_1k, col = "red", pch="|")
abline(v = mean(post_sims_si_plot$turns_deg), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")


#par(yaxs="i",las=1)
hist(post_sims_fr_plot$steps, prob=TRUE,col="grey80",border="white", xlab="Step length (km)",
	ylab = "Density", xlim = c(0, 2.5),ylim = c(0, 2.5), breaks =  seq(0, 4.5, .125), main=NULL)
#box(bty="l")
lines(density(post_sims_fr_plot$steps),col="black",lwd=2)
#points(steps_front_1k, col = "red", pch="|")
abline(v = mean(post_sims_fr_plot$steps), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")	
	
#par(yaxs="i",las=1)
hist(post_sims_fr_plot$turns_deg, prob=TRUE,col="grey80",border="white", xlab="Turn angle (deg)",
	ylab =NULL, xlim = c(-180, 180),  ylim = c(0, 0.005), breaks = 25, main=NULL)
#box(bty="l")
lines(density(post_sims_fr_plot$turns_deg),col="black",lwd=2)
#points(turns_front_1k, col = "red", pch="|")
abline(v = mean(post_sims_fr_plot$turns_deg), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")









#######
par(yaxs="i",las=1)
hist(post_sims_d_plot$steps, prob=TRUE,col="grey60",border="white", xlab="Step length (km)",
	ylab = "Density", xlim = c(0, 2.5), ylim = c(0, 6), breaks =  seq(0, 4.5, .125) , main=NULL)
box(bty="l")
lines(density(post_sims_d_plot$steps),col="black",lwd=2)
abline(v = mean(post_sims_d_plot$steps), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")	
hist(post_sims_su_plot$steps, prob=TRUE,col=rgb(0.8,0.8,0.8,0.6), border="white", xlab="Step length (km)",
	ylab = "Density", xlim = c(0, 2.5), ylim = c(0,6), breaks =  seq(0, 4.5, .125), main=NULL, add = T)
#box(bty="l")
lines(density(post_sims_su_plot$steps),col="black",lwd=2)
abline(v = mean(post_sims_su_plot$steps), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")	
	

par(yaxs="i",las=1)
hist(post_sims_d_plot$turns_deg, prob=TRUE,col="grey60",border="white", xlab="Turn angle (deg)",
	ylab = "Density", xlim = c(-180, 180), ylim = c(0, 0.004), breaks = 25, main=NULL)
#main="Distribution of Respirable Particle Concentrations")
box(bty="l")
lines(density(post_sims_d_plot$turns_deg),col="black",lwd=2)
abline(v = mean(post_sims_d_plot$turns_deg), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")
hist(post_sims_su_plot$turns_deg, prob=TRUE,col=rgb(0.8,0.8,0.8,0.6),border="white", xlab="Turn angle (deg)",
	ylab = "Density", xlim = c(-180, 180), breaks = 25, main=NULL, add = T)
#main="Distribution of Respirable Particle Concentrations")
box(bty="l")
lines(density(post_sims_su_plot$turns_deg),col="black",lwd=2)
abline(v = mean(post_sims_su_plot$turns_deg), lty = 2)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")




turns_cat <- post_sims_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) %>%
	mutate(turn_45 =  ifelse(abs_turns_deg <= 45.1, 1, 0)) %>%
	mutate(turn_90 = ifelse(abs_turns_deg <= 90.1 & abs_turns_deg > 45.1, 2, 0)) %>%
	mutate(turn_135 = ifelse(abs_turns_deg <= 135.1 & abs_turns_deg > 90.1, 3, 0)) %>%
	mutate(turn_180 = ifelse(abs_turns_deg > 135.1, 4, 0)) %>%
	mutate(turn_cat = rowSums(.[8:11]))

dive_turns_cat <- post_sims_d_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) %>%
	mutate(turn_45 =  ifelse(abs_turns_deg <= 45.1, 1, 0)) %>%
	mutate(turn_90 = ifelse(abs_turns_deg <= 90.1 & abs_turns_deg > 45.1, 2, 0)) %>%
	mutate(turn_135 = ifelse(abs_turns_deg <= 135.1 & abs_turns_deg > 90.1, 3, 0)) %>%
	mutate(turn_180 = ifelse(abs_turns_deg > 135.1, 4, 0)) %>%
	mutate(turn_cat = rowSums(.[8:11]))

surf_turns_cat <- post_sims_su_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) %>%
	mutate(turn_45 =  ifelse(abs_turns_deg <= 45.1, 1, 0)) %>%
	mutate(turn_90 = ifelse(abs_turns_deg <= 90.1 & abs_turns_deg > 45.1, 2, 0)) %>%
	mutate(turn_135 = ifelse(abs_turns_deg <= 135.1 & abs_turns_deg > 90.1, 3, 0)) %>%
	mutate(turn_180 = ifelse(abs_turns_deg > 135.1, 4, 0)) %>%
	mutate(turn_cat = rowSums(.[8:11]))

near_turns_cat <- post_sims_n_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) %>%
	mutate(turn_45 =  ifelse(abs_turns_deg <= 45.1, 1, 0)) %>%
	mutate(turn_90 = ifelse(abs_turns_deg <= 90.1 & abs_turns_deg > 45.1, 2, 0)) %>%
	mutate(turn_135 = ifelse(abs_turns_deg <= 135.1 & abs_turns_deg > 90.1, 3, 0)) %>%
	mutate(turn_180 = ifelse(abs_turns_deg > 135.1, 4, 0)) %>%
	mutate(turn_cat = rowSums(.[8:11]))

far_turns_cat <- post_sims_f_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) %>%
	mutate(turn_45 =  ifelse(abs_turns_deg <= 45.1, 1, 0)) %>%
	mutate(turn_90 = ifelse(abs_turns_deg <= 90.1 & abs_turns_deg > 45.1, 2, 0)) %>%
	mutate(turn_135 = ifelse(abs_turns_deg <= 135.1 & abs_turns_deg > 90.1, 3, 0)) %>%
	mutate(turn_180 = ifelse(abs_turns_deg > 135.1, 4, 0)) %>%
	mutate(turn_cat = rowSums(.[8:11]))

side_turns_cat <- post_sims_si_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) %>%
	mutate(turn_45 =  ifelse(abs_turns_deg <= 45.1, 1, 0)) %>%
	mutate(turn_90 = ifelse(abs_turns_deg <= 90.1 & abs_turns_deg > 45.1, 2, 0)) %>%
	mutate(turn_135 = ifelse(abs_turns_deg <= 135.1 & abs_turns_deg > 90.1, 3, 0)) %>%
	mutate(turn_180 = ifelse(abs_turns_deg > 135.1, 4, 0)) %>%
	mutate(turn_cat = rowSums(.[8:11]))

front_turns_cat <- post_sims_fr_plot %>%
	mutate(abs_turns_deg = abs(turns_deg)) %>%
	mutate(turn_45 =  ifelse(abs_turns_deg <= 45.1, 1, 0)) %>%
	mutate(turn_90 = ifelse(abs_turns_deg <= 90.1 & abs_turns_deg > 45.1, 2, 0)) %>%
	mutate(turn_135 = ifelse(abs_turns_deg <= 135.1 & abs_turns_deg > 90.1, 3, 0)) %>%
	mutate(turn_180 = ifelse(abs_turns_deg > 135.1, 4, 0)) %>%
	mutate(turn_cat = rowSums(.[8:11]))
	



turns_cat_new <-  turns_cat %>%
	dplyr::select(turn_cat) %>%
	mutate(type = "All")
far_turns_cat_new <-  far_turns_cat %>%
	dplyr::select(turn_cat) %>%
	mutate(type = "Far")
near_turns_cat_new <-  near_turns_cat %>%
	dplyr::select(turn_cat) %>%
	mutate(type = "Near")
surf_turns_cat_new <-  surf_turns_cat %>%
	dplyr::select(turn_cat) %>%
	mutate(type = "Short")
dive_turns_cat_new <-  dive_turns_cat %>%
	dplyr::select(turn_cat) %>%
	mutate(type = "Long")
front_turns_cat_new <- front_turns_cat %>%
	dplyr::select(turn_cat) %>%
	mutate(type = "Front")
side_turns_cat_new <- side_turns_cat %>%
	dplyr::select(turn_cat) %>%
	mutate(type = "Side")
	
counts_tmp <- rbind(turns_cat_new, far_turns_cat_new, near_turns_cat_new,
	surf_turns_cat_new, dive_turns_cat_new, front_turns_cat_new, side_turns_cat_new)
counts_tmp$type <- as.factor(counts_tmp$type)
levels(counts_tmp$type) <- c("All", "Short", "Long", "Near", "Far", "Front", "Side")

counts <- table(counts_tmp$turn_cat, counts_tmp$type)

par(xpd = T, mar = par()$mar + c(0,0,0,7))
barplot(counts, main=NULL, ylab = "Proportion of simulated turns",yaxt="n" ,
  xlab="Data subset category", col=c("grey90","grey70", "grey50", "grey30"))
 legend(9, 15000, c("+/- 45 deg", "+/- 90 deg", "+/- 135 deg","+/- 180 deg"),
	col=c("grey90","grey70", "grey50", "grey30"), pch = 15)
axis(2, at=seq(0, 15000, by=3750),  labels = c("0", "0.25", "0.5", "0.75", "1"))
	par(mar=c(5, 4, 4, 2) + 0.1)
	
	
	
	
	
	
	
	
	
	