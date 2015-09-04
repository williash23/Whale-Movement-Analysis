#  Load packages
	require(stats)
	library(plyr)
	library(dplyr)
	library("R2WinBUGS")
	library("lme4")
	
bd <- "C:/Program Files/WinBUGS14/"
working.directory = getwd()
bugs.directory = bd

#  Read and prep data	
temp <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Behavior-Analysis/data/Whale_Pts.csv")
temp2 <- subset(temp, temp$whale_behavior=="DF-Dive-fluke-up" | temp$whale_behavior=="DN-Dive-no-fluke" | 
temp$whale_behavior=="LF-Lunge-feed" | temp$whale_behavior=="RE-Resting" | temp$whale_behavior=="SA-Surface-active")
temp2$whale_behavior <- droplevels(temp2$whale_behavior)

#   Only use observations where there is more than one observation for comparison of behaviors
temp3 <- temp2[which(temp2$ObType=="MultiOb"),]
#   Only use first sigthing observations where 
temp4 <- temp3[which(temp3$ObOrder_Time==1),]
#   Make a new whale_behavior category that is numeric and attach to dataframe
#   0: transit type behaviors (blowing/ ive with no fluke, dive with fluke up) and 1: stationary type behaviors (lunge feed, resting, surface active)
new_beh <- revalue(temp4$whale_behavior, c("BL-Blowing"=0, "DF-Dive-fluke-up"=0, 
						 "DN-Dive-no-fluke"=0, "LF-Lunge-feed"=1, "RE-Resting"=1, "SA-Surface-active"=1))
y <- cbind(new_beh)
y <- arrange(x, SwB_Wpt_ID, ObOrder_Time)
 
transit <- y[which(y$new_beh==0),]
station <- y[which(y$new_beh==1),]

C.t <- nrow(transit)
C.s <- nrow(station) 
C <- c(C.t, C.s)
N <- nrow(y)
move <- new_beh
 
#  This model uses transit mode as baseline.
sink("Binomial.t.test.txt")
cat("
model {
 
#  Priors
     alpha[i] ~ dnorm(0, 0.01)
     beta[i] ~ dnorm(0, 0.01)

#  Likelihood
  for (i in 1:n) {
     C[i] ~ dbin(phi[i], N) 
     logit(phi[i]) <- alpha + beta * move[i]
  }
 
# Derived quantities
  occ.transit <- exp(alpha) / (1 + exp(alpha))
  occ.station <- exp(alpha + beta) / (1 + exp(alpha + beta))
  occ.diff <- occ.transit - occ.station
  
}
  ",fill=TRUE)
sink()
 
#  Bundle data
 win.data <- list(C = C, N = 526, move = move, n = length(C))
 
#  Inits function
 inits <- function(){ list(alpha=rlnorm(1), beta=rlnorm(1))}
 
#  Parameters to estimate
 params <- c("alpha", "beta", "occ.transit", " occ.station", "occ.diff")
 
#  MCMC settings
 nc <- 3
 ni <- 1200
 nb <- 200
 nt <- 2
 
#  Start Gibbs sampling
 out <- bugs(data=win.data, inits=inits, parameters.to.save=params, model.file="Binomial.t.test.txt", 
						n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, debug = TRUE)
