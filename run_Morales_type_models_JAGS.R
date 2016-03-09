# Sara Williams
# 12/8/2015; updated 2/1/2016
# Model variety of CRWs after code from Morales et al. 2004.
#   Run in JAGS.
################################################################################

#  Load packages
library(R2jags)
library(dplyr)
library(adehabitatLT)    
library(mcmcplots)

#  Load and prep data
dat <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/Whales_0615_locations_clean.csv")

#  Generate step lengths and turning angles using ADEpackage
locs2_a <- arrange(dat, same_whale_ID, ob_order_time)
locs3_a<- dplyr::select(locs2_a, same_whale_ID, X_whale_UTM, Y_whale_UTM)
n_uni_a <- length(unique(locs3_a$same_whale_ID))
locs_tmp_a <- locs3_a %>% 
            distinct(same_whale_ID) %>%
            mutate(Name = seq(1,n_uni_a))
locs_tmp2_a <- full_join(locs3_a, locs_tmp_a, by = "same_whale_ID")        
locs_a <- locs_tmp2_a %>%
        dplyr::select(Name, X_whale_UTM.x, Y_whale_UTM.x) %>%
        dplyr::rename(X = X_whale_UTM.x, Y = Y_whale_UTM.x)
#   Create ltraj object        
whale_traj <- as.ltraj(xy = locs_a[,c("X","Y")], id = locs_a$Name, typeII = FALSE)        
#   Convert into dataframe
traj_dat <- ld(whale_traj)        
#   Connect traj_dat to dataframe 
traj_dat_full_tmp <- cbind(traj_dat, dat)

#  Add location where whale was spotted (as grid ID) to data set
gridID <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Sighting-Density/data/adjusted_density_by_first_sighting.csv")
gridcovs <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Sighting-Density/data/env_ship_covs_sighting_density_by_gridID.csv")
grid_dat <- left_join(gridID, gridcovs, by = "grid_ID")

#  Join covariates from grid cell to whale observations, by the grid that the whale was first sighted in.
traj_dat_full <-  left_join(traj_dat_full_tmp, grid_dat, by = "same_whale_ID")
################################################################################

#  Sort only ID's that have 1+ observations. Used only "single", "double", 
#   "double covariate" models.

#  Select only needed variables and rename
tmp <- traj_dat_full %>%
         dplyr::select(id, x, y, dist, rel.angle, grid_ID, whale_dist_shore_m, ship_whale_dist, ship_speed_scaled,
                                    trk_length_sum_km, sst_clim, chlor_clim, bath, bath_buff_500, cell_area) 
names(tmp) <- c("ID", "X", "Y", "steps", "turns", "gridID", "shore_dist", "ship_dist", "ship_speed_scaled", 
                                    "trk_length_sum_km", "sst_clim", "chlor_clim", "bath", "bath_buff_500", "cell_area")
tmp2 <- filter(tmp,!is.na(steps))
tmp3 <- filter(tmp2,!is.na(turns), !is.na(gridID))
tmp4 <- filter(tmp3, steps < 5000)
tmp5 <- tmp4 %>%
                group_by(ID) %>%
                mutate(same_ID_indicator = ifelse(row_number() == 1, 1,0)) %>%
                as.data.frame()

#  Generate ordered ID variable
raw_ID <- as.numeric(tmp4$ID)
ID <- as.factor(raw_ID)
ID_new <- as.numeric(ID)
tmp6 <- cbind(ID_new, tmp5)
obs <- tmp6 %>%
             arrange(ID_new) %>%
             as.data.frame()

#   Indexing
npts_1 <- nrow(obs)
ID_1 <- obs$ID_new
same_1 <- obs$same_ID_indicator
nind_1 <- length(unique(ID))
#   Steps and turns
l_1 <-as.numeric((obs$steps)/1000)
theta_1 <- as.numeric(obs$turns)
#   Covariates
shore_dist_1 <-as.numeric(scale(obs$shore_dist))
ship_dist_1 <- as.numeric(scale(obs$ship_dist))
ship_speed_1 <-as.numeric(scale(obs$ship_speed_scaled))
ship_dens_1 <- as.numeric(scale(obs$trk_length_sum_km))
chlor_1 <- as.numeric(scale(obs$chlor_clim)) 
sst_1 <- as.numeric(scale(obs$sst_clim)) 
bath_1 <- as.numeric(scale(obs$bath)) 
bath_ave_1 <- as.numeric(scale(obs$bath_buff_500)) 
################################################################################

#  Sort only ID's that have 2+ observations. Used in "double switch" and 
#   "double switch with covariate" models.
tmp7 <- tmp4 %>%
                group_by(ID) %>%
                filter(n() > 1) %>%
                mutate(same_ID_indicator = ifelse(row_number() == 1, 1,0)) %>%
                as.data.frame()

#  Generate ordered ID variable
raw_ID <- as.numeric(tmp7$ID)
ID <- as.factor(raw_ID)
ID_new <- as.numeric(ID)
tmp8 <- cbind(ID_new, tmp7)
obs <- tmp8 %>%
             arrange(ID_new) %>%
             as.data.frame()

 #   Indexing
npts <- nrow(obs)
ID <- obs$ID_new
same <- obs$same_ID_indicator
nind <- length(unique(ID))
#   Steps and turns
l <- as.numeric((obs$steps)/1000)
theta <- as.numeric(obs$turns)
#   Covariates
shore_dist <-as.numeric(scale(obs$shore_dist))
ship_dist <- as.numeric(scale(obs$ship_dist))
ship_speed <-as.numeric(scale(obs$ship_speed_scaled))
ship_dens <- as.numeric(scale(obs$trk_length_sum_km))
chlor <- as.numeric(scale(obs$chlor_clim)) 
sst <- as.numeric(scale(obs$sst_clim)) 
bath <- as.numeric(scale(obs$bath)) 
bath_ave <- as.numeric(scale(obs$bath_buff_500)) 
################################################################################

#   MCMC settings
nc <- 3
ni <- 1000
nb <- 200
nt <- 2
na <- 200
load.module("glm")

#  Run "single" model
#   Bundle data
nstate <- 2
jags.dat <- list(npts = npts_1, l = l_1, theta = theta_1, nind = nind_1)

#   Inits function
inits <- function(){list(v0 = runif(1, 0.01,  5), 
                                             lambda0 = runif(1, 0.01, 5), 
                                             rho0 = runif(1, 0.01, 1), 
                                             mu0 = runif(1, -pi, pi))
                                             }

#   Parameters to monitor
params <- c("mean.v","mean.lambda", "mean.mu", "mean.rho", "scale")

out_single <- jags.model(data = jags.dat,
                                     file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models/single_test.txt", 
                                     inits = inits, 
                                     n.chains = nc, 
                                     n.adapt = na)

update(out_single, n.iter = nb)

single_fit <- coda.samples(out_single,
                                   variable.names= params, 
                                   n.iter = ni, 
                                   thin = nt)

mcmcplot(single_fit)
summary(single_fit)


#  Run "double" model
#   Bundle data
nstate <- 2
jags.dat <- list(npts = npts_1, l = l_1, theta = thet_1a, ID = ID_1, nind = nind_1)

#   Inits function
inits <- function(){list(v = runif(2, 0.01,  5), 
                                             lambda = c(NA, runif(1, 0.01, 5)), 
                                             eps=runif(1, 0.01, 5),
                                             rho = runif(2, 0.01, 1), 
                                             mu = runif(2, -pi, pi),
                                             beta0 = runif(1, -5, 5))
                                             }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "scale", "beta0", "mean.alpha", "prob1", "prob2")

out_double<- jags.model(data = jags.dat,
                                     file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models/double_test.txt", 
                                     inits = inits, 
                                     n.chains = nc, 
                                     n.adapt = na)

update(out_double, n.iter = nb)

double_fit <- coda.samples(out_double,
                                                         variable.names= params, 
                                                         n.iter = ni, 
                                                         thin = nt)

mcmcplot(double_fit)
summary(double_fit)





#  Run "double covariate" model
#   Bundle data
jags.dat <- list(npts = npts_1, l = l_1, theta = theta_1, ID = ID_1, nind = nind_1, sst = sst_1)

#   Inits function
inits <- function(){list(v = runif(2, 0.01,  5), 
                                              lambda = c(NA, runif(1, 0.01, 5)), 
                                              eps=runif(1, 0.01, 5),
                                              rho = runif(2, 0.01, 1), 
                                              mu = runif(2, -pi, pi),
                                              beta0 = runif(1, -5, 5),
                                              beta1 = runif(1, -5, 5))
                                              }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "scale", "beta0", "beta1")

out_double_cov <- jags.model(data = jags.dat,
                                                                file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models/double_cov_test.txt", 
                                                                inits = inits, 
                                                                n.chains = nc, 
                                                                n.adapt = na)

update(out_double_cov, n.iter = nb)

double_cov_fit <- coda.samples(out_double_cov,
                                                                   variable.names= params, 
                                                                   n.iter = ni, 
                                                                   thin = nt)

mcmcplot(double_cov_fit)
summary(double_cov_fit)




#  Run "double switch" model
#   Bundle data
nstate <- 2
jags.dat <- list(npts = npts, l = l, theta = theta, ID = ID, nstate = nstate, nind = nind)

#   Inits function
inits <- function(){list(v = runif(2, 0.01,  5), 
                                             lambda = c(NA, runif(1, 0.01, 5)), 
                                             eps=runif(1, 0.01, 5),
                                             rho = runif(2, 0.01, 1), 
                                             mu = runif(2, -pi, pi),
                                             phi = c(runif(1, 0.01, 1), NA),
                                             beta0 = runif(2, -5, 5))
                                             }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "scale", "beta0")

out_double_sw<- jags.model(data = jags.dat,
                                     file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models/double_switch_test.txt", 
                                     inits = inits, 
                                     n.chains = nc, 
                                     n.adapt = na)

update(out_double_sw, n.iter = nb)

double_sw_fit <- coda.samples(out_double_sw,
                                    variable.names= params, 
                                    n.iter = ni, 
                                    thin = nt)

mcmcplot(double_sw_fit)
summary(double_sw_fit)



#  Run "double switch covariate" model
#   Bundle data
nstate <- 2
jags.dat <- list(npts = npts, l = l, theta = theta, ID = ID, nstate = nstate, nind = nind, shipdens = ship_dens)

#   Inits function
inits <- function(){list(v = runif(2, 0.01,  5), 
                                             lambda = c(NA, runif(1, 0.01, 5)), 
                                             eps=runif(1, 0.01, 5),
                                             rho = runif(2, 0.01, 1), 
                                             mu = runif(2, -pi, pi),
                                             phi = c(runif(1, 0.01, 1), NA),
                                             beta0 = runif(2, -5, 5),
                                             beta1 = runif(2, -5, 5))
                                             }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "scale", "beta0", "beta1")

out_double_sw_cov<- jags.model(data = jags.dat,
                                                                       file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/models/double_switch_cov_test.txt", 
                                                                       inits = inits, 
                                                                       n.chains = nc, 
                                                                       n.adapt = na)

update(out_double_sw_cov, n.iter = nb)

double_sw_cov_fit <- coda.samples(out_double_sw_cov,
                                                                           variable.names= params, 
                                                                           n.iter = ni, 
                                                                           thin = nt)

mcmcplot(double_sw_cov_fit)
summary(double_sw_cov_fit)














#  Hierarchical switch model
model{

# Hyperpriors. These are priors on the parameters of the prior distributions for individual
# level parameters.

# Amu.h and Atau.h are priors for the parameters of the Normal distribution used as prior
# for scale parameter in step length during fast movement. eps.h and epstau.h are priors 
# for the censored Normal used for the difference between scale parameter of step length
# in fast and slow movement. 

Amu.h ~  dnorm(0.0,0.01)I(0,)
Atau.h ~ dnorm(0.0,0.01)I(0,)
eps.h ~ dnorm(0,0.01)I(0,)
epstau.h ~ dnorm(0,0.01)I(0,)

# priors on the parameters of the Normal distributions used as priors on the shape 
# parameters of the Weibul distributions used to model step length
Bmu.h[1] ~ dnorm(0.0,0.01)I(0,)
Bmu.h[2] ~ dnorm(0.0,0.01)I(0,)
Btau.h[1] ~ dnorm(0.0,0.01)I(0,)
Btau.h[2] ~ dnorm(0.0,0.01)I(0,)

# priors on the parameters of the Normal distributions used as priors on the mean
# direction of turning angles
mumean.h[1] ~ dunif(0,6.28318530717959)
mumean.h[2] ~ dunif(0,6.28318530717959)
mutau.h[1] ~ dnorm(0,0.01)I(0.0,)
mutau.h[2] ~ dnorm(0,0.01)I(0.0,)

# priors on the parameters of the Beta distribution used as priors on the mean 
# cosine of turning angles
arho.h[1] ~ dnorm(0,0.001)I(0,)
arho.h[2] ~ dnorm(0,0.001)I(0,)
brho.h[1] ~ dnorm(0,0.001)I(0,)
brho.h[2] ~ dnorm(0,0.001)I(0,)

# priors on the parameters of the Beta distribution used as priors on switching rates
qa.h[1] ~ dnorm(0,0.01)I(0,)
qa.h[2] ~ dnorm(0,0.01)I(0,)
qb.h[1] ~ dnorm(0,0.01)I(0,)
qb.h[2] ~ dnorm(0,0.01)I(0,)

Pi <- 3.14159265359   # define pi
# assign initial state
for(i in 1:npaths){
idx[1,i]~dcat(phi[])
}

# iterate over movement paths
  for(k in 1:npaths){
    # individual level priors
    A[k,2] ~ dnorm(Amu.h, Asigma.h)
    eps[k] ~ dnorm(eps.h, epstau.h)I(0,)
    A[k,1] <- A[k,2] + eps[k]
    B[k,1] ~ dnorm(Bmu.h[1], Btau.h[1]) I(0,)
    B[k,2] ~ dnrom(Bmu.h[2], Btau.h[2]) I(0,)
    RHO[k,1] ~ dbeta(arho.h[1], brho.h[1])
    RHO[k,2] ~ dbeta(arho.h[2], brho.h[2])
    MU[k,1] ~ dnorm(mumean.h[1], mutau.h[1])
    MU[k,2] ~ dnorm(mumean.h[2], mutau.h[2])
    Q[k,1] ~ dbeta(qa.h[1],qb.h[1])
    Q[k,2] ~ dbeta(qa.h[2],qb.h[2])

# iterate over observations
    for (t in 2:npts) {
      l[t,k] ~ dweib(B[k,idx[t,k]], A[k,idx[t,k]]) 

# use 'ones'trick to simulate the wrapped Cauchy distribution
      ones[t,k] <- 1
      ones[t,k] ~ dbern( wc[t,k] )
      wc[t,k] <- (1/(2*Pi)*(1-rho[t,k]*rho[t,k])/(1+rho[t,k]*rho[t,k]-2*rho[t,k]*cos(theta[t,k]-mu.t[t,k])))/ 300

      theta[t,k] ~ dunif(0,6.28318530717959)
      rho[t,k]<-RHO[k,idx[t,k]]
      mu.t[t,k]<- MU[k,idx[t,k]]
      # the probability of being in movement type 1
      p[t,k,1] <- Q[k,idx[t-1,k]]
      p[t,k,2] <- 1-Q[k,idx[t-1,k]]
      idx[t,k] ~ dcat(p[t,k,])
    }
  }
}










####################################################################################################

#  JAGS parameterization of weibull distribution.
#   dweib(v = shape, lambda)
#  R parameterization of weibull distribution.
#   rweibull(n, shape, scale)
#  Transform lambda to R scale parameter.
#    scale[1] <- (1/lambda[1])^ (1/v[1])
#  Find mean step length.
#   mean <- scale*gamma((1+1/shape))

#  Plots.
n <- 1000000

#  Step lengths.

#  Simulate step lengths.
#   Uses parameters estimated from Morales et al. code.
sim_steps1 <- as.data.frame(rweibull(n, 1.063, 600))
sim_steps2 <- as.data.frame(rweibull(n, 1.210,  300))
names(sim_steps1)[1] <- "dist"
names(sim_steps2)[1] <- "dist"
#   Create dataframe of simulated movement 1 and movement 2 steps.
sim_steps <- cbind(sim_steps1, sim_steps2)
names(sim_steps)[1] <- "1"
names(sim_steps)[2] <- "2"
#  Look at density plots of both on same figure.
df_steps_sim <- melt(sim_steps)
compare_steps_sim <- ggplot(df_steps_sim) + 
                         geom_density(aes(x = value, colour = variable, fill = variable), alpha = .25) +
                         xlim(0, 5000) +
                         xlab("Step length (m)") +
                         theme_bw()
 compare_steps_sim


rhoD <- est.rho(turns.dive)
muD<- circ.mean(turns.dive)

sim_turns1 <- rwrpcauchy(n, muD, rhoD)


data{
  #  Generate vector of 1's for 1's trick for wrapped Cauchy distribution
  for(i in 2:npts){
    ones[i] <- 1
    }  
  }
