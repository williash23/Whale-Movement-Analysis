################################################################################

#  Load packages
library(rjags)
library(mcmcplots)
library(coda)
#  Load "glm" module for JAGS
load.module("glm")
################################################################################

#   MCMC settings
nc <- 3
ni <- 40000
nb <- 10000
nt <- 2
na <- 5000
################################################################################
#  Run "single" model
#   Bundle data
jags.dat <- list(npts = npts_1, l = l_sing, theta = theta_sing)
 
#   Inits function
inits <- function(){list(v = runif(1, 0.01,  5), 
                                       lambda = runif(1, 0.01, 5), 
                                       rho = runif(1, 0.01, 1), 
                                       mu = runif(1, -3.14159265359, 3.14159265359))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho")

#  Initialize model and go through adaptation 
out_single <- jags.model(data = jags.dat,
                                          file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                                          inits = inits, 
                                          n.chains = nc, 
                                          n.adapt = na)

#  Burnin
update(out_single, n.iter = nb)

#  Sample posterior
single_fit <- coda.samples(out_single,
                                           variable.names= params, 
                                           n.iter = ni, 
                                           thin = nt)

#  Calculate Rhat
single_rhat <- gelman.diag(single_fit, multivariate = F)[[1]]
#  Calcualte DIC
single_dic <- dic.samples(out_single, n.iter = 5000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(single_fit)
#  Look at summary output
summary(single_fit)


#  Run "double not loop" model
#   Bundle data
jags.dat <- list(npts = npts_1, l = l_sing, theta = theta_sing)
 
#   Inits function
inits <- function(){list(v = runif(2, 0.01,  5), 
                                       lambda = c(NA, runif(1, 0.01, 5)), 
                                       rho = runif(2, 0.01, 1), 
                                       mu = runif(2, -3.14159265359, 3.14159265359),
                                       beta0 = runif(1, -5, 5))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "state_1")

#  Initialize model and go through adaptation 
out_double2 <- jags.model(data = jags.dat,
                                          file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/double2.txt", 
                                          inits = inits, 
                                          n.chains = nc, 
                                          n.adapt = na)

#  Burnin
update(out_double2, n.iter = nb)

#  Sample posterior
double2_fit <- coda.samples(out_double2,
                                           variable.names= params, 
                                           n.iter = ni, 
                                           thin = nt)

#  Calculate Rhat
double2_rhat <- gelman.diag(double2_fit, multivariate = F)[[1]]
#  Calcualte DIC
double2_dic <- dic.samples(out_double2, n.iter = 5000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(double2_fit)
#  Look at summary output
summary(double2_fit)

#  Run "double cov not loop" model
#   Bundle data
jags.dat <- list(npts = npts_1, l = l_sing, theta = theta_sing, ship = ship)
 
#   Inits function
inits <- function(){list(v = runif(2, 0.01,  5), 
                                       lambda = c(NA, runif(1, 0.01, 5)), 
                                       rho = runif(2, 0.01, 1), 
                                       mu = runif(2, -3.14159265359, 3.14159265359),
                                       beta0 = runif(1, -5, 5),
                                       beta1 = runif(1, -5, 5))
                                       }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "state_1", "beta1")

#  Initialize model and go through adaptation 
out_double2_cov <- jags.model(data = jags.dat,
                                                    file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/double2_cov.txt", 
                                                    inits = inits, 
                                                    n.chains = nc, 
                                                    n.adapt = na)

#  Burnin
update(out_double2_cov, n.iter = nb)

#  Sample posterior
double2_cov_fit <- coda.samples(out_double2_cov,
                                                       variable.names= params, 
                                                       n.iter = ni, 
                                                       thin = nt)

#  Calculate Rhat
double2_cov_rhat <- gelman.diag(double2_cov_fit, multivariate = F)[[1]]
#  Calcualte DIC
double2_cov_dic <- dic.samples(out_double2_cov, n.iter = 5000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(double2_cov_fit)
#  Look at summary output
summary(double2_cov_fit)





















#  Load packages
library(coda)
library(adehabitatLT)    
library(tidyr)
library(circular)
library(jagsUI)
library(reshape)
library(ggplot2)
library(dplyr)
################################################################################

#   MCMC settings
nc <- 3
ni <- 10000
nb <- 2000
nt <- 2
################################################################################

#  Run "single" model
#   Bundle data
jags.dat <- list(npts = npts_1, l = l_sing, theta = theta_sing)
 
#   Inits function
inits <- function(){list(v = runif(1, 0.01,  5), 
                                   lambda = runif(1, 0.01, 5), 
                                   rho = runif(1, 0.01, 1), 
                                   mu = runif(1, -3.14159265359, 3.14159265359))
                                   }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho")

single_fit <- jags(jags.dat,
                            inits,
                            params,
                            model.file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/single.txt", 
                            n.chains = nc, 
                            n.thin = nt,
                            n.iter = ni,
                            n.burnin = nb)
print(single_fit)
################################################################################

#  Run "double" model
#   Bundle data
jags.dat <- list(l = l_1, theta = theta_1, nind = nind_1, nocc = nocc_1)
 
#   Inits function
inits <- function(){list(v = runif(2, 0.01, 5), 
                                   lambda = c(NA, runif(1, 0.01, 5)), 
                                   eps=runif(1, 0.01, 5),
                                   rho = runif(2, 0.01, 1), 
                                   mu = runif(2, -3.14159265359, 3.14159265359),
                                   beta0 = runif(1, -5, 5))
                                   }

#   Parameters to monitor
params <- c("v","lambda", "mu", "rho", "beta0", "eps")

double_fit <- jags(jags.dat,
                            inits,
                            params,
                            model.file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/double_loop.txt", 
                            n.chains = nc, 
                            n.thin = nt,
                            n.iter = ni,
                            n.burnin = nb)
print(double_fit)
################################################################################

#  Single
lambda <- single_fit$sims.list$lambda
v <- single_fit$sims.list$v
mu <- single_fit$sims.list$mu
rho <- single_fit$sims.list$rho

sims_tmp <- as.data.frame(cbind(v, lambda, mu, rho))
sims <- sims_tmp %>%
            mutate(scale = (1/lambda))
steps <- numeric(length = nrow(sims))
turns <- numeric(length = nrow(sims))
          
for(i in 1:nrow(sims)){
            steps[i] <- rweibull(1, sims[i,1], sims[i,5])
            turns[i] <- rwrappedcauchy(1, sims[i,3],  sims[i,4])
            }
post_iters <- as.data.frame(cbind(steps, turns))
post_iters$turns[post_iters$turns>pi]=post_iters$turns[post_iters$turns>pi]-2*pi
post_iters_plot <- post_iters %>% 
                                      mutate(iter_num = 1:nrow(post_iters))


#  Double
lambda1 <- double_fit$sims.list$lambda[,1]
v1 <- double_fit$sims.list$v[,1]
mu1 <- double_fit$sims.list$mu[,1]
rho1 <- double_fit$sims.list$rho[,1]
lambda2 <- double_fit$sims.list$lambda[,2]
v2 <- double_fit$sims.list$v[,2]
mu2 <- double_fit$sims.list$mu[,2]
rho2 <- double_fit$sims.list$rho[,2]


sims_tmp1 <- as.data.frame(cbind(v1, lambda1, mu1, rho1))
sims1 <- sims_tmp1 %>%
            mutate(scale = (1/lambda1))
steps1 <- numeric(length = nrow(sims1))
turns1 <- numeric(length = nrow(sims1))
          
for(i in 1:nrow(sims1)){
            steps1[i] <- rweibull(1, sims1[i,1], sims1[i,5])
            turns1[i] <- rwrappedcauchy(1, sims1[i,3],  sims1[i,4])
            }
post_iters1 <- as.data.frame(cbind(steps1, turns1))
post_iters1$turns1[post_iters1$turns1>pi]=post_iters1$turns1[post_iters1$turns1>pi]-2*pi
post_iters_plot1 <- post_iters1 %>% 
                                      mutate(iter_num = 1:nrow(post_iters1))

sims_tmp2 <- as.data.frame(cbind(v2, lambda2, mu2, rho2))
sims2 <- sims_tmp2 %>%
            mutate(scale = (2/lambda2))
steps2 <- numeric(length = nrow(sims2))
turns2 <- numeric(length = nrow(sims2))
          
for(i in 1:nrow(sims2)){
            steps2[i] <- rweibull(1, sims2[i,1], sims1[i,5])
            turns2[i] <- rwrappedcauchy(1, sims2[i,3],  sims2[i,4])
            }
post_iters2 <- as.data.frame(cbind(steps2, turns2))
post_iters2$turns2[post_iters2$turns2>pi]=post_iters2$turns2[post_iters2$turns2>pi]-2*pi
post_iters_plot2 <- post_iters2 %>% 
                                      mutate(iter_num = 1:nrow(post_iters2))

#  Density plot
steps_plot1 <- ggplot(post_iters_plot1) + 
                      geom_density(aes(steps1), fill = "light grey", colour = "grey") +
                      xlab("Step length (km)") +
                      ylab("Frequency") +
                      xlim(c(0, 8)) +
                      ylim(c(0, 2)) +
                      theme_bw() +
                      theme(panel.grid.minor = element_blank())
steps_plot1
#  Point plot
steps_pt1 <- ggplot(post_iters_plot1, aes(steps1, iter_num)) +
                   geom_point(colour = "black", alpha = 0.3) +
                   xlab("Step length (km)") +
                   ylab("Iteration number") +
                   theme_bw()
steps_pt1


#  Density plot
steps_plot2 <- ggplot(post_iters_plot2) + 
                      geom_density(aes(steps2), fill = "light grey", colour = "grey") +
                      xlab("Step length (km)") +
                      ylab("Frequency") +
                      xlim(c(0, 4)) +
                      ylim(c(0, 2)) +
                      theme_bw() +
                      theme(panel.grid.minor = element_blank())
steps_plot12
#  Point plot
steps_pt2 <- ggplot(post_iters_plot2, aes(steps2, iter_num)) +
                   geom_point(colour = "black", alpha = 0.3) +
                   xlab("Step length (km)") +
                   ylab("Iteration number") +
                   theme_bw()
steps_pt2

#  Density plot
turns_plot1 <- ggplot(post_iters_plot1) + 
                      geom_density(aes(turns1), fill = "light grey", colour = "grey") +
                      xlab("Turn angle (rad)") +
                      ylab("Frequency") +
                      xlim(c(-4, 4)) +
                      theme_bw() +
                      theme(panel.grid.minor = element_blank())
turns_plot1
#  Point plot
turns_pt1 <- ggplot(post_iters_plot1, aes(turns1, iter_num)) +
                   geom_point(colour = "black", alpha = 0.3) +
                   xlab("Turn angle (rad)") +
                   ylab("Iteration number") +
                   theme_bw()
turns_pt1


#  Density plot
turns_plot2 <- ggplot(post_iters_plot2) + 
                      geom_density(aes(turns2), fill = "light grey", colour = "grey") +
                      xlab("Turn angle (rad)") +
                      ylab("Frequency") +
                      xlim(c(-4, 4)) +
                      theme_bw() +
                      theme(panel.grid.minor = element_blank())
turns_plot2
#  Point plot
turns_pt2 <- ggplot(post_iters_plot2, aes(turns2, iter_num)) +
                   geom_point(colour = "black", alpha = 0.3) +
                   xlab("Turn angle (rad)") +
                   ylab("Iteration number") +
                   theme_bw()
turns_pt2


            
rwcauchy <- function(n, mu = 0, rho = 0) {
    u = runif(n)
    V = cos(2 * pi * u)
    c = 2 * rho/(1 + rho^2)
    t <- (sign(runif(n) - 0.5) * acos((V + c)/(1 + c * V)) + mu)%%(2 * pi)
    return(t)
}

n <- length(l_sing)
n.sims <- single_fit$mcmc.info$n.samples  # number of simulations
pps.single <- matrix(NA, n.sims, n)  # matrix to save steps
ppt.single <- matrix(NA, n.sims, n)  # matrix to save turns         
            
for (i in 1:n.sims) {
    pps.single[i, ] <- rweibull(n, shape = b[i], scale = 1/a[i])  # in R we use 1/a for the scale of the Weibull because BUGS uses a different  fomrulation
    ppt.single[i, ] <- rwcauchy(n, mu[i], rho[i])
}


ppac.s = matrix(NA, n.sims, 61)
for (i in 1:n.sims) {
    ppac.s[i, ] = acf(log(pps.single[i, ]), lag.max = 60, plot = F)$acf
}

oac = acf(log(l_sing), lag.max = 60, plot = F)  #  observed acf

# matplot(oac$lag,t(ppac.s),lty=1,pch=1,col='gray')
plot(oac$lag, oac$acf, type = "b", lwd = 2, col = 2, xlab = "lag", ylab = "autocorrelation")
miquantile <- function(x) quantile(x, prob = c(0.025, 0.975))
qs = apply(ppac.s, 2, miquantile)
lines(oac$lag, qs[1, ])
lines(oac$lag, qs[2, ])