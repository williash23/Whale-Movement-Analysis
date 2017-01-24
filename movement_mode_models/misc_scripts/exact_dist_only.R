# Sara Williams
# 12/8/2015; updated 2/1/2016, 3/9/2016
# Data prep for model running script.
################################################################################



################################################################################
# Data prep
#  Load packages
library(dplyr)
library(adehabitatLT)    
library(tidyr)
library(stringr)
library(aspace)

#  Load data
dat_raw <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/Whales_0615_locations_clean.csv")
dat <- dat_raw %>%
            group_by(same_whale_ID) %>%
            filter(n() >1) %>%
            ungroup() %>%
            filter(count == 1) %>%
            filter(year > 2007) %>% ## per info from Karin, PTB (bearing) was kind of weird before 2010
            as.data.frame()

#  Using only data where first observation was a blow 
tmp1_blow <- dat %>%
                         group_by(same_whale_ID) %>%
                         filter(n() > 2) %>%
                         filter(first(whale_behavior)  == "BL-Blowing") %>%
                         filter(first(approx) == "n") %>%
                         ungroup() %>%
                         as.data.frame()
# Use below only if want only ALL observations to be blow or dive.
tmp2_blow <- tmp1_blow %>%
                        filter(whale_behavior  == "DF-Dive-fluke-up" | whale_behavior  == "DN-Dive-no-fluke" |
                        whale_behavior  == "BL-Blowing") %>%
                        group_by(same_whale_ID) %>%
                        filter(n() >1) %>%
                        ungroup()  %>%
                        as.data.frame()
tmp2_blow$whale_behavior <- droplevels(tmp2_blow$whale_behavior)
locs_a_blow <- arrange(tmp2_blow, same_whale_ID, ob_order_time)
locs_b_blow <- locs_a_blow %>%
                          dplyr::select(same_whale_ID, X_whale_UTM, Y_whale_UTM) %>%
                          dplyr::rename(X = X_whale_UTM, Y = Y_whale_UTM)
locs_b_blow$same_whale_ID <- droplevels(locs_b_blow$same_whale_ID)

#   Create ltraj object
whale_traj_blow <- as.ltraj(xy = locs_b_blow[,c("X","Y")], id = locs_b_blow$same_whale_ID, typeII = FALSE)

#   Convert into dataframe
traj_dat_blow <- ld(whale_traj_blow)
traj_dat_blow <- traj_dat_blow %>%
                                     dplyr::select(id, x, y, dist, rel.angle, abs.angle) 
names(traj_dat_blow) <- c("same_whale_ID", "X", "Y", "steps", "turns", "abs_angle")

#   First sighting is blow
obs_1 <- traj_dat_blow %>%
                group_by(same_whale_ID) %>%
                mutate(occ = 1:n()) %>%
                ungroup() %>%
                dplyr::mutate(ID_new = as.numeric(as.factor(as.character(same_whale_ID)))) %>%
                arrange(ID_new) %>%
                as.data.frame()

#   Indexing
npts_1 <- nrow(obs_1)
ind_1 <- obs_1$ID_new
nind_1 <- length(unique(obs_1$ID_new))
nocc_1 <- obs_1 %>%
                  group_by(ID_new) %>%
                  summarise(nocc = n()) %>%
                  .$nocc

#  Data 
l_single <- obs_1$steps/1000
theta_single <- obs_1$turns


#  Make data SQAURE - For double, double cov models
l_double <- obs_1 %>%
                    mutate(steps_km = steps/1000) %>%
                    dplyr::select(ID_new, occ, steps_km) %>%
                    spread(occ, steps_km, fill = NA, convert = FALSE)

theta_double <- obs_1 %>%
                           dplyr::select(ID_new, occ, turns) %>%
                           spread(occ, turns, fill = NA, convert = FALSE)


################################################################################
#  Run model
#  Load packages
library(rjags)
library(mcmcplots)
library(coda)
#  Load "glm" module for JAGS
load.module("glm")
load.module("dic")

#   MCMC settings
nc <- 3
ni <- 30000
nb <- 5000
nt <- 2
na <- 5000

#  Run "single" model
#   Bundle data
jags.dat <- list(npts = npts_1, l = l_single, theta = theta_single)
 
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
single_dic <- dic.samples(out_single, n.iter = 1000, thin = 1, type = "pD")
#  Look at simple MCMC plots
mcmcplot(single_fit)
#  Look at summary output
summary(single_fit)


#  Run "double" model
#   Bundle data
jags.dat <- list(l = l_double, theta = theta_double, nind = nind_1, nocc = nocc_1)

#   Inits function
 inits <- function(){list(v = runif(2, 0.01, 5), 
                                        lambda = c(NA, runif(1, 0.01, 5)), 
                                        eps=runif(1, 0.01, 5),
                                        rho = runif(2, 0.01, 1), 
                                        #mu = runif(2, -3.14159265359, 3.14159265359),
                                        mu = runif(2, -2, 2),
                                        beta0 = runif(1, -5, 5))
                                        }

#   Parameters to monitor
params <- c("v","lambda", "eps", "mu", "rho", "beta0")

#  Run model in rjags.
out_double <- jags.model(data = jags.dat,
                                           file = "C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/models/double_loop.txt", 
                                           inits = inits, 
                                           n.chains = nc, 
                                           n.adapt = na)
update(out_double, n.iter = nb)
double_fit <-coda.samples(out_double,
                                             variable.names= params, 
                                             n.iter = ni, 
                                             thin = nt)
                                             
#  Calculate Rhat
double_rhat <- gelman.diag(double_fit, multivariate = F)[[1]]
#  Calcualte DIC
double_dic <- dic.samples(out_double, n.iter = 1000, thin = 1, type = "pD")

mcmcplot(double_fit)
summary(double_fit)

#################################################################
#  Generate figures
#  Load packages
library(CircStats)
library(reshape)
library(ggplot2)
library(dplyr)
library(rjags)
library(emojifont)


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


#  Density plot
steps_plot <- ggplot(post_sims_plot) + 
                          geom_density(aes(steps), fill = "#EBCC2A", colour = "#EBCC2A") +
                          xlab("Step length (km)") +
                          ylab("Frequency") +
                          xlim(c(0, 5)) +
                          ylim(c(0, 3)) +
                          theme_bw() +
                          theme(panel.grid.minor = element_blank())
steps_plot

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

