# Sara Williams
# 2/22/2016; updated 1/19/2017
# Generate figures from parameter estimates from movement state models
#   and run in JAGS. For data seperated into surface interval and deep dive by time between
#   observations.
################################################################################

#  Load packages
library(CircStats)
library(reshape)
library(ggplot2)
library(dplyr)
library(rjags)
#library(emojifont)

################################################################################
# Load output from MCMC posterior distribtuion iterations (JAGS object)
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/double_fit.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/double_cov_fit.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit_near.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit_far.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit_side.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit_front.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit_dive.RData")
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_mode_models/model_output/single_fit_surf.RData")
################################################################################

# NOTE: due to differences in formulation of Weibull distribution between JAGS and R:
###### scale = (1/lambda)^(1/v) ######
niter <- 11000
nsamp <- 8000

rwcauchy <- function(n, mu = 0, rho = 0) {
  u = runif(n)
  V = cos(2 * pi * u)
  c = 2 * rho/(1 + rho^2)
  t <- (sign(runif(n) - 0.5) * acos((V + c)/(1 + c * V)) + mu)%%(2 * pi)
  return(t)
}

#geom_vline(aes(xintercept=mean(weight)),
#color="blue", linetype="dashed", size=1)
################################################################################



#  All data in single movement mode model
#   Select rows from posterior distribution for each parameter estimate.
#   model_fit_object[rows, columns, sheets]
#   Insert model output object as "x"
x <- single_fit
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

for(i in 1:nrow(sims)){
  steps[i] <- rweibull(1, sims[i,3], (1/sims[i,4])^(1/sims[i,3]))
  turns[i] <- rwcauchy(1, sims[i,1], sims[i,2])
  #rwrappedcauchy(1, sims[i,1],  sims[i,2])
 }

post_sims <- as.data.frame(cbind(steps, turns))
post_sims$turns[post_sims$turns>pi]=post_sims$turns[post_sims$turns>pi]-2*pi

post_sims_plot <- post_sims %>% 
                              mutate(iter_num = 1:nrow(post_sims))
################################################################################      

#  Density plots
steps_all <- ggplot(post_sims_plot, aes(steps)) + 
                    geom_density(size = 1.25) +
                    geom_vline(aes(xintercept=mean(steps)), color="black", linetype="dashed", size=1) +
                    xlab("\n Step length (km)") +
                    ylab("Frequency \n") +
                    xlim(c(0, 5)) +
                    #ylim(c(0, 2.5)) +
                    theme_bw() +
                    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                    axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
                    #theme(legend.position="none")
steps_all

turns_all <- ggplot(post_sims_plot, aes(turns)) + 
                    geom_density(size = 1.25) +
                    geom_vline(aes(xintercept=mean(turns)), color="black", linetype="dashed", size=1) +
                    xlab("\n Turn angle (km)") +
                    ylab("Frequency \n") +
                    #xlim(c(0, 5)) +
                    #ylim(c(0, 2.5)) +
                    theme_bw() +
                    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                    axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2))) 
                    #theme(legend.position="none")
turns_all
################################################################################      



#  All data - two movement mode
x_all_d <- double_fit
keep_1 <- sample(1:niter, nsamp, replace = F)
keep_2 <- sample(1:niter, nsamp, replace = F)
keep_3 <- sample(1:niter, nsamp, replace = F)

chain_1_all_d1 <- x_all_d[[1]]
sims_1_all_d1 <- chain_1_all_d1[keep_1, c(5, 7, 9, 3)]
chain_1_all_d2<- x_all_d[[1]]
sims_1_all_d2 <- chain_1_all_d2[keep_1, c(6, 8, 10, 4)]

chain_2_all_d1 <- x_all_d[[2]]
sims_2_all_d1 <- chain_2_all_d1[keep_2, c(5, 7, 9, 3)]
chain_2_all_d2<- x_all_d[[2]]
sims_2_all_d2 <- chain_2_all_d2[keep_2, c(6, 8, 10, 4)]

chain_3_all_d1 <- x_all_d[[3]]
sims_3_all_d1 <- chain_3_all_d1[keep_3, c(5, 7, 9, 3)]
chain_3_all_d2<- x_all_d[[3]]
sims_3_all_d2 <- chain_3_all_d2[keep_3, c(6, 8, 10, 4)]

sims_all_d1 <- rbind(sims_1_all_d1, sims_2_all_d1, sims_3_all_d1)
sims_all_d2 <- rbind(sims_1_all_d2, sims_2_all_d2, sims_3_all_d2)

steps_all_d1 <- numeric(length = nrow(sims_all_d1))
turns_all_d1 <- numeric(length = nrow(sims_all_d1))
steps_all_d2 <- numeric(length = nrow(sims_all_d2))
turns_all_d2 <- numeric(length = nrow(sims_all_d2))
         
for(i in 1:nrow(sims_all_d1)){
  steps_all_d1[i] <- rweibull(1, sims_all_d1[i,3], (1/sims_all_d1[i,4])^(1/sims_all_d1[i,3]))
  turns_all_d1[i] <- rwcauchy(1, sims_all_d1[i,1], sims_all_d1[i,2])
  #rwrappedcauchy(1, sims_all_d1[i,1],  sims_all_d1[i,2])
  }
for(i in 1:nrow(sims_all_d2)){
  steps_all_d2[i] <- rweibull(1, sims_all_d2[i,3], (1/sims_all_d2[i,4])^(1/sims_all_d2[i,3]))
  turns_all_d2[i] <- rwcauchy(1, sims_all_d2[i,1], sims_all_d2[i,2])
  #rwrappedcauchy(1, sims_all_d2[i,1],  sims_all_d2[i,2])
  }

post_sims_all_d1 <- as.data.frame(cbind(steps_all_d1, turns_all_d1))
post_sims_all_d1$turns_all_d1[post_sims_all_d1$turns_all_d1>pi]=post_sims_all_d1$turns_all_d1[post_sims_all_d1$turns_all_d1>pi]-2*pi
post_sims_all_d2 <- as.data.frame(cbind(steps_all_d2, turns_all_d2))
post_sims_all_d2$turns_all_d2[post_sims_all_d2$turns_all_d2>pi]=post_sims_all_d2$turns_all_d2[post_sims_all_d2$turns_all_d2>pi]-2*pi

post_sims_all_d1_plot <- post_sims_all_d1 %>% 
                                          mutate(iter_num = 1:nrow(post_sims_all_d1)) %>%
                                          mutate(type = "Two movement modes - 1")  %>%
                                          dplyr::rename(steps = steps_all_d1, turns = turns_all_d1)
                                        
post_sims_all_d2_plot <- post_sims_all_d2 %>% 
                                          mutate(iter_num = 1:nrow(post_sims_all_d2)) %>%
                                          mutate(type = "Two movement modes - 2") %>%
                                          dplyr::rename(steps = steps_all_d2, turns = turns_all_d2) 
                                  
all_double <- rbind(post_sims_all_d1_plot, post_sims_all_d2_plot)
all_double$type <- as.factor(all_double$type)
################################################################################      

#  Density plots
steps_all_double_plot <- ggplot(all_double, aes(steps, linetype = type)) + 
                                         geom_density(size = 1.25) +
                                         xlab("\n Step length (km)") +
                                         ylab("Frequency \n") +
                                         xlim(c(0, 15)) +
                                         #ylim(c(0, 2.5)) +
                                         theme_bw() +
                                         theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                         axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2))) +
                                         scale_linetype_manual(values = c(1, 2))
                                         #theme(legend.position="none")
steps_all_double_plot

turns_all_double_plot <- ggplot(all_double, aes(turns, linetype = type)) + 
                                         geom_density(size = 1.25) +
                                         xlab("\n Turn angle (km)") +
                                         ylab("Frequency \n") +
                                         #xlim(c(0, 5)) +
                                         #ylim(c(0, 2.5)) +
                                         theme_bw() +
                                         theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                         axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2))) +
                                         scale_linetype_manual(values = c(1, 2))
                                         #theme(legend.position="none")
turns_all_double_plot
################################################################################   



#  Deep dive vs surface interval
x_d <- single_fit_dive
x_su <- single_fit_surf
x_all <- single_fit
keep_1 <- sample(1:niter, nsamp, replace = F)
keep_2 <- sample(1:niter, nsamp, replace = F)
keep_3 <- sample(1:niter, nsamp, replace = F)

chain_1_d <- x_d[[1]]
sims_1_d <- chain_1_d[keep_1, c(2, 3, 4, 1)]
chain_1_su<- x_su[[1]]
sims_1_su <- chain_1_su[keep_1, c(2, 3, 4, 1)]
chain_1_all<- x_all[[1]]
sims_1_all <- chain_1_all[keep_1, c(2, 3, 4, 1)]

chain_2_d <- x_d[[2]]
sims_2_d <- chain_2_d[keep_2, c(2, 3, 4, 1)]
chain_2_su <- x_su[[2]]
sims_2_su <- chain_2_su[keep_2, c(2, 3, 4, 1)]
chain_2_all<- x_all[[2]]
sims_2_all <- chain_2_all[keep_2, c(2, 3, 4, 1)]

chain_3_d <- x_d[[3]]
sims_3_d <- chain_3_d[keep_3, c(2, 3, 4, 1)]
chain_3_su <- x_su[[3]]
sims_3_su <- chain_3_su[keep_3, c(2, 3, 4, 1)]
chain_3_all<- x_all[[3]]
sims_3_all <- chain_3_all[keep_3, c(2, 3, 4, 1)]

sims_d <- rbind(sims_1_d, sims_2_d, sims_3_d)
sims_su <- rbind(sims_1_su, sims_2_su, sims_3_su)
sims_all <- rbind(sims_1_all, sims_2_all, sims_3_all)

steps_d <- numeric(length = nrow(sims_d))
turns_d <- numeric(length = nrow(sims_d))
steps_su <- numeric(length = nrow(sims_su))
turns_su <- numeric(length = nrow(sims_su))
steps_all <- numeric(length = nrow(sims_all))
turns_all <- numeric(length = nrow(sims_all))

for(i in 1:nrow(sims_d)){
  steps_d[i] <- rweibull(1, sims_d[i,3], (1/sims_d[i,4])^(1/sims_d[i,3]))
  turns_d[i] <- rwcauchy(1, sims_d[i,1], sims_d[i,2])
  #rwrappedcauchy(1, sims_d[i,1],  sims_d[i,2])
  }
for(i in 1:nrow(sims_su)){
  steps_su[i] <- rweibull(1, sims_su[i,3], (1/sims_su[i,4])^(1/sims_su[i,3]))
  turns_su[i] <- rwcauchy(1, sims_su[i,1], sims_su[i,2])
  #rwrappedcauchy(1, sims_su[i,1],  sims_su[i,2])
  }
for(i in 1:nrow(sims_all)){
  steps_all[i] <- rweibull(1, sims_all[i,3], (1/sims_all[i,4])^(1/sims_all[i,3]))
  turns_all[i] <- rwcauchy(1, sims_all[i,1], sims_all[i,2])
  #rwrappedcauchy(1, sims_all[i,1],  sims_all[i,2])
  }

post_sims_d <- as.data.frame(cbind(steps_d, turns_d))
post_sims_d$turns_d[post_sims_d$turns_d>pi]=post_sims_d$turns_d[post_sims_d$turns_d>pi]-2*pi
post_sims_su <- as.data.frame(cbind(steps_su, turns_su))
post_sims_su$turns_su[post_sims_su$turns_su>pi]=post_sims_su$turns_su[post_sims_su$turns_su>pi]-2*pi
post_sims_all <- as.data.frame(cbind(steps_all, turns_all))
post_sims_all$turns_all[post_sims_all$turns_all>pi]=post_sims_all$turns_all[post_sims_all$turns_all>pi]-2*pi

post_sims_d_plot <- post_sims_d %>% 
                                  mutate(iter_num = 1:nrow(post_sims_d)) %>%
                                  mutate(type = "Deep dive data")  %>%
                                  dplyr::rename(steps = steps_d, turns = turns_d)
post_sims_su_plot <- post_sims_su %>% 
                                   mutate(iter_num = 1:nrow(post_sims_su)) %>%
                                   mutate(type = "Surface interval data") %>%
                                   dplyr::rename(steps = steps_su, turns = turns_su)
post_sims_all_plot <- post_sims_all %>% 
                                    mutate(iter_num = 1:nrow(post_sims_all)) %>%
                                    mutate(type = "All data") %>%
                                    dplyr::rename(steps = steps_all, turns = turns_all)
                                    
all_surf_dive <- rbind(post_sims_d_plot, post_sims_su_plot, post_sims_all_plot)
all_surf_dive$type <- as.factor(all_surf_dive$type)
################################################################################      

#  Density plots
steps_all_surf_dive_plot <- ggplot(all_surf_dive, aes(steps, linetype = type)) + 
                                             geom_density(size = 1.25) +
                                             xlab("\n Step length (km)") +
                                             ylab("Frequency \n") +
                                             xlim(c(0, 5)) +
                                             #ylim(c(0, 2.5)) +
                                             theme_bw() +
                                             theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                             axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2))) +
                                             scale_linetype_manual(values = c(1, 2, 3))
                                             #theme(legend.position="none")
steps_all_surf_dive_plot

turns_all_surf_dive_plot <- ggplot(all_surf_dive, aes(turns, linetype = type)) + 
                                             geom_density(size = 1.25) +
                                             xlab("\n Turn angle (km)") +
                                             ylab("Frequency \n") +
                                             #xlim(c(0, 5)) +
                                             #ylim(c(0, 2.5)) +
                                             theme_bw() +
                                             theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                             axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2))) +
                                             scale_linetype_manual(values = c(1, 2, 3))
                                             #theme(legend.position="none")
turns_all_surf_dive_plot
################################################################################



#  First observation near to ship or far from ship
x_n <- single_fit_near
x_f <- single_fit_far
x_all <- single_fit
keep_1 <- sample(1:niter, nsamp, replace = F)
keep_2 <- sample(1:niter, nsamp, replace = F)
keep_3 <- sample(1:niter, nsamp, replace = F)

chain_1_n <- x_n[[1]]
sims_1_n <- chain_1_n[keep_1, c(2, 3, 4, 1)]
chain_1_f<- x_f[[1]]
sims_1_f <- chain_1_f[keep_1, c(2, 3, 4, 1)]
chain_1_all<- x_all[[1]]
sims_1_all <- chain_1_all[keep_1, c(2, 3, 4, 1)]

chain_2_n <- x_n[[2]]
sims_2_n <- chain_2_n[keep_2, c(2, 3, 4, 1)]
chain_2_f <- x_f[[2]]
sims_2_f <- chain_2_f[keep_2, c(2, 3, 4, 1)]
chain_2_all<- x_all[[2]]
sims_2_all <- chain_2_all[keep_2, c(2, 3, 4, 1)]

chain_3_n <- x_n[[3]]
sims_3_n <- chain_3_n[keep_3, c(2, 3, 4, 1)]
chain_3_f <- x_f[[3]]
sims_3_f <- chain_3_f[keep_3, c(2, 3, 4, 1)]
chain_3_all<- x_all[[3]]
sims_3_all <- chain_3_all[keep_3, c(2, 3, 4, 1)]

sims_n <- rbind(sims_1_n, sims_2_n, sims_3_n)
sims_f <- rbind(sims_1_f, sims_2_f, sims_3_f)
sims_all <- rbind(sims_1_all, sims_2_all, sims_3_all)

steps_n <- numeric(length = nrow(sims_n))
turns_n <- numeric(length = nrow(sims_n))
steps_f <- numeric(length = nrow(sims_f))
turns_f <- numeric(length = nrow(sims_f))
steps_all <- numeric(length = nrow(sims_all))
turns_all <- numeric(length = nrow(sims_all))

for(i in 1:nrow(sims_n)){
  steps_n[i] <- rweibull(1, sims_n[i,3], (1/sims_n[i,4])^(1/sims_n[i,3]))
  turns_n[i] <- rwcauchy(1, sims_n[i,1], sims_n[i,2])
  #rwrappedcauchy(1, sims_n[i,1],  sims_n[i,2])
  }
for(i in 1:nrow(sims_f)){
  steps_f[i] <- rweibull(1, sims_f[i,3], (1/sims_f[i,4])^(1/sims_f[i,3]))
  turns_f[i] <- rwcauchy(1, sims_f[i,1], sims_f[i,2])
  #rwrappedcauchy(1, sims_f[i,1],  sims_f[i,2])
  }
for(i in 1:nrow(sims_all)){
  steps_all[i] <- rweibull(1, sims_all[i,3], (1/sims_all[i,4])^(1/sims_all[i,3]))
  turns_all[i] <- rwcauchy(1, sims_all[i,1], sims_all[i,2])
  #rwrappedcauchy(1, sims_all[i,1],  sims_all[i,2])
  }

post_sims_n <- as.data.frame(cbind(steps_n, turns_n))
post_sims_n$turns_n[post_sims_n$turns_n>pi]=post_sims_n$turns_n[post_sims_n$turns_n>pi]-2*pi
post_sims_f <- as.data.frame(cbind(steps_f, turns_f))
post_sims_f$turns_f[post_sims_f$turns_f>pi]=post_sims_f$turns_f[post_sims_f$turns_f>pi]-2*pi
post_sims_all <- as.data.frame(cbind(steps_all, turns_all))
post_sims_all$turns_all[post_sims_all$turns_all>pi]=post_sims_all$turns_all[post_sims_all$turns_all>pi]-2*pi

post_sims_n_plot <- post_sims_n %>% 
                                  mutate(iter_num = 1:nrow(post_sims_n)) %>%
                                  mutate(type = "Near data")  %>%
                                  dplyr::rename(steps = steps_n, turns = turns_n)
post_sims_f_plot <- post_sims_f %>% 
                                 mutate(iter_num = 1:nrow(post_sims_f)) %>%
                                 mutate(type = "Far data") %>%
                                 dplyr::rename(steps = steps_f, turns = turns_f)
post_sims_all_plot <- post_sims_all %>% 
                                    mutate(iter_num = 1:nrow(post_sims_all)) %>%
                                    mutate(type = "All data") %>%
                                    dplyr::rename(steps = steps_all, turns = turns_all)
                                    
all_near_far <- rbind(post_sims_n_plot, post_sims_f_plot, post_sims_all_plot)
all_near_far$type <- as.factor(all_near_far$type)
################################################################################      

#  Density plots
steps_all_near_far_plot <- ggplot(all_near_far, aes(steps, linetype = type)) + 
                                             geom_density(size = 1.25) +
                                             xlab("\n Step length (km)") +
                                             ylab("Frequency \n") +
                                             xlim(c(0, 5)) +
                                             #ylim(c(0, 2.5)) +
                                             theme_bw() +
                                             theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                             axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2))) +
                                             scale_linetype_manual(values = c(1, 2, 3))
                                             #theme(legend.position="none")
steps_all_near_far_plot

turns_all_near_far_plot <- ggplot(all_near_far, aes(turns, linetype = type)) + 
                                             geom_density(size = 1.25) +
                                             xlab("\n Turn angle (km)") +
                                             ylab("Frequency \n") +
                                             #xlim(c(0, 5)) +
                                             #ylim(c(0, 2.5)) +
                                             theme_bw() +
                                             theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                             axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2))) +
                                             scale_linetype_manual(values = c(1, 2, 3))
                                             #theme(legend.position="none")
turns_all_near_far_plot
################################################################################



#  First observation at side of ship or at front of ship
x_si <- single_fit_side
x_fr <- single_fit_front
x_all <- single_fit
keep_1 <- sample(1:niter, nsamp, replace = F)
keep_2 <- sample(1:niter, nsamp, replace = F)
keep_3 <- sample(1:niter, nsamp, replace = )

chain_1_si <- x_si[[1]]
sims_1_si <- chain_1_si[keep_1, c(2, 3, 4, 1)]
chain_1_fr<- x_fr[[1]]
sims_1_fr <- chain_1_fr[keep_1, c(2, 3, 4, 1)]
chain_1_all<- x_all[[1]]
sims_1_all <- chain_1_all[keep_1, c(2, 3, 4, 1)]

chain_2_si <- x_si[[2]]
sims_2_si <- chain_2_si[keep_2, c(2, 3, 4, 1)]
chain_2_fr <- x_fr[[2]]
sims_2_fr <- chain_2_fr[keep_2, c(2, 3, 4, 1)]
chain_2_all<- x_all[[2]]
sims_2_all <- chain_2_all[keep_2, c(2, 3, 4, 1)]

chain_3_si <- x_si[[3]]
sims_3_si <- chain_3_si[keep_3, c(2, 3, 4, 1)]
chain_3_fr <- x_fr[[3]]
sims_3_fr <- chain_3_fr[keep_3, c(2, 3, 4, 1)]
chain_3_all<- x_all[[3]]
sims_3_all <- chain_3_all[keep_3, c(2, 3, 4, 1)]

sims_si <- rbind(sims_1_si, sims_2_si, sims_3_si)
sims_fr <- rbind(sims_1_fr, sims_2_fr, sims_3_fr)
sims_all <- rbind(sims_1_all, sims_2_all, sims_3_all)

steps_si <- numeric(length = nrow(sims_si))
turns_si <- numeric(length = nrow(sims_si))
steps_fr <- numeric(length = nrow(sims_fr))
turns_fr <- numeric(length = nrow(sims_fr))
steps_all <- numeric(length = nrow(sims_all))
turns_all <- numeric(length = nrow(sims_all))

for(i in 1:nrow(sims_si)){
  steps_si[i] <- rweibull(1, sims_si[i,3], (1/sims_si[i,4])^(1/sims_si[i,3]))
  turns_si[i] <- rwcauchy(1, sims_si[i,1], sims_si[i,2])
  #rwrappedcauchy(1, sims_si[i,1],  sims_si[i,2])
  }
for(i in 1:nrow(sims_fr)){
  steps_fr[i] <- rweibull(1, sims_fr[i,3], (1/sims_fr[i,4])^(1/sims_fr[i,3]))
  turns_fr[i] <- rwcauchy(1, sims_fr[i,1], sims_fr[i,2])
  #rwrappedcauchy(1, sims_fr[i,1],  sims_fr[i,2])
  }
for(i in 1:nrow(sims_all)){
  steps_all[i] <- rweibull(1, sims_all[i,3], (1/sims_all[i,4])^(1/sims_all[i,3]))
  turns_all[i] <- rwcauchy(1, sims_all[i,1], sims_all[i,2])
  #rwrappedcauchy(1, sims_all[i,1],  sims_all[i,2])
  }

post_sims_si <- as.data.frame(cbind(steps_si, turns_si))
post_sims_si$turns_si[post_sims_si$turns_si>pi]=post_sims_si$turns_si[post_sims_si$turns_si>pi]-2*pi
post_sims_fr <- as.data.frame(cbind(steps_fr, turns_fr))
post_sims_fr$turns_fr[post_sims_fr$turns_fr>pi]=post_sims_fr$turns_fr[post_sims_fr$turns_fr>pi]-2*pi
post_sims_all <- as.data.frame(cbind(steps_all, turns_all))
post_sims_all$turns_all[post_sims_all$turns_all>pi]=post_sims_all$turns_all[post_sims_all$turns_all>pi]-2*pi

post_sims_si_plot <- post_sims_si %>% 
                                   mutate(iter_num = 1:nrow(post_sims_si)) %>%
                                   mutate(type = "Side approach data")  %>%
                                   dplyr::rename(steps = steps_si, turns = turns_si)
post_sims_fr_plot <- post_sims_fr %>% 
                                   mutate(iter_num = 1:nrow(post_sims_fr)) %>%
                                   mutate(type = "Front approach data") %>%
                                   dplyr::rename(steps = steps_fr, turns = turns_fr)
post_sims_all_plot <- post_sims_all %>% 
                                    mutate(iter_num = 1:nrow(post_sims_all)) %>%
                                    mutate(type = "All data") %>%
                                    dplyr::rename(steps = steps_all, turns = turns_all)
                                    
all_side_front <- rbind(post_sims_si_plot, post_sims_fr_plot, post_sims_all_plot)
all_side_front$type <- as.factor(all_side_front$type)
################################################################################      

#  Density plots
steps_all_side_front_plot <- ggplot(all_side_front, aes(steps, linetype = type)) + 
                                               geom_density(size = 1.25) +
                                               xlab("\n Step length (km)") +
                                               ylab("Frequency \n") +
                                               #xlim(c(0, 5)) +
                                               #ylim(c(0, 2.5)) +
                                               theme_bw() +
                                               theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                               axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2))) +
                                               scale_linetype_manual(values = c(1, 2, 3))
                                               #theme(legend.position="none")
steps_all_side_front_plot

turns_all_side_front_plot <- ggplot(all_side_front, aes(turns, linetype = type)) + 
                                               geom_density(size = 1.25) +
                                               xlab("\n Turn angle (km)") +
                                               ylab("Frequency \n") +
                                               #xlim(c(0, 5)) +
                                               #ylim(c(0, 2.5)) +
                                               theme_bw() +
                                               theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                                               axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2))) +
                                               scale_linetype_manual(values = c(1, 2, 3))
                                               #theme(legend.position="none")
turns_all_side_front_plot
################################################################################






# #  Point plots
# steps_pt_both <- ggplot(both, aes(steps, iter_num, colour = type)) +
                             # geom_point(alpha = 0.3) +
                             # xlab("Step length (km)") +
                             # ylab("Iteration number") +
                             # theme_bw() +
                             # theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                             # axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
# steps_pt_both_col <- steps_pt_both +
                                   # scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                                   # scale_fill_manual(values = c("grey47","grey87"), guide = FALSE)
# turns_pt_both <- ggplot(both, aes(turns, iter_num, colour = type)) +
                             # geom_point(alpha = 0.3) +
                             # xlab("Turn angle (rad)") +
                             # ylab("Iteration number") +
                             # theme_bw() +
                             # theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                             # axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(2)))
# turns_pt_both_col <- turns_pt_both +
                                   # scale_colour_manual(values = c("grey47","grey87"), guide = FALSE) +
                                   # scale_fill_manual(values = c("grey47","grey87"), guide = FALSE)
################################################################################

