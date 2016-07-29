# Sara Williams
# 7/26/2016; 
# Generate plots from MCMC posterior distribution object (iterations)
#   and run in JAGS.
################################################################################

#  Load packages
library(ggplot2)
library(plyr)
library(dplyr)
library(packcircles)
library(ggmcmc)
library(rjags)
library(gridExtra)
library(ggthemes)
################################################################################

#  Load simulated step length and turning angle values generated from MCMC posterior distribtions
load("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/movement_modes/model_out/post_sims.RData")
################################################################################

#  Generate vectors to hold data with package 'packcircles'
r <- post_sims$steps
id <- seq(1, length(r), 1)
x <- rep(0, length(r)) 
y <- rep(0, length(r)) 
points <- data.frame(x, y, r, id)
points <- dplyr::arrange(points, by = r)

#  Generate circle vertices for each iterations step length
#   Here, each iteration's step length is the radius of the circle
dat <- circlePlotData(points, npoints = 100, xyr.cols = 1:3, id.col = 4)

#  Plot all of thesee circles (45000 circles from the 45000 simulated values from 45000 MCMC iterations)
p <- ggplot() +
        geom_polygon(data=dat, aes(x, y, group=id), fill = "grey47", alpha = 0.1, linetype = 0) +
        coord_equal() +
        theme_bw() +
        xlab("X (km)") +
        ylab("Y (km)") +
        theme(panel.grid.minor = element_blank(), axis.text = element_text(size = rel(2)), 
        axis.title = element_text(size = rel(2)))
p

################################################################################
library(gridExtra)

#  MCMC plots with package 'ggmcmc'
s_1 <- ggs(single_fit, family = "^mu|^rho")
s_2 <- ggs(single_fit, family = "^lambda|^v")
d_1 <- ggs(double_fit, family = "^mu|^rho")
d_2 <- ggs(double_fit, family = "^lambda|^v")
d_3 <- ggs(double_fit, family = "^beta|^eps")

trace_s1 <- ggs_traceplot(s_1) + 
                    theme_bw() +
                    theme(panel.grid.minor = element_blank(), 
                                strip.text = element_text(size = rel(2)), 
                                axis.text = element_text(size = rel(2)), 
                                axis.title = element_text(size = rel(2)),
                                legend.position="none")

dens_s1 <- ggs_density(s_1) + 
                   theme_bw() +
                   theme(panel.grid.minor = element_blank(), 
                               strip.text = element_text(size = rel(2)), 
                               axis.text = element_text(size = rel(2)), 
                               axis.title = element_text(size = rel(2)),
                               legend.position="none")

trace_s2 <- ggs_traceplot(s_2) + 
                    theme_bw() +
                    theme(panel.grid.minor = element_blank(), 
                                strip.text = element_text(size = rel(2)), 
                                axis.text = element_text(size = rel(2)), 
                                axis.title = element_text(size = rel(2)),
                                legend.position="none")

dens_s2 <- ggs_density(s_2) + 
                   theme_bw() +
                   theme(panel.grid.minor = element_blank(), 
                               strip.text = element_text(size = rel(2)), 
                               axis.text = element_text(size = rel(2)), 
                               axis.title = element_text(size = rel(2)),
                               legend.position="none")

grid.arrange(trace_s1, dens_s1, trace_s2, dens_s2, ncol=2, nrow=2)

trace_d1 <- ggs_traceplot(d_1) + 
                    theme_bw() +
                    theme(panel.grid.minor = element_blank(), 
                                strip.text = element_text(size = rel(2)), 
                                axis.text = element_text(size = rel(2)), 
                                axis.title = element_text(size = rel(2)),
                                legend.position="none")

dens_d1 <- ggs_density(d_1) + 
                   theme_bw() +
                   theme(panel.grid.minor = element_blank(), 
                               strip.text = element_text(size = rel(2)), 
                               axis.text = element_text(size = rel(2)), 
                               axis.title = element_text(size = rel(2)),
                               legend.position="none")

trace_d2 <- ggs_traceplot(d_2) + 
                    theme_bw() +
                    theme(panel.grid.minor = element_blank(), 
                                strip.text = element_text(size = rel(2)), 
                                axis.text = element_text(size = rel(2)), 
                                axis.title = element_text(size = rel(2)),
                                legend.position="none")

dens_d2 <- ggs_density(d_2) + 
                   theme_bw() +
                   theme(panel.grid.minor = element_blank(), 
                               strip.text = element_text(size = rel(2)), 
                               axis.text = element_text(size = rel(2)), 
                               axis.title = element_text(size = rel(2)),
                               legend.position="none")

trace_d3 <- ggs_traceplot(d_3) + 
                    theme_bw() +
                    theme(panel.grid.minor = element_blank(), 
                                strip.text = element_text(size = rel(2)), 
                                axis.text = element_text(size = rel(2)), 
                                axis.title = element_text(size = rel(2)),
                                legend.position="none")

dens_d3 <- ggs_density(d_3) + 
                   theme_bw() +
                   theme(panel.grid.minor = element_blank(), 
                               strip.text = element_text(size = rel(2)), 
                               axis.text = element_text(size = rel(2)), 
                               axis.title = element_text(size = rel(2)),
                               legend.position="none")

grid.arrange(trace_d1, dens_d1, ncol=2, nrow=1)
grid.arrange(trace_d2, dens_d2, ncol=2, nrow=1)
grid.arrange(trace_d3, dens_d3, ncol=2, nrow=2)




cor_plot <- ggs_autocorrelation(s)+ 
                   theme_fivethirtyeight() +
                   theme(panel.grid.minor = element_blank(), axis.text = element_text(size = rel(2)), 
                   axis.title = element_text(size = rel(2)))
cor_plot

