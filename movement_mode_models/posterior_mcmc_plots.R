# Sara Williams
# 7/26/2016; 1/20/2017
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

#  MCMC plots with package 'ggmcmc'
s_1 <- ggs(single_fit, family = "^mu|^rho")
s_2 <- ggs(single_fit, family = "^lambda|^v")

dens_s2 <- ggs_density(s_2) + 
                   theme_bw() +
                   theme(panel.grid.minor = element_blank(), 
                               strip.text = element_text(size = rel(1.5)), 
                               axis.text = element_text(size = rel(1.5)), 
                               axis.title = element_text(size = rel(1.5)),
                               legend.position="none")

dens_s1 <- ggs_density(s_1) + 
                   theme_bw() +
                   theme(panel.grid.minor = element_blank(), 
                               strip.text = element_text(size = rel(1.5)), 
                               axis.text = element_text(size = rel(1.5)), 
                               axis.title = element_text(size = rel(1.5)),
                               legend.position="none")


grid.arrange(dens_s1, dens_s2, ncol=2, nrow=2)

trace_s1 <- ggs_traceplot(s_1) + 
                    theme_bw() +
                    theme(panel.grid.minor = element_blank(), 
                                strip.text = element_text(size = rel(1.5)), 
                                axis.text = element_text(size = rel(1.5)), 
                                axis.title = element_text(size = rel(1.5)),
                                legend.position="none")

trace_s2 <- ggs_traceplot(s_2) + 
                    theme_bw() +
                    theme(panel.grid.minor = element_blank(), 
                                strip.text = element_text(size = rel(1.5)), 
                                axis.text = element_text(size = rel(1.5)), 
                                axis.title = element_text(size = rel(1.5)),
                                legend.position="none")

grid.arrange(trace_s1, dens_s1, trace_s2, dens_s2, ncol=2, nrow=2)



su_1 <- ggs(single_fit_surf, family = "^mu|^rho")
su_2 <- ggs(single_fit_surf, family = "^lambda|^v")
d_1 <- ggs(single_fit_dive, family = "^mu|^rho")
d_2 <- ggs(single_fit_dive, family = "^lambda|^v")


dens_su2 <- ggs_density(su_2) + 
                     theme_bw() +
                     theme(panel.grid.minor = element_blank(), 
                                 strip.text = element_text(size = rel(1.5)), 
                                 axis.text = element_text(size = rel(1.5)), 
                                 axis.title = element_text(size = rel(1.5)),
                                 legend.position="none")

dens_su1 <- ggs_density(su_1) + 
                     theme_bw() +
                     theme(panel.grid.minor = element_blank(), 
                                 strip.text = element_text(size = rel(1.5)), 
                                 axis.text = element_text(size = rel(1.5)), 
                                 axis.title = element_text(size = rel(1.5)),
                                 legend.position="none")

grid.arrange(dens_su1, dens_su2, ncol=1, nrow=2)

dens_d2 <- ggs_density(d_2) + 
                   theme_bw() +
                   theme(panel.grid.minor = element_blank(), 
                               strip.text = element_text(size = rel(1.5)), 
                               axis.text = element_text(size = rel(1.5)), 
                               axis.title = element_text(size = rel(1.5)),
                               legend.position="none")

dens_d1 <- ggs_density(d_1) + 
                   theme_bw() +
                   theme(panel.grid.minor = element_blank(), 
                               strip.text = element_text(size = rel(1.5)), 
                               axis.text = element_text(size = rel(1.5)), 
                               axis.title = element_text(size = rel(1.5)),
                               legend.position="none")

grid.arrange(dens_d1, dens_d2, ncol=1, nrow=2)


trace_su1 <- ggs_traceplot(su_1) + 
                      theme_bw() +
                      theme(panel.grid.minor = element_blank(), 
                                  strip.text = element_text(size = rel(1.5)), 
                                  axis.text = element_text(size = rel(1.5)), 
                                  axis.title = element_text(size = rel(1.5)),
                                  legend.position="none")

trace_su2 <- ggs_traceplot(su_2) + 
                      theme_bw() +
                      theme(panel.grid.minor = element_blank(), 
                                  strip.text = element_text(size = rel(1.5)), 
                                  axis.text = element_text(size = rel(1.5)), 
                                  axis.title = element_text(size = rel(1.5)),
                                  legend.position="none")

grid.arrange(trace_s1, dens_s1, trace_s2, dens_s2, ncol=2, nrow=2)














d_1 <- ggs(double_fit, family = "^mu|^rho")
d_2 <- ggs(double_fit, family = "^lambda|^v")
d_3 <- ggs(double_fit, family = "^beta|^eps")

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




