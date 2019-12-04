

#  Load packages.
library(Distance)
library(dplyr)

#  Read in final data set. NOTE that 85th percentil trucation is done WITHIN the ds() function, 
#   so we do NOT use final_dat_truncated85
final_dat <- read.csv("C:/Users/saraw/Desktop/Whale/final_dat_2015.csv")


#  Truncate to the 85th percentile so functions only go to the range of those estimated
#   through detection function models.
percentile85 <- quantile(final_dat$distance, probs=0.85, na.rm=TRUE)
percentile85
final_dat_truncated85 <- final_dat[final_dat$distance < 4565,]


#  Generate values from each model to plot across distance range (0 to 4565)

#  Distance values (m)
x <- seq(min(final_dat_truncated85$distance), 
         max(final_dat_truncated85$distance))

#  Detection probability values creation.

#  Visibility
se_vis_no <- 0.04457347
ci_vis_no <- 1.96*se_vis_no

hr.vis <- function(x, scale = exp(6.860), shape = exp(0.780)){
               1 - exp(-(x/scale)^(-shape))
      } # excellent - baseline
      hr.vis_l <- function(x, scale = exp(6.860-ci_vis_no), shape = exp(0.780)){
               1 - exp(-(x/scale)^(-shape))
      }
      hr.vis_u <- function(x, scale = exp(6.860+ci_vis_no), shape = exp(0.780)){
               1 - exp(-(x/scale)^(-shape))
      }

vis_dat <- as.data.frame(cbind(x, hr.vis(x), hr.vis_l(x), hr.vis_u(x)))


vis_plot <- ggplot(vis_dat, aes(x=x, y=V2)) +
                   geom_line(size = 1) +
                   geom_ribbon(data=vis_dat, aes(x=x, ymin=V3,ymax=V4),
                   alpha=0.4,  fill="grey70") +
                   xlab("\nDistance from ship (m)") +
                   ylab("Detection probability under excellent visibility \n") +
                   theme_bw() +
                   theme(panel.grid.major.x=element_blank(),
                   panel.grid.minor=element_blank()) +
                   theme(axis.title=element_text(size=15), 
                  axis.text=element_text(size=12))  +
                  scale_x_continuous(expand = c(0, 0)) +
                  scale_y_continuous(expand = c(0, 0.01))
vis_plot
