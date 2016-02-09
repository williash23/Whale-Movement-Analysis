library(crawl)
library(dplyr)

temp1 <- read.csv("C:/Users/sara.williams/Documents/GitHub/Whale-Movement-Analysis/data/Whales_0615_locations_clean.csv")
temp2 <- temp1 %>%
				   group_by(same_whale_ID) %>%
				   filter(n() > 1) %>%
				   ungroup(.) %>%
				   arrange(same_whale_ID, ob_order_time) %>%
				   as.data.frame(.)
temp3 <- temp2 %>%
				   dplyr::select(X_whale_UTM, Y_whale_UTM, same_whale_ID) 
Time <- seq(0, (nrow(temp3)-1), 1)
dat <- cbind(temp3, Time)				   


## Initial state values
init <- list(
  a1.x=c(dat$X_whale_UTM[1],0),
  a1.y=c(dat$Y_whale_UTM[1],0),
  P1.x=diag(c(1,1)),
  P1.y=diag(c(1,1))
)

tst <- dat[3:6,]

fit1 <- crwMLE(mov.model=~1, err.model=NULL, drift.model=FALSE,
           data=tst, coord=c("X_whale_UTM", "Y_whale_UTM"), polar.coord=FALSE,
           Time.name="Time", initial.state=init, 
           control=list(maxit=2000,trace=1, REPORT=10))

pred_time <- seq(6,9,1)
#seq((nrow(temp3)), (nrow(temp3)+5), 1)		   
pred <- crwPredict(fit1, pred_time)
crwPredictPlot(pred)		   

set.seed(123)
simObj <- crwSimulator(fit1, pred_time, parIS=100, df=20, scale=18/20)
w <- simObj$thetaSampList[[1]][,1]
dev.new()
hist(w*100, main='Importance Sampling Weights', sub='More weights near 1 is desirable')	

##Approximate number of independent samples
round(100/(1+(sd(w)/mean(w))^2))

dev.new(bg=gray(0.75))
jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
crwPredictPlot(pred, 'map')

## Sample 20 tracks from posterior predictive distribution
iter <- 20
cols <- jet.colors(iter)
for(i in 1:iter){
  samp <- crwPostIS(simObj)
  lines(samp$alpha.sim.x[,'mu'], samp$alpha.sim.y[,'mu'],col=cols[i])
}
	   
		   
out <- dat %>%
			 group_by(same_whale_ID) %>%
		     do(fit <- crwMLE(mov.model=~1, err.model=NULL, drift.model=FALSE,
             data=dat, coord=c("X_whale_UTM", "Y_whale_UTM"), polar.coord=FALSE,
             Time.name="Time", initial.state=init, 
             control=list(maxit=2000,trace=1, REPORT=10)))
