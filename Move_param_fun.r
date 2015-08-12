# Code from Josh Nowak - adjusted to remove need for dates and time interval
# 8/6/2015
# Function to generate step length and turning angles
################################################################################

#  This function takes the telemetry data and calculates bearing, step length, 
#  time between fixes in hours and turn angles.  Each calculation should be 
#  interpreted as the step length it took to get to where you are now.  In other
#  words, the data on line two is the result of what happened between row 1 and 
#  row 2.  The function does not accept lists for the loc.data.




move.fun <- function(loc.data){

        # Create three datasets
        len <- nrow(loc.data)
        d1 <- loc.data[-c(len-1, len),c("x","y")]
        d2 <- loc.data[-c(1, len),c("x","y")]
        d3 <- loc.data[-c(1:2),c("x","y")]
        
        #  Function for calculating steps and turns       
       anglefun <- function(xx,yy,bearing=TRUE,as.deg=FALSE){

        ## calculates the compass bearing of the line between two points
        ## xx and yy are the differences in x and y coordinates between two points
        ## Options:
        ## bearing = FALSE returns +/- pi instead of 0:2*pi
        ## as.deg = TRUE returns degrees instead of radians
        c = 1
        if (as.deg){
          c = 180/pi
        }
  
        b<-sign(xx)
        b[b==0]<-1  #corrects for the fact that sign(0) == 0
        tempangle = b*(yy<0)*pi+atan(xx/yy)
        if(bearing){
        #return a compass bearing 0 to 2pi
        #if bearing==FALSE then a heading (+/- pi) is returned
        tempangle[tempangle<0]<-tempangle[tempangle<0]+2*pi
        }
        return(tempangle*c)
        }

        bearing.ta <- function(loc1,loc2,loc3,as.deg=FALSE){
        ## calculates the bearing and length of the two lines
        ##    formed by three points
        ## the turning angle from the first bearing to the
        ##    second bearing is also calculated
        ## locations are assumed to be in (X,Y) format.
        ## Options:
        ## as.deg = TRUE returns degrees instead of radians
        if (length(loc1) != 2 | length(loc2) != 2 | length(loc3) !=2){
        print("Locations must consist of either three vectors, length == 2, or three two-column dataframes")
        return(NaN)
        }
        c = 1
        if (as.deg){
        c = 180/pi
        }
  
        locdiff1 <- loc2-loc1
        if(any(locdiff1$x == 0 | locdiff1$y == 0)){
            locdiff1[locdiff1$x == 0 | locdiff1$y == 0,] <- 0.01
        }
        locdiff2 <- loc3-loc2
        if(any(locdiff2$x == 0 | locdiff2$y == 0)){
            locdiff2[locdiff2$x == 0 | locdiff2$y == 0,] <- 0.01
        }
        bearing1<-anglefun(locdiff1[1],locdiff1[2],bearing=F)
        bearing2<-anglefun(locdiff2[1],locdiff2[2],bearing=F)

        if(is.data.frame(locdiff1)){
        dist1<-sqrt(rowSums(locdiff1^2))
        dist2<-sqrt(rowSums(locdiff2^2))
        }else{
        dist1<-sqrt(sum(locdiff1^2))
        dist2<-sqrt(sum(locdiff2^2))
        }
  
        ta= (bearing2-bearing1)
        
        ta[ta < -pi,1] <- ta[ta < -pi,1] + 2*pi
        ta[ta > pi,1] <- ta[ta > pi,1] - 2*pi
        
        return(list(bearing1=unlist(bearing1*c),bearing2=unlist(bearing2*c),
        ta=unlist(ta*c),dist1=unlist(dist1),dist2=unlist(dist2)))
        }
        
        xxx <- bearing.ta(d1, d2, d3)
        loc.data$turn.angle <- c(NA, NA, xxx$ta)
        loc.data$heading <- c(NA, xxx$bearing1, 
                              xxx$bearing2[length(xxx$bearing2)])
        loc.data$dist <- c(NA, xxx$dist1, xxx$dist2[length(xxx$dist2)])
              
        return(loc.data)
        }
        