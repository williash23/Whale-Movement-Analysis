#  Create Random points by building and sampling from an empirical distribution

# random empirical sample
r_emp = function(data.out, mids.lower.upper = "upper", b = 500, nout){
    n = length(data.out)
    a.r = runif(nout,0,n)
    a.h = hist(data.out, plot = FALSE, breaks = b)
    cumsum.bin = cumsum(a.h$counts)
    selected = sapply(a.r, function(x) sum(x >= cumsum.bin)+1)
    if (mids.lower.upper == "mids") {
         out = a.h$mids[selected] # These are the random distances drawn from the empirical distribution
       } else if (mids.lower.upper == "lower" ) {
         out = a.h$breaks[selected]
       } else if (mids.lower.upper == "upper" ) {
         out = a.h$breaks[selected+1]
       }
    out
}

###  Turn angles
samp.turn <- function(data.out, mids.lower.upper = "mids", b = 100, nout){
    n = length(data.out)
    a.r = runif(nout,0,n)
    a.h = hist(data.out, plot = FALSE, breaks = b)
    cumsum.bin = cumsum(a.h$counts)
    selected = sapply(a.r, function(x) sum(x >= cumsum.bin)+1)
    if (mids.lower.upper == "mids") {
         out = a.h$mids[selected] # These are the random distances drawn from the empirical distribution
       } else if (mids.lower.upper == "lower" ) {
         out = a.h$breaks[selected]
       } else if (mids.lower.upper == "upper" ) {
         out = a.h$breaks[selected+1]
       }
    out
}

destination <-function(x,y,dist,bearing,as.deg=FALSE){
  if(as.deg){
    ##if bearing is in degrees, convert to radians
    bearing=bearing*pi/180
  }
  newx<-x+dist*sin(bearing)  ##X
  newy<-y+dist*cos(bearing)  ##Y
  return(data.frame(x = newx, y = newy))
}
################################################################################
#  To run the above

    create.rand <- function(data.out, as.deg=FALSE, turn = F){
          tmp1 <- r_emp(data.out$dist)
          if(turn == T){
            tmp2 <- samp.turn(data.out$turn.angle, mids.lower.upper = "mids",
                            length.out = 25)
            tmp3 <- destination(data.out$x, data.out$y, tmp1, tmp2, as.deg = F)
          }else(tmp3 <- destination(data.out$x, data.out$y, tmp1, 
                          runif(nrow(data.out), -3.14, 3.14), as.deg = F) )
          
          c.rand <- data.frame(x = tmp3$x, y = tmp3$y)
      return(c.rand)
      }
      
