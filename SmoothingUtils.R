library(KernSmooth)

# Smooth a plasma time/activity curve by applying locpoly to overlapping
# intervals and stitching the results together. Degree and bandwidth parameters
# passed to locpoly depend on the interval.
SmoothPlasmaTac <- function(time, activity){
  ca <- activity  
  tau <- time # renaming because I'm pasting code from scripts
  i.max <- which.max(ca)
  lead <- tau < tau[i.max] & ca < 0.05*max(ca)
  i.onset <- which.max(tau[lead])
  # Smooth from start to just beyond bolus onset:
  idx <- 1:(i.onset+3)
  temp <- locpoly(tau[idx],ca[idx],degree=1,bandwidth=.05,gridsize=50)
  idx <- temp$x <= tau[i.onset+1]
  temp$x <- temp$x[idx]
  temp$y <- temp$y[idx]
  t.last <- tail(temp$x,1)
  ans <- cbind(temp$x,temp$y)
  # The peak. Smooth from just beyond bolus onset to well beyond the peak,
  # but discard left and right extremes.
  idx <- (i.onset):(i.max+10)
  temp <- locpoly(tau[idx],ca[idx],degree=2,bandwidth=.05,gridsize=50)
  idx <- temp$x > t.last & temp$x < tau[i.max+5]
  temp$x <- temp$x[idx]
  temp$y <- temp$y[idx]
  t.last <- tail(temp$x,1)
  ans <- rbind(ans,cbind(temp$x,temp$y))
  # From the just beyond peak to somewhat later
  t.max <- tau[i.max]
  idx <- tau > t.max & tau <= 12  
  temp <- locpoly(tau[idx],ca[idx],degree=3,bandwidth=.5,gridsize=50)
  idx <- temp$x > t.last & temp$x <= 6 # 9
  temp$x <- temp$x[idx]
  temp$y <- temp$y[idx]
  t.last <- tail(temp$x,1)
  ans <- rbind(ans,cbind(temp$x,temp$y))
  # From 2 to 20 minutes, keeping from 3 to 15.
#   idx <- tau > 2 & tau <=20
  idx <- tau > t.last - 2 & tau <=20
  temp <- locpoly(tau[idx],ca[idx],degree=1,bandwidth=2,gridsize=50)
  idx <- temp$x > t.last & temp$x <= 15
  temp$x <- temp$x[idx]
  temp$y <- temp$y[idx]
  t.last <- tail(temp$x,1)
  ans <- rbind(ans,cbind(temp$x,temp$y))
  # From 10 to 45, keeping from 15 to 20
  idx <- tau > 10 & tau <= 45
  temp <- locpoly(tau[idx],ca[idx],degree=1,bandwidth=5,gridsize=50)
  idx <- temp$x > t.last & temp$x <= 20
  temp$x <- temp$x[idx]
  temp$y <- temp$y[idx]
  t.last <- tail(temp$x,1)
  ans <- rbind(ans,cbind(temp$x,temp$y))
  # From 2 to the end, using a variable bandwidth, using data from 2
  # minutes on for context, but discarding 2 thru 20.
  dt = (tail(tau,1)-2)/50
  times <- c(2,20,40, 50,tail(tau,1))
  bandwidths <- c(2,5,12,16,30)
  bws <- interpolate_bw(times,bandwidths,min(times),max(times),dt)
  idx <- tau >= 2
  temp <- locpoly(tau[idx],ca[idx],degree=1,bandwidth=bws,gridsize=length(bws))
  idx <- temp$x > t.last
  temp$x <- temp$x[idx]
  temp$y <- temp$y[idx]
  ans <- rbind(ans,cbind(temp$x,temp$y))
  return(ans)
}

# Interpolate a bandwidth schedule between t1 and t2 at sampling
# interval dt. Parameters times and bandwidths are vectors of
# the same length.
interpolate_bw <- function(times,bandwidths,t1,t2,dt){
  gridsize <- dt2gridsize(t1,t2,dt)
  xout <- seq(from=t1,to=t2,length.out=gridsize)
  temp <- approx(times,bandwidths,xout)
  return(temp$y)
}


# Calculate gridsize given t1, t2, and a sampling interval dt.
dt2gridsize <- function(t1,t2,dt){
  gridsize <- round((t2-t1)/dt) 
  gridsize <- gridsize + 1 - gridsize%%2 # always odd
  return(gridsize)
}
