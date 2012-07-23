SmoothTracerInTissue <- function(t, tac) {
  
  tac.max <- max(tac)
  i.max <- which.max(tac)
  t.max <- t[i.max]
  
  temp <- tac > 0.05 * tac.max
  start <- t[temp][1]
  
  temp <- tac > 0.25 * tac.max
  quartermax <- t[temp][1]
  
  temp <- tac > 0.5 * tac.max
  halfmax <- t[temp][1]
  
  temp <- tac > 0.75 * tac.max
  threequartermax <- t[temp][1]
  
  temp <- tac == tac.max
  max <- t[temp]
  
  cuts <- c(0, halfmax, max, tail(t, 1))
  greed <- 0.2
  
  boolvec <-SliceTAC(t, cuts[1:2], greed)[[1]]
  timeslice <- t[boolvec]
  slice <- tac[boolvec]
  
  arg.t <- timeslice
  arg.times <- c(arg.t[1], exp(c(-2, -1.5)), tail(arg.t, 1))
  arg.schedule <- c(0.05, 0.1, 0.15, 0.2)
  
  arg.bandwidth <- GetBandwidths(arg.t, arg.times,
                                 arg.schedule , 51)
  
  smoothed0 <- locpoly(timeslice, slice, degree=2, bandwidth=0.5, gridsize=51)
  drvsmoothed0 <- locpoly(timeslice, slice, drv=1, degree=2,
                          bandwidth=arg.bandwidth,
                          gridsize=51)
  
  
  
  boolvec <-SliceTAC(t, cuts[2:3], greed)[[1]]
  timeslice <- t[boolvec]
  slice <- tac[boolvec]
  
  arg.t <- timeslice
  arg.times <- c(arg.t[1], tail(arg.t, 1))
  arg.schedule <- c(0.6, 2)
  
  arg.bandwidth <- GetBandwidths(arg.t, arg.times,
                                 arg.schedule , 51)
  
  smoothed1 <- locpoly(timeslice, slice, degree=2, bandwidth=1, gridsize=51)
  drvsmoothed1 <- locpoly(timeslice, slice, drv=1, degree=2,
                          bandwidth=arg.bandwidth, gridsize=51)
  
  
  boolvec <-SliceTAC(t, cuts[3:4], greed)[[1]]
  timeslice <- t[boolvec]
  slice <- tac[boolvec]
  
  smoothed2 <- locpoly(timeslice, slice, degree=2, bandwidth=4, gridsize=51)
  drvsmoothed2 <- locpoly(timeslice, slice, drv=1, degree=2, bandwidth=4, gridsize=51)
  
  
  out <- list(smoothed0, smoothed1, smoothed2)
  output <- StitchTAC(out, cuts)
  
  outdrv <- list(drvsmoothed0, drvsmoothed1, drvsmoothed2)
  outputdrv <- StitchTAC(outdrv, cuts)
  
  data <- approx(t, tac, xout=output[,1])$y
  
  outputlist <- list(cbind(output, data), cbind(outputdrv[,1], outputdrv[,2]))
  
  return(outputlist)
}

SmoothTracerInTissuePeak <- function(t, tac) {
  
  tac.max <- max(tac)
  i.max <- which.max(tac)
  t.max <- t[i.max]
  
  temp <- tac > 0.50 * tac.max
  halfmax <- t[temp][1]
  
  temp <- tac > 0.75 * tac.max
  threequartermax <- t[temp][1]
  
  temp <- t > 1.50 * t.max
  endnearpeak <- t[temp][1]
  
  temp <- tac == tac.max
  max <- t[temp]
  
  cuts <- c(halfmax, endnearpeak, tail(t, 1))
  greed <- 0.2
  
  #ActOnSlice (time, tac, t1, t2, greed, c(curvedegree, drvdegree),
  #           c(curvebandwidth, drvbandwidth), gridsize)
  
  smooth0 <- ActOnSlice(t, tac, cuts[1], cuts[2], 0.2, c(2, 2), c(4, 5), 51)
  
  #smooth1 <- ActOnSlice(t, tac, cuts[2], cuts[3], 0.2, c(2, 2), c(4, 5), 51)
  
  out <- list(smooth0[,1], smooth1[,1])
  output <- StitchTAC(out, cuts)
  
  outdrv <- list(smooth0[,2], smooth1[,2])
  outputdrv <- StitchTAC(outdrv, cuts)
  
  data <- approx(t, tac, xout=output[,1])$y
  
  outputlist <- list(cbind(output, data), cbind(outputdrv[,1], outputdrv[,2]))
  
  return(outputlist)
}

ActOnSlice <- function(time, tac, t1, t2, greed, degree, bandwidth, grid) {
  boolvec <-SliceTAC(time, c(t1, t2), greed)[[1]]
  timeslice <- time[boolvec]
  slice <- tac[boolvec]
  
  smoothed <- locpoly(timeslice, slice, degree=degree[1], bandwidth=bandwidth[1],
                      gridsize=grid)
  drvsmoothed <- locpoly(timeslice, slice, drv=1, degree=degree[2],
                         bandwidth=bandwidth[2], gridsize=grid)
  
  return(cbind(smoothed, drvsmoothed))
}

GetBandwidths <- function(t, times, schedule, gridsize) {
  
  xout <- (t[1] + (0:(gridsize-1)) * (tail(t, 1) - t[1])/(gridsize-1))
  
  xout[1] <- t[1]
  xout[gridsize] <- tail(t,1)
  
  return(approx(times, schedule, 
                xout=xout)$y)
}