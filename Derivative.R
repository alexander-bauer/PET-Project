#returns: (time, ours, modelmemb, ideal)
SmoothCurve <- function(t, tissue, plasma, V0) {
  tac <- tissue - V0 * plasma
  
  tac.max <- max(tac)
  i.max <- which.max(tac)
  t.max <- t[i.max]
  
  temp <- tac > 0.50 * tac.max
  halfmax <- t[temp][1]
  
  temp <- tac > 0.90 * tac.max
  ninetymax <- t[temp][1]
  endninetymax <- tail(t[temp], 1)  
  
  temp <- tac == tac.max
  max <- t[temp]
  
  cuts <- c(halfmax, ninetymax, endninetymax, endninetymax+20, tail(t, 1))
  
  head <- SmoothHead(t, tissue, plasma, V0, cuts[1], cuts[2])
  peak <- SmoothPeak(t, tissue, plasma, V0, cuts[2], cuts[3])
  mid  <- SmoothMid (t, tissue, plasma, V0, cuts[3], cuts[4])
  end  <- SmoothEnd (t, tissue, plasma, V0, cuts[4], cuts[5])
  
  ans <- head[head[,1] > cuts[1] & head[,1] <= cuts[2],]
  ans <- rbind(ans, peak[peak[,1] > cuts[2] & peak[,1] <= cuts[3],])
  ans <- rbind(ans, mid[mid[,1] > cuts[3] & mid[,1] <= cuts[4],])
  ans <- rbind(ans, end[end[,1] > cuts[4] & mid[,1] <= cuts[5],])
  
  return(ans)
}

SmoothHead <- function(t, tissue, plasma, V0, t1, t2) {
  return(GetDerivative(t, tissue, plasma, V0, t1, t2, 0.4, 2, 1, 51))
}

SmoothPeak <- function(t, tissue, plasma, V0, t1, t2) {
  return(GetDerivative(t, tissue, plasma, V0, t1, t2, 0.4, 2, 5, 51))
}

SmoothMid <- function(t, tissue, plasma, V0, t1, t2) {
  return(GetDerivative(t, tissue, plasma, V0, t1, t2, 0.4, 2, 2, 51))
}

SmoothEnd <- function(t, tissue, plasma, V0, t1, t2) {
  return(GetDerivative(t, tissue, plasma, V0, t1, t2, 0.4, 1, 7, 51))
}

GetDerivative <- function(t, tissue, plasma, V0, t1, t2,
                          greed, degree, bandwidth, gridsize=51L) {
  idx <- SliceTAC(t,c(t1,t2),greed)[[1]]
  brain.rac <- pmax(0, tissue - V0 * plasma)
  d.dt <- locpoly(t[idx],brain.rac[idx],drv=1,degree=degree,
                  bandwidth=bandwidth, gridsize=gridsize)
  
  return(cbind(d.dt$x, d.dt$y))
}