SliceTAC <- function(t, times, greed) {
  answer <- list()
  for(n in 1:(length(times)-1)) {
    extra <- 0.2 * times[n]
    answer[[n]] <- SliceSegment(t, (1-greed)*times[n], (1+greed)*times[n+1])
  }
  
  return(answer)
}

SliceSegment <- function(t, t1, t2) {
  return(t>=t1&t<=t2)
}

StitchTAC <- function(tac, cuts) {
  for(n in 1:(length(tac))) {
    boolvec <- tac[[n]]$x > cuts[n] & tac[[n]]$x <= cuts[n+1]
    tac[[n]]$x <- tac[[n]]$x[boolvec]
    tac[[n]]$y <- tac[[n]]$y[boolvec]
  }
  
  output <- cbind(tac[[1]]$x, tac[[1]]$y)
  for(n in 2:(length(tac))) {
    output <- rbind(output, cbind(tac[[n]]$x, tac[[n]]$y))
  }
  return(output)
}