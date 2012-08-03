SolveBBB <- function(tm, ca.star, m.star, V0, K1, k2){
  # Estimates me, and mb for the specified tissue tac, m.star, plasma tac,
  # ca.star, sample times, tm, and constants V0, K1, k2, using the BBB
  # and mass balance equations only. Specifically, it uses
  #              mb = m*-V0xca* + d/dt(m*-V0xca*)/k2 - ca*xK1/k2,
  # where the derivative is estimated numerically by localized polynomial
  # approximation. Also returns the numerical derivative of tracer in tissue 
  # and intermediate results of interest.
  #
  # Args:
  #   tm -- times at which ca and m are sampled
  #   ca.star -- plasma activity curve
  #   m.star -- activity curve for the region of interest
  #   V0 -- apparent vascular volume of distribution for the region of interest
  #   K1, k2 -- clearance and rate constant for the blood/brain barrier
  #
  # Returns a list containing:
  #   mb -- activity of bound tracer as estimated by 
  #         numerical solution of, a vector
  #   me -- activity of exchangeable tracer, tracer in tissue minus me, a vector
  #   d.dt -- estimated derivative of tracer in tissue, a vector
  #   in.tissue -- cumulative integral of d.dt, with appropriate constant 
  #         of integration, a smooth approximation to m*-V0ca*
  #   raw.in.tissue -- m*-V0ca* interpolated to time base t
  #   ca -- ca.star interpolated to time base t
  #   t  -- time base of me, mb, d.dt, in.tissue, raw.in.tissue, ca

  # Estimate the derivative of m*-V0ca*. Result is a matrix whose
  # first column is time and whose second column is activity/min.
  d.dt <- SmoothCurve(tm, m.star, ca.star, V0)
  # Put tracer in tissue on the same time base as its derivative
  raw.in.tissue <- approx(tm,
                          m.star-V0*ca.star,
                          xout = d.dt[,1])$y
  # Reintegrate d.dt and adjust constant of integration
  in.tissue <- qdint(d.dt[,1],d.dt[,2])
  in.tissue <- in.tissue + mean(raw.in.tissue-in.tissue)
  # Put plasma on the same time base
  ca <- approx(tm, ca.star, xout = d.dt[,1])$y  
  # Calculate mb
  mb <- in.tissue + d.dt[,2]/k2 - ca*K1/k2
  # Calculate me
  me <- in.tissue - mb
  return(list(mb = mb,
              me = me,
              t  = d.dt[,1],
              d.dt = d.dt[,2],
              in.tissue = in.tissue,
              raw.in.tissue = raw.in.tissue,
              ca = ca))
}

LoessBBB <- function(bbb.solution, span=.5){
  # Applies loess smoothing to bbb.solution$in.tissue, bbb.solution$mb,
  # and bbb.solution$me. Returns the named list, bbb.solution, with
  # smoothed curves added as named elements in.tissue.loess, mb.loess, 
  # and me.loess, respectively.
  #
  # Args:
  #  bbb.solution -- output of BBBSolve or similar function
  #  span -- span argument for loess
  #
  # Returns the input list augmented by:
  #  in.tissue.loess
  #  me.loess
  #  mb.loess
  tm <- bbb.solution$t
  for(n in list("in.tissue", "mb", "me")){
    nout <- paste(n,"loess",sep=".")
    temp <- bbb.solution[[n]]
    bbb.solution[[nout]] <- loess(temp ~ tm, span=span)$fitted
  }
  return(bbb.solution)
}
