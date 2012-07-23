# Provisional methods for estimating Ve from cerebellum and, given Ve,
# for estimating V0, K1, and k2 (= K1/Ve) in regions of interest. The
# main functions are V0K1k2inCerebellum(t, ca.star, m.star), which returns
# V0, K1, and k2 for cerebellum from which Ve = K1/k2 may be calculated, and
# K1andV0(t, ca.star, m.star, Ve), which returns estimates of K1 and V0 for
# a region of interest.
# 
# Both methods operate by finding coefficients of the integral form of
# a standard model consisting of one vascular and one tissue compartment.
# In the case of regions of interest, in which binding can occur, the
# model is applied only to data at early times when binding is negligible
# and the model approximately holds. In the case of cerebellum, the model 
# is exact and is applied to all data.
# 
# The methods are provisional in the sense that they are reasonable approaches,
# though not necessarily the best. They tend, for instance, to favor values of
# V0 near the maximum possible values (determined by m.star > V0*ca.star,)
# hence may suffer from statistical bias of some sort.

K1andV0 <- function(t, ca.star, m.star, Ve){
  # Estimates K1 and V0, assuming negligible binding between bolus onset
  # and the time at which plasma activity declines to half its peak value.
  #
  # Args:
  #   t: time
  #   c.star: plasma activity
  #   m.star: activity in the region of interest
  #   Ve: equilibrium volume of distribution of exchangeable tracer as estimated
  #   from analysis of cerebellum TAC
  #
  # Returns:
  #   a vector of values V0 and K1, respectively.
  temp <- ErrorAndK1vsV0(t, ca.star, m.star, Ve)
  n <- which.min(temp[,"error"])
  return(temp[n,c("V0","K1")])
}

V0K1k2inCerebellum <- function(t, ca.star, m.star){
  # Estimates V0, K1, and k2 for cerebellum by fitting
  # the integral form of a standard compartmental model
  # to data.
  #
  # Args:
  #   t: time
  #   ca.star: plasma activity
  #   m.star: cerebellum activity
  #
  # Returns:
  #   a vector of values V0, K1, and k2.
  temp <- CerebellumAnalysis(t,ca.star,m.star)
  i.min <- which.min(temp[,"error"])
  return(temp[i.min,c("V0","K1","k2")])
}

ErrorAndK1vsV0 <- function(t, ca.star, m.star, Ve){
  # Estimates K1 and associated fitting error for a range of V0 assumptions.
  #
  # Args:
  #   t: time
  #   c.star: plasma activity
  #   m.star: activity in the region of interest
  #   Ve: equilibrium volume of distribution of exchangeable tracer as estimated
  #   from analysis of cerebellum TAC
  #
  # Returns:
  #   a matrix of 3 columns, V0, error, and K1, respectively. The row for
  #   which error is minimum contains the best estimate of V0 and K1.
  rng <- m.star > 0 & ca.star > 0
  # m.star can not exceed V0*ca.star, since m.star=V0*ca.star+me+mb
  max.V0 <- min(na.omit(m.star[rng]/ca.star[rng]))
  min.V0 <- .1*max.V0  # arbitrary, hopefully realistic
  # 50 values of V0 between min.V0 and max.V0
  V0 <- seq(from=min.V0,to=max.V0,length.out=50)
  error <- numeric()
  K1 <- numeric()
  # Compute best values of K1 and associated errors for the 50 values of V0
  for(n in 1:length(V0)){
    ans <- K1andError(t,ca.star,m.star,V0[n],Ve)
    error <- c(error,ans$error)
    K1 <- c(K1,ans$K1)
  }
  return(cbind(V0,error,K1))
}

K1andError <- function(t,ca.star,m.star,V0,Ve){
  # Computes the value of K1 which best fits the data, assuming
  # the given value of V0 is correct.
  #
  # Args:
  #   t: time
  #   ca.star: plasma activity
  #   m.star: activity in the region of interest
  #   V0: an assumed value for apparent vascular volume of distribution
  #   Ve: equilibrium volume of distribution of exchangeable tracer as estimated
  #   from analysis of cerebellum TAC
  #
  # Returns:
  #   a list, list(K1=K1,error=error), of the best K1 and its associated error
  #   assuming the given value of V0 is correct.
  rng <- HalfMaxRange(t,ca.star)
  in.tissue <- (m.star - V0*ca.star)[rng]
  ca.int <- qdint(t[rng],ca.star[rng])
  me.int <- qdint(t[rng],in.tissue)
  x <- ca.int - me.int/Ve
  mdl <- lm(in.tissue ~ 0 + x)
  K1 <- mdl$coefficients[1]
  error <- sqrt(sum(mdl$residuals^2))
  return(list(K1=K1,error=error))
}


HalfMaxRange <- function(t,ca){
  # Finds the times between bolus onset and decline to half its peak value.
  #
  # Args:
  #   t: time
  #  ca: plasma activity
  #
  # Returns:
  #   a logical vector which is TRUE in the appropriate range. 
  t.max <- t[which.max(ca)]
  t.halfmax <- max(t[ t > t.max & ca > 0.5*max(ca)])
  t.onset <- min(t[ t<t.max & ca > .05*max(ca)])
  return( t >= t.onset & t <= t.halfmax)
}

CerebellumAnalysis <- function(t, ca.star, m.star){
  # Estimates K1, k2, and error for a reasonable range of V0 values,
  # in cerebellum.
  #
  # Args:
  #   t: time
  #   ca.star: plasma activity
  #   m.star: cerebellum activity
  #
  # Returns:
  #   a matrix with columns V0, error, K1, and k2. The row with
  #   smallest error constains best estimates of V0, K1, and k2.
  int.ca <- qdint(t,ca.star)
  rng <- m.star > 0 & ca.star > 0
  # m.star can not exceed V0*ca.star, since m.star=V0*ca.star+me+mb
  max.V0 <- min(na.omit(m.star[rng]/ca.star[rng]))
  min.V0 <- .1*max.V0  # arbitrary, hopefully realistic
  # 50 values of V0 between min.V0 and max.V0
  V0 <- seq(from=min.V0,to=max.V0,length.out=50)
  K1 <- numeric()
  k2 <- numeric()
  error <- numeric()
  for(n in 1:length(V0)){
    me <- m.star-V0[n]*ca.star
    int.me <- qdint(t,me)
    mdl <- lm(me ~ 0 + int.ca + int.me)
    K1 <- c(K1,as.numeric(mdl$coefficients["int.ca"]))
    k2 <- c(k2,-as.numeric(mdl$coefficients["int.me"]))
    error <- c(error, sqrt(sum(mdl$residuals^2)))
  }
  return(cbind(V0,error,K1,k2))
}

# Originally From matlab_utils
# Quick and dirty numerical integration of irregularly sampled
# data. Uses trapezoidal approximation.
qdint <- function(t,y){
  n <- length(y)
  temp <- diff(t)*(y[1:(n-1)]+y[2:n])/2
  return(c(0,cumsum(temp)))
}