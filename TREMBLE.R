TREMBLE <- function(path, ab, region, debug=FALSE) {
  # path is the path to the directory containing the scans
  # ab is the prefix of the scans (a or b)
  # region is the brain region to look at
  hsa <- MergeTimes(ReadScan(path, c(paste(ab, "HSA", sep="")))) #get HSA data
  lsa <- MergeTimes(ReadScan(path, c(paste(ab, "LSA", sep="")))) #get LSA data
  
  shsa <- AnalyzeScan(hsa, region)
  slsa <- AnalyzeScan(lsa, region)
  
  phsa <- GetRatios(hsa, shsa)
  plsa <- GetRatios(lsa, slsa)
  
  points <- rbind(phsa, plsa)
  rownames(points) <- c("HSA", "LSA")
  
  return(list(
    subject=basename(path),
    type=ab,
    region=region,
    result=points))
}

GetRatios <- function(scan, analysis) {
  # Args:
  #  scan -- output of MergeTimes
  #  analysis -- the return value from AnalyzeScan
  #
  # Returns a two-wide vector containing pB and mbmax/SA
  imaxtracerintissue <- which.max(analysis$me + analysis$mb)
  idx <- analysis$t >= analysis$t[imaxtracerintissue]
  
  imaxmb <- which.max(analysis$mb[idx]) + imaxtracerintissue
  #TODO: improve the max-finding method
  
  
  pB <- analysis$mb[imaxmb]/analysis$me[imaxmb]
  mbmax.over.SA <- analysis$mb[imaxmb]/scan$SA
  
  return(c(pB=pB, mbmax.over.SA=mbmax.over.SA))
}

AnalyzeScan <- function(scan,roi){
  # Estimates V0, K1, k2, me, and mb for the given scan and region.
  # Also returns the time scale, t, and the derivative of tracer in tissue
  #
  # Args:
  #   scan -- a list of TACs and SA on a common time base, 
  #           as given by MergeTimes(ReadScan...)
  #   roi -- the region of interest, a string, 'LCN', 'RCN', 'LPu', 'RPu'
  # Returns a list containing
  #   V0, K1, k2 -- for the region of interest
  #   me -- activity of exchangable tracer as estimated by 
  #         numerical solution of BBB equation, a vector
  #   mb -- activity of bound tracer, tracer in tissue minus me, a vector
  #   d.dt -- estimated derivative of tracer in tissue, a vector
  #   t  -- time base of me, mb, and d.dt
  
  # Estimate V0, K1, and k2 for cerebellum
  Cb.constants <- V0K1k2inCerebellum(scan$Cb[,1], scan$ca[,2], scan$Cb[,2])
  # Using Ve = K1/k2 for cerebellum, estimate V0, K1, and k2 for the region
  Ve <- Cb.constants["K1"] / Cb.constants["k2"]
  roi.constants <- GetConstants(scan[[roi]][,1], 
                                scan$ca[,2], scan[[roi]][,2], Ve)
  # Smooth tracer in tissue, m*-V0ca*.
  rough.tracer.in.tissue <- scan[[roi]][,2] - roi.constants["V0"]*scan$ca[,2]
  #    NOTE: SmoothTracerInTissue returns a list, whose first item is smoothed
  #    input. The next item is a deprecated derivative, ignored here.
  tracer.in.tissue <- SmoothTracerInTissue(scan[[roi]][,1],
                                           rough.tracer.in.tissue)[[1]]
  # Estimate the derivative of m*-V0ca*.
  deriv <- SmoothCurve(scan[[roi]][,1], scan[[roi]][,2],
                       scan$ca[,2], roi.constants["V0"])
  # Put smoothed tracer in tissue on the same time base as its derivative
  tracer.in.tissue <- approx(tracer.in.tissue[,1],
                             tracer.in.tissue[,2],
                             xout = deriv[,1])$y
  # Put plasma on the same time base
  ca <- approx(scan$ca[,1], scan$ca[,2], xout = deriv[,1])$y
  # Calculate me
  me <- (roi.constants["K1"]*ca - deriv[,2])/roi.constants["k2"]
  # Calculate mb
  mb <- tracer.in.tissue - me
  return(list(V0 = roi.constants["V0"],
              K1 = roi.constants["K1"],
              k2 = roi.constants["k2"],
              me = me,
              mb = mb,
              t  = deriv[,1],
              d.dt = deriv[,2]))
}

ScanIsPresent <- function(path,condition){
  # Returns TRUE if all data for the given condition, 'aLSA', 'aHSA', 'bLSA',
  # or 'bHSA" is present. This will only work for data sets 1 and 2.
  #
  # Args:
  #   path: a string, path to the data directory
  #   condition: one of 'aHSA', 'aLSA', 'bHSA', or 'bLSA'
  ls.dir <- dir(path)
  present <- ls.dir[grep(condition,ls.dir)]
  has.plasma <- length(grep('plasma',present)) > 0
  has.SA <- length(grep('[[:punct:]]SA',present)) > 0
  has.5.TACs <- length(grep('[[:punct:]]TAC',present)) == 5
  return(has.plasma && has.SA && has.5.TACs)
}