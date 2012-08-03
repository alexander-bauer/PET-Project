RoughTREMBLE <- function(subject, ab, region, plot=FALSE) {
  # subject - String such as "Subject01"
  # ab is the prefix of the scans (a or b)
  # region is the brain region to look at
  raw.path <- "../../securedata/RAC_material_1_Sept_2008/"
  subject.path <- paste(raw.path, subject, sep="")
  
  if(!ScanIsPresent(subject.path, paste(ab, "LSA", sep=""))) {return(NULL)}
  
  s <- Load2008HSAModel(subject, region, paste(ab, "HSA", sep=""))
  
  lsa <- MergeTimes(ReadScan(subject.path,
                             paste(ab, "LSA", sep=""))) #get LSA data
  
  lsa.bbb <- LoessBBB(SolveBBB(lsa$ca[,1], lsa$ca[,2], lsa[[region]][,2],
                               s$V0, s$K1, s$k2))
  
  idx <- lsa.bbb$t > 5 & lsa.bbb$t < 25
  shift <- sum((lsa.bbb$t <= 5), na.rm=TRUE) #count the number of TRUEs
  
  i.max <- which.max(lsa.bbb$mb.loess[idx]) + shift
  
  lsa.mb.max <- lsa.bbb$mb[i.max]/lsa$SA
  lsa.pB <- lsa.bbb$mb[i.max]/lsa.bbb$me[i.max]
  
  hsa.mb.max <- max(s$mb/s$SA)
  hsa.pB <- s$k3/s$k4
  
  pB <- c(hsa.pB, lsa.pB)
  mb.max <- c(hsa.mb.max, lsa.mb.max)
  mdl <- lm(mb.max ~ pB)
  KD <- as.numeric(-mdl$coefficients[2])
  Bmax <- as.numeric(mdl$coefficients[1])
  
  
  if(plot)
  {
    matplot(lsa.bbb$t, cbind(lsa.bbb$me.loess, lsa.bbb$mb.loess),
            type='l', col=c(3, 2))
    abline(v=5, col=2)
    abline(v=25, col=2)
    abline(v=lsa.bbb$t[i.max], col=1)
  }
  
  output <- data.frame(hsa.pB=hsa.pB,
                       hsa.mb.max=hsa.mb.max,
                       lsa.pB=lsa.pB,
                       lsa.mb.max=lsa.mb.max,
                       KD=KD,
                       Bmax=Bmax)
    
  return(output)
}


PlotRoughTREMBLE <- function(a.segments, b.segments, col, main){
  # Given data frames whose rows are RoughTremble outputs, plot
  # the a.segments in color col[1], the b.segments in color col[2]
  
  # Change the column names to permit following rbind
  names(a.segments) <- names(b.segments) <- c("pB", "mb.max", "pB", "mb.max",
                                              "KD", "Bmax")
  
  # Set up the scale and labels
  plot(rbind(a.segments[,1:2],
             a.segments[,3:4],
             b.segments[,1:2],
             b.segments[,3:4]),
       type='n',
       xlab='Binding Potential', ylab='Max mb (pmol/ml)', main=main)
  # Plot the a segments
  segments(x0=a.segments[,1], y0=a.segments[,2],
           x1=a.segments[,3], y1=a.segments[,4],
           col=col[1],lty=1)
  # Plot the b segments
  segments(x0=b.segments[,1], y0=b.segments[,2],
           x1=b.segments[,3], y1=b.segments[,4],
           col=col[2],lty=1)
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
  ca <- SmoothPlasmaTac(scan$ca[,1], scan$ca[,2])
  ca <- approx(ca[,1], ca[,2], xout = deriv[,1])$y
  
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