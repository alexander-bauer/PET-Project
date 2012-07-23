TREMBLE <- function(path, ab, region) {
  # path is the path to the directory containing the scans
  # ab is the prefix of the scans (a or b)
  # region is the brain region to look at
  hsa <- MergeTimes(ReadScan(path, c(paste(ab, "HSA", sep="")))) #get HSA data
  lsa <- MergeTimes(ReadScan(path, c(paste(ab, "LSA", sep="")))) #get LSA data
  
  hsa.Cb.constants <- V0K1k2inCerebellum(hsa$Cb[,1], hsa$ca[,2], hsa$Cb[,2])
  lsa.Cb.constants <- V0K1k2inCerebellum(lsa$Cb[,1], lsa$ca[,2], lsa$Cb[,2])
  
  hsa.constants <- GetConstants(hsa[[region]][,1], hsa$ca[,2], hsa[[region]][,2],
                                hsa.Cb.constants["K1"] / hsa.Cb.constants["k2"])
  lsa.constants <- GetConstants(lsa[[region]][,1], lsa$ca[,2], lsa[[region]][,2],
                                lsa.Cb.constants["K1"] / lsa.Cb.constants["k2"])
  
  hsa.smoothderiv <- SmoothCurve(hsa[[region]][,1], hsa[[region]][,2],
                             hsa$ca[,2], hsa.constants["V0"])
  lsa.smoothderiv <- SmoothCurve(lsa[[region]][,1], lsa[[region]][,2],
                                 lsa$ca[,2], lsa.constants["V0"])
  
  hsa.k1ca <- approx(hsa[[region]][,1], (hsa.constants["K1"] * hsa$ca[,2]),
                     xout=hsa.smoothderiv[,1])$y
  lsa.k1ca <- approx(lsa[[region]][,1], (lsa.constants["K1"] * lsa$ca[,2]),
                     xout=lsa.smoothderiv[,1])$y
  
  hsa.me <- (hsa.k1ca - hsa.smoothderiv[,2])/hsa.constants["k2"]
  lsa.me <- (lsa.k1ca - lsa.smoothderiv[,2])/lsa.constants["k2"]
  
  hsa.mb <- hsa.smoothderiv[,2] - hsa.me
  lsa.mb <- lsa.smoothderiv[,2] - lsa.me
  
  hsa.ipeakmb <- which.max(hsa.mb)
  lsa.ipeakmb <- which.max(lsa.mb)
  
  hsa.pB <- hsa.mb[hsa.ipeakmb]/hsa.me[hsa.ipeakmb]
  lsa.pB <- lsa.mb[lsa.ipeakmb]/lsa.me[lsa.ipeakmb]
  
  hsa.mboverSA <- hsa.mb[hsa.ipeakmb]/hsa[["SA"]]
  lsa.mboverSA <- lsa.mb[lsa.ipeakmb]/lsa[["SA"]]
  
  
  
#   +2. For each scan and one specific region, say LCN
#   +----Estimate V0, K1, and k2 using K1_V0_utils.R
#   +----Using V0, form and differentiate tracer in tissue
#   +----Using K1 and k2, solve for me
#   +----Subtract me from tracer in tissue to form mb
#   +----Find the peak value of mb, and the ratio mb/me at this point. The ratio, me/mb, is called binding potential, pB.
#   +----Divide mb.peak by specific activity. (Converts from counts of radioactive tracer, to total tracer, hot plus cold)
#   3. Plot the two (mb.peak/SA, pB) pairs on a graph and connect them with a line. (In principle, not literally.)
#   ---- The NEGATIVE of the slope of that line is K_D, the Michaelis half-saturation concentration (If I remember correctly. It may be K_M, the half-sat constant.)  
#   ---- The intercept of that line is B'max, which estimates the available receptors. It's called the binding capacity after blockade by inhibitors.
}