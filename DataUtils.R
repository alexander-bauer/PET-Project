path.raw.2008 <- '../securedata/RAC_material_1_Sept_2008'

ReadScan <- function(path,scan){
  # Reads a scan, such as 'aHSA', for a subject, such as 1 for Sept 2008 data.
  # 
  # Args:
  #   subject: path to directory containing raw tac files etc.
  #   scan: 'aHSA', 'bHSA', 'aLSA', or 'bLSA'
  # 
  # Returns:
  #   a labeled list of time activity curves, identified by region, and specific
  #   activity. Each TAC is a 2-column matrix, time and activity, resp. Plasma (ca)
  #   and tissue will be sampled at different times.
  tac.files <- ListTacFiles(path)
  scan.files <- tac.files[grep(scan,tac.files)]
  ans <- list()
  for(n in 1:length(scan.files)){
    temp <- ReadTacFile(scan.files[n])
    bname <- basename(scan.files[n])
    if(Contains('.plasma',bname)){
      ans[['ca']] <- temp[,c(1,3)]
      colnames(ans[['ca']]) <- c("time","activity")
    } else if(substr(bname,nchar(bname)-1,nchar(bname))=='SA'){
      ans[['SA']] <- temp
    } else {
      roi <- tail(strsplit(bname,'_')[[1]],1)
      roi <- substr(roi,1,nchar(roi)-4)
      colnames(temp) <- c("time","activity")
      ans[[roi]] <- temp
    }
  }
  return(ans)
}

MergeTimes <- function(tacs){
  # Puts a list of time-activity curves on a common time scale.
  #
  # Args:
  #   tacs: a list of TACS from function ReadScan
  #
  # Returns:
  #   a similar list with common times and linearly interpolated values,
  #   which takes into account that tissue cannot respond prior to the
  #   onset of bolus.
  ans <- list()
  names <- setdiff(names(tacs),"SA")  # Remove specific activity
  if(!Contains("ca",names))stop("No plasma!") # Abort if there is no plasma
  # Merge and dedupe all times
  common.times <- numeric()
  for(n in 1:length(names)){
    common.times <- c(common.times,tacs[[names[n]]][,1])
  }
  common.times <- sort(common.times)
  common.times <- common.times[!duplicated(common.times)]
  # Find the time of bolus onset, defined as 5% of bolus max.
  # We'll artificially zeroize all activity prior to this time.
  max.ca <- max(tacs$ca[,2])
  bolus.on <- tacs$ca[,2] >= 0.05*max.ca
  t.start <- min(tacs$ca[bolus.on,1])
  # Interpolate all tacs to the common time, with the bolus onset restriction
  # and keeping track of NA positions.
  na.s <- rep(FALSE,length(common.times))
  for(n in 1:length(names)){
    tac <- tacs[[names[n]]]
    # create an artificial lead which is zero before bolus onset
    lead.idx <- tac[,1] <= t.start
    lead.t <- c(tac[lead.idx, 1],t.start)
    lead.t <- lead.t[!duplicated(lead.t)]  # times at which activity is presumed 0
    lead.a <- rep(0,times=length(lead.t))  # activity at those times
    t <- c(lead.t,tac[!lead.idx,1])
    a <- c(lead.a,tac[!lead.idx,2])
    temp <- approx(t,a,xout=common.times)
    na.s <- na.s | is.na(temp$y)
    ans[[names[n]]] <- cbind(temp$x,temp$y)
  }
  for(n in 1:length(names)){
    ans[[n]] <- ans[[n]][!na.s,]
  }
  ans$SA <- tacs$SA
  return(ans)
}

PlotTacs <- function(tacs, regions){
  # Plots the named regions of the given list of tacs on a single panel.
  # 
  # Args:
  #   tacs: a list of tacs such as produced by ReadScan
  #   regions: a vector of names of regions to be plotted
  #
  # Returns:
  #   a matplot of the regions given. NOTE: plasma, if included, will
  #   be scaled by a nominal V0 of 5%.
  temp <- MergeTimes(tacs)
  m <- as.matrix(temp[["ca"]][,1]) # common times
  for(n in 1:length(regions)){
    a <- temp[[regions[n]]][,2]
    if(regions[n]=="ca"){
      a <- 0.05*a
    }
    m <- cbind(m,a)
  }
  matplot(log(m[,1]),m[,2:ncol(m)],type='l',lwd=2,
          xlab='log time',
          ylab='activity',
          main=paste(regions,collapse=", "))
}

ListTacFiles <- function(path) {
  return(dir(path, full.names=TRUE)) 
}

# ListTacFiles <- function(subject){
#   raw.tac.dir <- '../securedata/RAC_material_1_Sept_2008'
#   n <- as.character(subject)
#   if(subject<10){
#     n <- paste('0',n,sep='')
#   }
#   directory <- paste(raw.tac.dir,'/Subject',n,sep='')
#   return(dir(directory,full.names=TRUE))
# }

Contains <- function(pattern,x){
  return(length(grep(pattern,x)) != 0)
}

ReadTacFile <- function(path){
  library(stringr) # for str_trim()
  temp <- readLines(path) # a string (character vector) per line
  if(length(temp)==1 && Contains('.SA',basename(path))){
    # this is a specific activity file with precisely one entry
    return(as.numeric(temp))
  }
  # Split each line on spaces, obtaining a list of "words." Each item
  # in the list represents a line.
  temp <- strsplit(str_trim(temp),' +')
  # If the first line contains only one "word", this word is a line
  # count useful only in the days of yore and FORTRAN. Drop it!
  if(length(temp[[1]]==1)){
    temp <- temp[-1]
  }
  # Convert the list to a matrix of numbers and return it.
  return(t(sapply(temp,as.numeric,simplify="array")))
}

Load2008HSAModel <- function(subject,region,condition){
  # Loads the results of a standard compartmental analysis of HSA
  # data from the Sept 2008 data set. The analyses were performed
  # in Matlab, and give estimates of rate constants, vascular volume
  # of distribution, and bound and exchangeable tracer TACs. The
  # raw data is also returned.
  #
  # Args:
  # subject -- a string, 'Subject01', ..., 'Subject21'
  # region -- region of interest, 'LCN', 'RCN', 'LPu', or 'RPu'
  # condition -- 'aHSA' or 'bHSA'
  #
  # Returns a list containing:
  # K1, k2, k3, k4, and V0 -- clearance, rate constants, and vascular volume of distribution
  # mb -- bound tracer activity
  # me -- exchangeable tracer activity
  # SA -- specific activity (e.g., mb/SA is in pmol/ml)
  # m.star -- activity in region of interest
  # ca.star -- plasma activity
  # t -- common time base for all TACs
  data.dir <- '../securedata/HSA_analyses_for_Sept_2008'
  bname <- paste(paste(subject,condition,region,sep='_'),'.RData',sep="")
  load(paste(data.dir,bname,sep='/'))
  return(s)
}