
## INPUTS

binvec <-  seq(1, to = 2000, by = 30)
binYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)])) #i.e. from 1.5 to 1999.5, by 1
binAges <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))
ageVar = "year"
alignInterpDirection = FALSE
spread = TRUE
duration = 50
searchRange = c(0,2000)
normalizeVariance = FALSE
spreadBy = abs(mean(diff(binvec)))/10 # From simpleBinTS
gaussianizeInput= FALSE # From simpleBinTS
minAge = -69
ts <- effMoistTS
sp <- spreadPaleoData # FROM spreadPaleoData_TEST.r 
binMatR <- as.matrix(as.data.frame(simpleBinTs))
binMatR[is.nan(binMatR)] <- NA
minN = 2 
thresh = 0.01

ages <- binAges
pdm <- binMatR

start <- scaledPdm
  
  #remove records that failed (this should potentially cause an error in the future)
  start[!is.finite(start)] <- NA
  colsds <- apply(binMatR,2,sd,na.rm =TRUE)
  filledBins <- apply(binMatR,2,function(x) sum(is.finite(x))) ## apply function over matrix columns (2 means apply over columns, not rows)
  
  bad <- which(colSums(is.finite(start))==0 | colsds < 0.1 | filledBins < minN)
  
  if(length(bad) >= (NCOL(start)-2)){
    stop("No good columns after standardization")
  }
  
  
  
  #set up while loop
  delta <- thresh+1
  meanAllRMSE <- 1000 #pick something large
  upmat <- start
  
  stanFun <- list()
  
  while(delta > thresh){
    allRMSE <- matrix(NA,nrow = ncol(start))
    
    rvec <- sample(1:ncol(start)) #randomize the order of columns in the caled paleodata matrix
    
    for(j in rvec){
      #optimize
      #r <- optim(par = mean(upmat[,j],na.rm = TRUE),fn = optimizeMean,td = upmat[,j],palMat = upmat)
      es <- recordRMSE(upmat[,j],upmat) #calculating RMSE of the standardised mean for the column, compared with equivalent rows in the entire pdm
      meanBias <- es$totalBias
      #update values
      upmat[,j] <- upmat[,j]-meanBias #get rid of this bias i.e. bring it down to match the others in the equivalent interval
      allRMSE[j] <- es$totalRMSE
    }
    old <- meanAllRMSE
    meanAllRMSE <- mean(allRMSE,na.rm = TRUE)
    delta <- old-meanAllRMSE
  }
