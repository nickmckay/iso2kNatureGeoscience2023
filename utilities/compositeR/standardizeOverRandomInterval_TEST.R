# INPUTS FROM spreadPaleoData_TEST.r, simpleBinTS_TEST.r and other

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

ages <- binAges
pdm <- binMatR

scaledPdm <- matrix(NA,nrow = nrow(pdm),ncol = ncol(pdm))

  for(i in 1:(ncol(pdm))){
    
    #get the possible range
    
    goodAges <-   ages[!is.na(pdm[,i])] # The ages for which there is binned data in this record
    
    pStart <- min(goodAges)
    pEnd <- max(goodAges)
    
    if(pEnd-pStart < duration){ #only run if the record in total last long enough to warrant inclusion (50 years)
      warning(paste("column",i,"doesnt have the required duration"))
      scaledPdm[,i] <- NA
      next
    }
    
    #find any places where goodAges has data that is a) larger than the lower two of 0 or the first age and also b) less than the lowest of either 2000 or the last age, minus the acceptable duration. this is going to be most of goodAges. 
    
    startOptions <- goodAges[goodAges>= max(min(c(searchRange,pStart))) & goodAges<=min(c(max(searchRange),pEnd)-duration)]
    
    if(length(startOptions) < 1){
      warning(paste("column",i,"doesnt have an overlap between the required duration and the searchRange"))
      scaledPdm[,i] <- NA
      next
    }
    
    nVal <- 0
    tt <- 0
    posStarts <- sample(startOptions) # samples startOptions, without replacement. Basically a shuffle? I think here we're looking across the whole vector, using any point that is less than the set miniumm duration (50 years) then getting that little duration chunk and shuffling it into place. 
    for(ps in posStarts){
      iStart <- ps
      iEnd <- iStart+duration
      
      #check nvalues
      pass <- which(ages>=iStart & ages<iEnd) #find which section of the binAges vector is covered by this particular record, so we can compare this record with the rest of the recrods in these bins
      nVal <- sum(is.finite(removeConsecutiveDuplicates(pdm[pass,i]))) #passes the spread data to the rCD function
      
      if(nVal >= minN){
        break #move on
      }
    }
    
    if(nVal < minN){ #here we have looked over all of the possible start dates in the particular dataset, to see if any have a long enough period of data! This may be what I need to change. 
      warning(paste("cant find an interval in column",i,"that has enough (at least minN (",minN,") observations"))
      scaledPdm[,i] <- NA
      next
    }
    
    #subset the matrix by the interval
    spdm <- as.matrix(pdm[pass,i])
    
    sm <- mean(spdm,na.rm = TRUE) #find mean
    ss <- sd(removeConsecutiveDuplicates(spdm),na.rm = TRUE) #find SD
    #scale the matrix
    if(normalizeVariance){
      scaledPdm[,i] <- scale(pdm[,i],center = sm,scale = ss)
    }else{
      scaledPdm[,i] <- scale(pdm[,i],center = sm,scale = FALSE)
    }
  }
    

