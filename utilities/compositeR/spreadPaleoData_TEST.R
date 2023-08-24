
####################
## INPUTS
####################

binvec <-  seq(1, to = 2000, by = 30)
binYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)])) #i.e. from 1.5 to 1999.5, by 1
binFun = simpleBinTs
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


spreadPaleoData <- list()

for (i in 1:(length(ts))) {

  age = ts[[i]]$year #from simpleBinTS
  value = ts[[i]]$paleoData_values # From simpleBinTS
  maxGap = as.numeric(quantile(abs(diff(ts[[i]]$year)),probs = .75,na.rm = TRUE)) # From simpleBinTS. Finds 75th percentile of gaps.

  if(length(age)==0){
    #what's happening
    return(list(spreadAge = age,spreadVal = value))
  }
  
  #remove nonfinite ages
  goodAge <- which(is.finite(age))
  if(length(goodAge) < 2){
    return(list(spreadAge = matrix(NA,ncol = length(age)),spreadVal = matrix(NA,ncol = length(age))))
  }
  age <- age[goodAge] # subset age and value to only move forward with good data 
  value <- value[goodAge]
  
  
  hasNas <- FALSE
   #it doesn't handle NAs appropriately, so lets fix that
   if(any(is.na(value))){
     hasNas <- TRUE
     value[is.na(value)] <- -999999
   }
  
  #spread and interpolate
  newAge <- seq(ceiling(min(age)),floor(max(age)),by = spreadBy) #rounds up the minimum, and rounds down the maximum, then makes a sequence between those with step = spreadBy 
  newVals <- pracma::interp1(as.vector(age),as.vector(value),xi = newAge,method = "nearest") #interpolates along newAge, using nearest neighbour interpolation
  
   if(hasNas){
     newVals[which(newVals == -999999)] <- NA
   }
  
  
  #add on to the beginning
  newValsStart <- which(newVals == value[1]) # first point where newVals = the first value in the original data vector
  if(any(diff(newValsStart) != 1)){
    newValsStart <- newValsStart[1:(min(which(diff(newValsStart) != 1)))] # the first value where they're equal
  }
  
  
  f <- min(newAge)-min(c(0,length(newValsStart)-1))*spreadBy
  t <- min(newAge)-spreadBy
  
  if(f < t){
    begAge <- seq(from = f,to = t,by = spreadBy) #closes the gap between f and t, by time step spreadBy
    begVal <- rep(newVals[1],times = length(begAge)) #fills the first few entries of begVal with the first entry of newVals
  }else{
    begAge <- c()
    begVal <- c()
  }
  
  #add on to the end
  end <- length(newAge)
  newValsEnd <- which(newVals == value[length(value)]) # this pulls out the last value of value, and loooks for where that occurs in newVals
  if(any(diff(newValsEnd) != 1)){
    newValsEnd <- newValsEnd[max(which(diff(newValsEnd) != 1)):length(newValsEnd)]
  }
  
  
  f <- max(newAge)+spreadBy # end of the bins, plus one spreader
  t <- max(newAge)+max(c(0,length(newValsEnd)-1))*spreadBy
  if(f < t){
    endAge <- seq(from = f,to = t,by = spreadBy)
    endVal <- rep(newVals[length(newVals)],times = length(endAge))
  }else{
    endAge <- c()
    endVal <- c()
  }
  
  
  
  #append them
  newAgeOut <- c(begAge,newAge,endAge) # sticks those new values onto the ends of the interpolated ages (and values)
  newValsOut <- c(begVal,newVals,endVal)
  
  
  #distance to nearest
  d2n <- purrr::map_dbl(newAgeOut,function(x) min(abs(x-age))) #for each of the new age vector, catches the original chronned sample that's closest and calculates its distance away
  
  #remove values that exceed maxGap
  newValsOut[d2n > maxGap] <- NA #gets rid of values that are too far away from a 'real' data point
  
  #remove values that are too young
  tooYoung <- which(newAgeOut < minAge)
  if(length(tooYoung) > 0){
    newAgeOut <- newAgeOut[-tooYoung]
    newValsOut <- newValsOut[-tooYoung]
  }
  
  
  spreadPaleoData[[i]] <- list(spreadAge = newAgeOut,spreadVal = newValsOut)
}
