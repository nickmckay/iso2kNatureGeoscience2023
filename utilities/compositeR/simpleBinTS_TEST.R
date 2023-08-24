
# INPUTS FROM spreadPaleoData_TEST.r and other

binvec <-  seq(1, to = 2000, by = 30)
binYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)])) #i.e. from 1.5 to 1999.5, by 1
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


simpleBinTs <- list()

for(i in 1:(length(ts))) {
  
  #from spreadPaleoData.R
    age <- sp[[i]]$spreadAge
    vals <- sp[[i]]$spreadVal
  
  #gaussianize?
  if(gaussianizeInput){
    vals <- geoChronR::gaussianize(vals)
  }
  
  
  if(alignInterpDirection){
    #check for direction
    din <- names(ts[[i]])[stringr::str_detect("_interpDirection",string = names(ts))]
    di <- unlist(magrittr::extract(ts[[i]],din))
    
    sin <- names(ts[[i]])[stringr::str_detect("_scope",string = names(ts[[i]]))]
    si <- unlist(magrittr::extract(ts[[i]],sin))
    
    di <- di[grepl(pattern = "climate",x = si)]
    
    if(length(di)>0){
      if(all(grepl(di,pattern = "negative",ignore.case = TRUE))){
        vals <- sp$spreadVal * -1
      }
    }
  }
  
  bd <- geoChronR::bin(age,vals,binvec = binvec)[,2] #binvec is the bin EDGES
  simpleBinTs[[i]] <- bd
  
}