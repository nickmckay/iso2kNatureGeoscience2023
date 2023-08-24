#composite with uncertainty
compositeEnsembles <- function(fTS,binvec,spread = TRUE,stanFun = standardizeMeanIteratively,ageVar = "age",gaussianizeInput = FALSE,alignInterpDirection = TRUE,binFun = sampleEnsembleThenBinTs,minN = 8,...){
  binAges <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

  binMatR <- as.matrix(purrr::map_dfc(fTS,binFun,binvec,ageVar = ageVar,spread = spread,gaussianizeInput = gaussianizeInput,alignInterpDirection = alignInterpDirection))
  binMatR[is.nan(binMatR)] <- NA

  compMat <- stanFun(ages = binAges,pdm = binMatR,minN = minN,...) #one column for every dataset that went in; one row for every bin

  #which records contributed?
  gm <- is.finite(compMat) # logical; did that record (i.e. column) contribute data to a particular bin (i.e. row)
  comp <- rowMeans(compMat,na.rm = TRUE)
  count <- rowSums(gm)
  contributed <- gm
  colnames(contributed) <- lipdR::pullTsVariable(fTS,variable = "paleoData_TSid")

  return(list(composite = comp, count = count, contributed = contributed))

}






