
## ---------------------------------------------------------------------------
##
## Script name: Figure 2
##
## Purpose of script: Recreate Figure 2 from Konecky et al (2023) Nature Geoscience
##
## Script authors: Georgina Falster & Jonathan Tyler
##
## Date create/updated: 2023-08-23
##
## Email: georgina.falster@anu.edu.au / georgina.falster@gmail.com

## ----------------------------------------------------------------------------
## 1. Perform PCA on Iso2k records (for years 850-1840)
## 2. Calculate trends in individual Iso2k records over the same time period
## 3. Plot the results
## ----------------------------------------------------------------------------

## ----------------------------------------------------------------------------
## You will need to set your working directory to the folder containing the auxiliary data and functions
## ----------------------------------------------------------------------------

# =============================================================================
# set display options
# =============================================================================

options(scipen = 10, digits = 4) # I'm not a fan of scientific notation

# =============================================================================
# load required packages
# =============================================================================

library(geoChronR) # instructions for installation https://nickmckay.github.io/GeoChronR/
library(lipdR) # instructions for installation https://github.com/nickmckay/lipdR
library(magrittr)
library(rgeos)
library(rgdal)
library(data.table)
library(vegan)
library(tidyverse)
library(RSpectra)
library(patchwork)

# =============================================================================
# Load required functions
# =============================================================================

# -------------------------------------------------------------
# you will need the folder 'functions' and then edit the filepaths here so R can read them in
# -------------------------------------------------------------

source("utilities/R/dineof.r")
source("utilities/R/iso2kPCA_fun.r")

# =============================================================================
# Read in v1 of the Iso2k database
# =============================================================================

# -------------------------------------------------------------
# download to a location of your choice (edit path passed to 'destfile')
# -------------------------------------------------------------

download.file("https://lipdverse.org/iso2k/1_0_0/iso2k1_0_0.RData",
              destfile = "data/iso2k1_0_0.RData", method = 'curl')

load("data/iso2k1_0_0.RData")

rm(D, TS)

# =============================================================================
# Read in PAGES 2k Global Mean Temperature (GMT) reconstruction
# =============================================================================

# Paper here: https://www.nature.com/articles/s41561-019-0400-0
# Data here: https://www.ncei.noaa.gov/pub/data/paleo/pages2k/neukom2019temp/recons/

# If you have trouble with curl, then you can download the GMT reconstructions directly from the link above

# -------------------------------------------------------------
# download to a location of your choice (edit path passed to 'destfile')
# -------------------------------------------------------------

download.file(url = "https://www.ncei.noaa.gov/pub/data/paleo/pages2k/neukom2019temp/recons/Full_ensemble_median_and_95pct_range.txt",
              destfile = "data/Full_ensemble_median_and_95pct_range.txt", method = 'curl')

# -------------------------------------------------------------
# read in the downloaded GMST reconstruction
# -------------------------------------------------------------

p2k <- read.table("data/Full_ensemble_median_and_95pct_range.txt",
                  sep = "\t", header = TRUE)

p2k <- p2k[, -c(2:6)] %>%
  setNames(c("year", "median_31", "lwr_31", "upr_31"))

# =============================================================================
# Database filtering
# =============================================================================

# -------------------------------------------------------------
# extract variables for filtering
# -------------------------------------------------------------

variableName <- pullTsVariable(sTS,variable = "paleoData_variableName")
primaryTS <- pullTsVariable(sTS,variable = "paleoData_iso2kPrimaryTimeseries")

# -------------------------------------------------------------
# filter database for primary timeseries
# -------------------------------------------------------------

isd18O <- which(variableName == "d18O" & primaryTS == "TRUE")
isd2H <- which(variableName == "d2H" & primaryTS == "TRUE")
is.iso <- c(isd18O,isd2H)

# ---------------------
# contains all records that have O or H isotopes as the variable and are the primary timeseries
# ---------------------

primary.isoTS <- sTS[is.iso]

# -------------------------------------------------------------
# split by isotope interpretation group
# -------------------------------------------------------------

interpGroup <- pullTsVariable(primary.isoTS,variable = "isotopeInterpretation1_variableGroup")

precIsoTS <- primary.isoTS[which(interpGroup == "P_isotope")]
tempTS <- primary.isoTS[which(interpGroup == "Temperature")]
effMoistTS <- primary.isoTS[which(interpGroup == "EffectiveMoisture")]

# =============================================================================
# Create database summary table
# =============================================================================

int_names <- c("precIso", "temp", "effMoist")

for(int_name in int_names) {

  this_ts <- get(paste0(int_name, "TS"))

  data_summary <- matrix(NA, length(this_ts), 11)
  colnames(data_summary) <- c("firstdate", "lastdate", "recordlength", "resolution", "recordName",
                              "lat", "lon", "archive", "name","interp", "direction")

  for(i in 1:length(this_ts)){
    thisIsoVec <- as.numeric(this_ts[[i]]$paleoData_values)
    thisYearVec <- as.numeric(this_ts[[i]]$year)
    if(length(thisIsoVec) == length(thisYearVec))
    {theseVec <- na.omit(cbind(thisIsoVec, thisYearVec))}
    iso <- theseVec[,1]
    year <- theseVec[,2]
    data_summary[i,1] <- min(year)
    data_summary[i,2] <- max(year)
    data_summary[i,3] <- max(year) - min(year)
    data_summary[i,4] <- (max(year) - min(year)) / length(iso)
    data_summary[i,5] <- this_ts[[i]]$paleoData_TSid
    data_summary[i,6] <- this_ts[[i]]$geo_latitude
    data_summary[i,7] <- this_ts[[i]]$geo_longitude
    data_summary[i,8] <- this_ts[[i]]$archiveType
    data_summary[i,9] <- this_ts[[i]]$dataSetName
    data_summary[i,10] <- this_ts[[i]]$isotopeInterpretation1_variable
    data_summary[i,11] <- if
    (is.null(this_ts[[i]]$isotopeInterpretation1_direction)) {
      NA
    } else{
      this_ts[[i]]$isotopeInterpretation1_direction
    }

  }

  data_summary <- as.data.frame(data_summary) %>%
    mutate_at(c("lat", "lon", "recordlength"), as.numeric)

  # ---------------------
  # make new object
  # ---------------------

  assign(paste0(int_name, "_df"), data_summary)

  rm(data_summary, this_ts)

}

# =============================================================================
# Bin data to fixed range
# =============================================================================

# -------------------------------------------------------------
# common parameters
# -------------------------------------------------------------

XMIN = 850
XMAX = 1840
RES = 30
YEAR <- seq(XMIN, XMAX, RES) #bin edges
binYear <- seq((XMIN + RES/2), (XMAX - RES/2), RES) #bin middles

# -------------------------------------------------------------
# go through each isotope interpretation group
# -------------------------------------------------------------

# there is some legacy code here which isn't actually relevant any more

for(int_name in int_names) {

  this_ts <- get(paste0(int_name, "TS"))

  isoInterpolate <- matrix(NA, length(YEAR), length(this_ts))
  isoBin <- matrix(NA, length(YEAR)-1, length(this_ts))
  isoOriginal <- vector("list")
  isoOriginalDate <- vector("list")
  isoIntName <- matrix(NA, 1, length(this_ts))
  isoIntTSid <- matrix(NA, 1, length(this_ts))


  for(i in 1:length(this_ts)){
    thisIsoVec <- as.numeric(this_ts[[i]]$paleoData_values)
    thisYearVec <- as.numeric(this_ts[[i]]$year)
    thisRecordName <- this_ts[[i]]$paleoData_TSid
    TSid <- this_ts[[i]]$paleoData_TSid

    if(length(thisIsoVec) == length(thisYearVec)){
      theseVec <- na.omit(cbind(thisIsoVec, thisYearVec))
      iso <- theseVec[,1]
      year <- theseVec[,2]
      if(nrow(theseVec) > 0 && length(unique(year)) > 1){
        isoInterpolate[,i] <- approx(year, iso, YEAR)$y
        isoBin[,i] <- as.vector(tapply(iso, cut(year, YEAR), mean))
        isoOriginal[[i]] <- iso
        isoOriginalDate[[i]] <- year
        isoIntName[i] <- thisRecordName
      }
    }
  }

  ## screen datasets for just those with >80% data coverage

  # First assign isoInt to be either the interpolated or binned data (binned, in this case)
  isoInt <- isoBin

  isoInt.keepCol <- matrix(NA, length(this_ts), 1)
  isoOriginal.keep <- vector("list")
  isoOriginalDate.keep <- vector("list")

  for(i in 1:ncol(isoInt)){
    if( (length(na.omit(isoInt[,i])) / length(isoInt[,i])) >= .85
        && length(unique(isoInt[,i])) > 2
        && mean(na.omit(isoInt[,i])) > -1000)
    {
      isoInt.keepCol[i,1] <- i}
  }

  for(i in 1:length(isoInt.keepCol)){
    if(is.numeric(isoInt.keepCol[i])){
      isoOriginal.keep[[i]] <- isoOriginal[[i]]
      isoOriginalDate.keep[[i]] <- isoOriginalDate[[i]]
    }
  }


  isoInt.keepCol.onlykeep <- c(na.omit(isoInt.keepCol[,1]))

  isoInt.keep <- isoInt[, isoInt.keepCol.onlykeep]
  isoIntName.keep <- isoIntName[, isoInt.keepCol.onlykeep]

  # ---------------------
  # make new objects
  # ---------------------

  assign(paste0(int_name, ".isoInt"), isoInt.keep)
  assign(paste0(int_name, ".isoIntName"), isoIntName.keep)

  rm(this_ts, isoInterpolate, isoOriginal, isoOriginalDate, isoIntName, isoIntTSid, thisIsoVec, thisYearVec, TSid, theseVec, iso, year, isoInt, isoBin,
     isoInt.keepCol, isoOriginal.keep, isoOriginalDate.keep, isoInt.keepCol.onlykeep, isoInt.keep, isoIntName.keep)

}

# =============================================================================
# Find site means if necessary
# =============================================================================

# -------------------------------------------------------------
# this is actually not necessary here, but would be useful for sites with multiple records. here we just go through the motions...
# -------------------------------------------------------------

for(int_name in int_names) {

  isoIntName <- get(paste0(int_name, ".isoIntName"))
  isoInt <- get(paste0(int_name, ".isoInt"))

  uniqueNames <- unique(isoIntName)
  isoSiteMeans <- matrix(NA, nrow(isoInt), length(uniqueNames))

  for(i in 1:length(uniqueNames)){
    thisSite <- which(isoIntName == uniqueNames[i])
    thisSiteDat <- data.frame(isoInt[,thisSite])
    thisSiteMean <- apply(thisSiteDat, 1, mean)
    isoSiteMeans[,i] <- thisSiteMean
  }

  # ---------------------
  # make new objects
  # ---------------------

  assign(paste0(int_name, ".uniqueNames"), uniqueNames)
  assign(paste0(int_name, ".isoSiteMeans"), isoSiteMeans)

}

# =============================================================================
# Perform Principal Component Analysis
# =============================================================================

for(int_name in int_names) {

  isoSiteMeans <- get(paste0(int_name, ".isoSiteMeans"))
  uniqueNames <- get(paste0(int_name, ".uniqueNames"))
  data_summary <- get(paste0(int_name, "_df"))

  iso2kPCA(data=isoSiteMeans, recordName=uniqueNames, data_summary=data_summary)

  PCA_results.means <- PCA_results
  PCA_ts <- PCA_timeseries
  PCA_out.means <- pca_out

  # ---------------------
  # make new objects
  # ---------------------

  assign(paste0(int_name, ".PCA_results.means"), PCA_results.means)
  assign(paste0(int_name, ".PCA_ts"), PCA_ts)
  assign(paste0(int_name, ".PCA_out.means"), PCA_out.means)
}

# -------------------------------------------------------------
# at this point, if you like, you can e.g. inspect the broken stick plots. For example:
# -------------------------------------------------------------

screeplot(precIso.PCA_out.means, bstick=TRUE, type="lines", main="")

# =============================================================================
# Prepare to plot the PC1s, with the PAGES 2k Global Mean Temperature (GMT) reconstruction
# =============================================================================

# -------------------------------------------------------------
# some preparation for plotting
# -------------------------------------------------------------

binvec <-  seq(XMIN, XMAX, by = RES)
binYear <- seq((XMIN + RES/2), (XMAX - RES/2), RES) # bin middles, for plotting

for(int_name in int_names) {

 these_pcs <- get(paste0(int_name, ".PCA_ts")) %>%
   as.data.frame()

 these_pcs$binMid <- binYear
 these_pcs$binStart <- binvec[1:(length(binvec)-1)]

 # ---------------------
 # write out same object
 # ---------------------

 assign(paste0(int_name, ".PCA_ts"), these_pcs)

}

# =============================================================================
# Plot the PC1s, with the PAGES 2k Global Mean Temperature (GMT) reconstruction (i.e., Fig. 2a)
# =============================================================================

# -------------------------------------------------------------
# Note that ggplot does NOT like different y axes (because one of the creators, Hadley Wickham, does not 'believe' in dual axes) -
# see replies to https://stackoverflow.com/questions/3099219/ggplot-with-2-y-axes-on-each-side-and-different-scales
#
# So, we have to do a bit of mucking around to make this work nicely
# -------------------------------------------------------------

# -------------------------------------------------------------
# Also note that there is a small element of randomness in the PCA, and that the sign of a PC is arbitrary.
# So you may need to flip one or more of your PCs so make them correlate positively with the GMST reconstruction
# -------------------------------------------------------------

allPCaPlot <- ggplot() +
  # show PAGES 2k temperature reconstruction
  geom_ribbon(data = p2k, aes(x = year, ymin = lwr_31*5.5+1.8,
                              ymax = upr_31*5.5+1.8), fill = "grey") +
  # PC1: Effective moisture
  geom_line(data = effMoist.PCA_ts, aes(x = binMid, y = PC1*-1),colour = "seagreen4",  size = 1) +
  # PC1: P_isotope
  geom_line(data = precIso.PCA_ts, aes(x = binMid, y = PC1*-1), colour = "mediumpurple4",  size = 1) +
  # PC1: Temperature
  geom_line(data = temp.PCA_ts, aes(x = binMid, y = PC1),colour = "indianred4",   size = 1) +
  # axis scaling etc
  scale_x_continuous(name = "Year CE",oob = scales::squish, breaks = seq(900, 1900, 100),
                     minor_breaks = seq(900, 1900, 50), limits = c(XMIN, XMAX), expand = c(0, 0)) +
  scale_y_continuous(name = "PC1",limits = c(-2.6, 2.6), oob = scales::squish, sec.axis = sec_axis(~./5.5-0.327, name = "GMST anomaly"))+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 10, colour = "black"),
    plot.margin = unit(c(0, 1, 0.5, 1), "cm"),
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_line(colour = "grey75"),
    panel.grid.minor.x = element_line(colour = "grey92"),
    axis.title.y.right = element_text(colour = "grey"),
    axis.line.y.right = element_line(colour = "grey")
  )

allPCaPlot

###############################################################################
###############################################################################

# From here, we are doing the necessaries for the maps in Figure 2 (Fig. 2b-d)

###############################################################################
###############################################################################

# =============================================================================
# Smooth over the data frames with the PCA loadings (for later matching)
# =============================================================================

for(int_name in int_names) {

  this_df <- get(paste0(int_name, ".PCA_results.means")) %>%
    dplyr::rename(uID = recordName)

  assign(paste0(int_name, "_loadings"), this_df)

  rm(this_df)
}

# =============================================================================
# Now we actaully go back to the full database, calculate trends, then match with records used in the PCA
# =============================================================================

# I recognise that this looks wildly inefficient, but this figure combines what was originally two separate analyses
# In the interests of reproducibility, I am therefore showing exactly how this was done originally, hence there is quite a bit of unnecessary work

# Just persevere, and you will eventually come to some nice plots

# =============================================================================
# Some database filtering
# =============================================================================

variableName <- pullTsVariable(sTS,variable = "paleoData_variableName")
primaryTS <- pullTsVariable(sTS,variable = "paleoData_iso2kPrimaryTimeseries")
units <- pullTsVariable(sTS,variable = "paleoData_units")

isd18O <- which(variableName == "d18O" & primaryTS == "TRUE")
isd2H <- which(variableName == "d2H" & primaryTS == "TRUE")

d18O.ts <- sTS[isd18O]
d2H.ts <- sTS[isd2H]

# -------------------------------------------------------------
# further filtering, to only keep sufficiently long datasets
# -------------------------------------------------------------

startYear <- 850
endYear <- 1840
binSize <- 50
nptsReq <- 20
nptsEndBins <- 1

d18O.filt <- lapply(1:(length(d18O.ts)), function(i) {
  yearVec <- d18O.ts[[i]]$year
  intYears <- yearVec[yearVec >= startYear & yearVec <= endYear]
  startBin <- intYears[intYears >= startYear & intYears <= (startYear + binSize)]
  endBin <- intYears[intYears <= endYear & intYears >= (endYear - binSize)]
  if (length(startBin) >= nptsEndBins &&
      length(endBin) >= nptsEndBins &&
      length(intYears) >= nptsReq) {
    d18O.ts[[i]]
  } else {
    NA
  }
})

d2H.filt <- lapply(1:(length(d2H.ts)), function(i) {
  yearVec <- d2H.ts[[i]]$year
  intYears <- yearVec[yearVec >= startYear & yearVec <= endYear]
  startBin <- intYears[intYears >= startYear & intYears <= (startYear + binSize)]
  endBin <- intYears[intYears <= endYear & intYears >= (endYear - binSize)]
  if (length(startBin) >= nptsEndBins &&
      length(endBin) >= nptsEndBins &&
      length(intYears) >= nptsReq) {
    d2H.ts[[i]]
  } else {
    NA
  }
})

d18O.filt <- Filter(Negate(anyNA), d18O.filt)
d2H.filt <- Filter(Negate(anyNA), d2H.filt)

# =============================================================================
# Like we did above, make summary data frames
# =============================================================================

for(varname in c("d18O", "d2H")) {

  ts_filt <- get(paste0(varname, ".filt"))

  df <- matrix(NA, length(ts_filt), 11)
  colnames(df) <- c("uID", "elevation", "lat", "lon", "interp", "interpDir", "infMat", "res", "dur", "archive", "var")

  for(i in 1:length(ts_filt)){
    df[i,1] <- ts_filt[[i]]$paleoData_TSid
    df[i,2] <- if
    (is.null(ts_filt[[i]]$geo_elevation)) {
      NA
    } else {
      ts_filt[[i]]$geo_elevation
    }
    df[i,3] <- ts_filt[[i]]$geo_latitude
    df[i,4] <- ts_filt[[i]]$geo_longitude
    df[i,5] <- if
    (is.null(ts_filt[[i]]$isotopeInterpretation1_variableGroup)) {
      NA
    } else {
      ts_filt[[i]]$isotopeInterpretation1_variableGroup
    }
    df[i,6] <- if
    (is.null(ts_filt[[i]]$isotopeInterpretation1_direction)) {
      NA
    } else {
      ts_filt[[i]]$isotopeInterpretation1_direction
    }
    df[i,7] <- if
    (is.null(ts_filt[[i]]$paleoData_inferredMaterialGroup)) {
      NA
    } else {
      ts_filt[[i]]$paleoData_inferredMaterialGroup
    }
    df[i,8] <- (max(ts_filt[[i]]$year[ts_filt[[i]]$year <= endYear]) -
                        min(ts_filt[[i]]$year[ts_filt[[i]]$year >= startYear]))/length(ts_filt[[i]]$year[ts_filt[[i]]$year >= startYear & ts_filt[[i]]$year <= endYear])
    df[i,9] <- (max(ts_filt[[i]]$year[ts_filt[[i]]$year <= endYear]) -
                        min(ts_filt[[i]]$year[ts_filt[[i]]$year >= startYear]))
    df[i,10] <- ts_filt[[i]]$archiveType
    df[i,11] <- ts_filt[[i]]$paleoData_variableName
  }

  df <- data.frame(df) %>%
    mutate_at(c("elevation", "lat", "lon", "res", "dur"), as.numeric)

  # ---------------------
  # make a new object
  # ---------------------

  assign(paste0(varname, ".sum"), df)

  rm(df, ts_filt)
}

# =============================================================================
# Calculate the trends
# =============================================================================

divideD = 8
pval = 0.05

d18O.slope <- list()
d18O.r2 <- list()
d18O.p <- list()

d2H.slope <- list()
d2H.r2 <- list()
d2H.p <- list()

for(i in 1:length(d18O.filt)){
  dataVec <- as.numeric(d18O.filt[[i]]$paleoData_values)
  yearVec <- as.numeric(d18O.filt[[i]]$year)
  comVec <- as.data.frame(na.omit(cbind(dataVec, yearVec)))
  toKeep <- subset(comVec, yearVec >= startYear & yearVec <= endYear)
  if (nrow(toKeep) >= nptsReq &&
      length(na.omit(toKeep$dataVec)) == length(na.omit(toKeep$yearVec))) {
    mod <- lm(toKeep$dataVec ~ toKeep$yearVec)
    d18O.slope[i] <- mod$coefficients[2]
    d18O.r2[i] <- summary(mod)$r.squared
    d18O.p[i] <- summary(mod)$coefficients[, 4][2]
  }
  else {
    d18O.slope[[i]] <- NA
    d18O.r2[[i]] <- NA
    d18O.p[[i]] <-NA
  }
}

for(i in 1:length(d2H.filt)){
  dataVec <- as.numeric(d2H.filt[[i]]$paleoData_values)
  yearVec <- as.numeric(d2H.filt[[i]]$year)
  comVec <- as.data.frame(na.omit(cbind(dataVec, yearVec)))
  toKeep <- subset(comVec, yearVec >= startYear & yearVec <= endYear)
  if (nrow(toKeep) >= nptsReq &&
      length(na.omit(toKeep$dataVec)) == length(na.omit(toKeep$yearVec))) {
    mod <- lm(toKeep$dataVec ~ toKeep$yearVec)
    d2H.slope[i] <- mod$coefficients[2]
    d2H.r2[i] <- summary(mod)$r.squared
    d2H.p[i] <- summary(mod)$coefficients[, 4][2]
  }
  else {
    d2H.slope[[i]] <- NA
    d2H.r2[[i]] <- NA
    d2H.p[[i]] <-NA
  }
}

# =============================================================================
# Arrange and combine the trend dataframes
# =============================================================================

# -------------------------------------------------------------
# Unlist the results, and add to the summary data frames. For d2H records, divide by 8.
# -------------------------------------------------------------

d18O.sum$slope <- unlist(d18O.slope)
d18O.sum$r2 <- unlist(d18O.r2)
d18O.sum$pval <- unlist(d18O.p)

d2H.sum$slope <- unlist(d2H.slope)
d2H.sum$r2 <- unlist(d2H.r2)
d2H.sum$pval <- unlist(d2H.p)

trends.sum <- rbind(d18O.sum, d2H.sum)
trends.sum <- na.omit(trends.sum, cols = "slope") # Removes rows for which there is no slope value
trends.sum$slope <- trends.sum$slope*10 # Change to per mil per decade

trends.sum <- trends.sum %>%
  mutate(slope = ifelse(var =="d2H", slope/divideD, slope)) %>%
  mutate(sig = ifelse(trends.sum$pval < 0.05, "sig", "notSig"))

# -------------------------------------------------------------
# make an index of positivity
# -------------------------------------------------------------

trends.sum$slopeDir <- ifelse(trends.sum$slope < 0, "negative trend", "positive trend")

# -------------------------------------------------------------
# finally, split trend dataframes by isotope interpretation
# -------------------------------------------------------------

effMoist_trends <- trends.sum %>%
  filter(interp == "EffectiveMoisture" & sig == "sig")

temp_trends <- trends.sum %>%
  filter(interp == "Temperature" & sig == "sig")

precIso_trends <- trends.sum %>%
  filter(interp == "P_isotope" & sig == "sig")

# =============================================================================
# Combine trend and loading data frames
# =============================================================================

effMoist_comb <- full_join(effMoist_trends, effMoist_loadings, by = "uID") %>%
  select(lat.x, lat.y, lon.x, lon.y, archive.x, archive.y,
                slope, slopeDir, PC1) %>%
  mutate(lat = coalesce(lat.x, lat.y), lon = coalesce(lon.x, lon.y),
         archive = coalesce(archive.x, archive.y)) %>%
  select(lat, lon, archive,slope, slopeDir, PC1)

temp_comb <- full_join(temp_trends, temp_loadings, by = "uID") %>%
  select(lat.x, lat.y, lon.x, lon.y, archive.x, archive.y,
         slope, slopeDir, PC1) %>%
  mutate(lat = coalesce(lat.x, lat.y), lon = coalesce(lon.x, lon.y),
         archive = coalesce(archive.x, archive.y)) %>%
  select(lat, lon, archive,slope, slopeDir, PC1)

precIso_comb <- full_join(precIso_trends, precIso_loadings, by = "uID") %>%
  select(lat.x, lat.y, lon.x, lon.y, archive.x, archive.y,
         slope, slopeDir, PC1) %>%
  mutate(lat = coalesce(lat.x, lat.y), lon = coalesce(lon.x, lon.y),
         archive = coalesce(archive.x, archive.y)) %>%
  select(lat, lon, archive,slope, slopeDir, PC1)

# =============================================================================
# Fill in gaps where necessary
# =============================================================================

effMoist_med <- median(effMoist_comb$slope, na.rm = TRUE)
temp_med <- median(temp_comb$slope, na.rm = TRUE)
precIso_med <- median(precIso_comb$slope, na.rm = TRUE)

effMoist_comb <- effMoist_comb %>%
  mutate(slopeDir = replace_na(slopeDir, "No significant trend"),
         slope = replace_na(slope, effMoist_med),
         PC1 = replace_na(PC1, 0))

temp_comb <- temp_comb %>%
  mutate(slopeDir = replace_na(slopeDir, "No significant trend"),
         slope = replace_na(slope, temp_med),
         PC1 = replace_na(PC1, 0))

precIso_comb <- precIso_comb %>%
  mutate(slopeDir = replace_na(slopeDir, "No significant trend"),
         slope = replace_na(slope, precIso_med),
         PC1 = replace_na(PC1, 0))


# =============================================================================
# Some preparation for mapping (make basemap, project data to Pacific-centred Robinson)
# =============================================================================

proj_string_data <- "+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# -------------------------------------------------------------
# create world base maps at this point, head over to the 'create_world_basemap_PCs.R' script and run the whole thing
# -----------------------------------------
source('utilities/R/create_world_basemap_PCs.R')


for(int_name in int_names) {

  trans <- get(paste0(int_name, "_comb"))
  coordinates(trans) <- c("lon", "lat")
  proj4string(trans) <- CRS("+proj=longlat +datum=WGS84")
  trans <- spTransform(trans, CRS(proj_string_data))

  assign((paste0(int_name, "_rob")), as.data.frame(trans))

}

# =============================================================================
# Set-up for mapping
# =============================================================================

# =============================================================================
# set up shapes
# =============================================================================

shapes <- c("negative trend" = 25, "positive trend" = 24, "Coral" = 0, "MolluskShells" = 13, "GlacierIce" = 5, "GroundIce" = 8,
            "LakeSediment" = 6, "MarineSediment" = 14, "TerrestrialSediment" = 7, "Speleothem" = 9,
            "Sclerosponge" = 3, "Wood" = 2, "No significant trend" = 21)

# -------------------------------------------------------------
# define vector for colourbar (matches the MATLAB red-blue gradient)
# -------------------------------------------------------------

matlab_cols <- c("#0000FF", "#1A1AFF", "#3333FF", "#4D4DFF", "#6666FF", "#8080FF", "#9999FF", "#B3B3FF", "#CDCDFF", "#E6E6FF", "#FFFFFF",
                 "#FFFFFF", "#FFE6E6", "#FFCDCD", "#FFB3B3", "#FF9999", "#FF8080", "#FF6666", "#FF4D4D", "#FF3333", "#FF1A1A", "#FF0000")

matlab_ramp <- colorRampPalette(matlab_cols)

# =============================================================================
# Maps showing each record's loading on PC1, as well as the record's trend across the pre-industrial last millennium (Fig. 2b-d)
# =============================================================================

# -------------------------------------------------------------
# VERY IMPORTANT!! Make sure to flip the PC1s (or not) as necessary to match what you did when plotting the timeseries above
# -------------------------------------------------------------

precIso_map <- world_basemap_rob +
  # Outer points showing loading, and direction of slope
  geom_point(data = precIso_rob, aes(lon, lat, fill = PC1*-1, shape = slopeDir), size = 6) +
  # Inner points showing archive type
  geom_point(data = precIso_rob, aes(lon, lat, pch = archive), size = 3) +
  # Some scalings
  scale_shape_manual(values = shapes) +
  scale_fill_stepsn(colours = matlab_ramp(19), breaks = seq(-0.8, 0.8, 0.1), limits = c(-0.9, 0.9)) +
  ggtitle("P_isotope", subtitle = "30 year bins 850-1840 CE") +
  labs(fill = "PC1") +
  theme(
    legend.text = element_text(size = 12),
    plot.subtitle = element_text(hjust = 0.5),
    plot.title = element_text(hjust = 0.5)) +
  # tinker with this as desired, to show whichever legends you choose
  guides(fill = guide_colourbar(), shape = "none") +
  guides(fill = "none", shape = "none")

effMoist_map <- world_basemap_rob +
  # Outer points showing loading, and direction of slope
  geom_point(data = effMoist_rob, aes(lon, lat, fill = PC1*-1, shape = slopeDir), size = 6) +
  # Inner points showing archive type
  geom_point(data = effMoist_rob, aes(lon, lat, pch = archive), size = 3) +
  # Some scalings
  scale_shape_manual(values = shapes) +
  scale_fill_stepsn(colours = matlab_ramp(19), breaks = seq(-0.8, 0.8, 0.1), limits = c(-0.9, 0.9)) +
  ggtitle("Effective moisture", subtitle = "30 year bins 850-1840 CE") +
  labs(fill = "PC1") +
  theme(
    legend.text = element_text(size = 12),
    plot.subtitle = element_text(hjust = 0.5),
    plot.title = element_text(hjust = 0.5)) +
  # tinker with this as desired, to show whichever legends you choose
  guides(fill = "none", shape = "none")

temp_map <- world_basemap_rob +
  # Outer points showing loading, and direction of slope
  geom_point(data = temp_rob, aes(lon, lat, fill = PC1, shape = slopeDir), size = 6) +
  # Inner points showing archive type
  geom_point(data = temp_rob, aes(lon, lat, pch = archive), size = 3) +
  # Some scalings
  scale_shape_manual(values = shapes) +
  scale_fill_stepsn(colours = matlab_ramp(19), breaks = seq(-0.8, 0.8, 0.1), limits = c(-0.9, 0.9)) +
  ggtitle("Temperature", subtitle = "30 year bins 850-1840 CE") +
  labs(fill = "PC1") +
  theme(
    legend.text = element_text(size = 12),
    plot.subtitle = element_text(hjust = 0.5),
    plot.title = element_text(hjust = 0.5)) +
  # tinker with this as desired, to show whichever legends you choose
  guides(fill = guide_colourbar(), shape = "none") +
  theme(legend.position = "bottom", legend.key.width = unit(0.04, "npc"))

# ---------------------
# show maps together
# ---------------------

precIso_map + effMoist_map + temp_map + plot_layout(ncol = 1)

