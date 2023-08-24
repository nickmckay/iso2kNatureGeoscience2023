
## ---------------------------------------------------------------------------
##
## Script name: Figure 1
##
## Purpose of script: Recreate Figure 1 from Konecky et al (2023) Nature Geoscience
##
## Script author: Georgina Falster & Nick McKay
##
## Date create/updated: 2023-08-23
##
## Email: georgina.falster@anu.edu.au / georgina.falster@gmail.com

## ----------------------------------------------------------------------------
## 1. Calculate global composites from Iso2k records
## 2. Calculate correlation of individual records with the composites
## 3. Plot the composites, along with a reconstruction of global mean temperature
## 4. Plot record-composite correlations on maps
## ----------------------------------------------------------------------------

## ----------------------------------------------------------------------------
##
## Note that this figure was made using functions written by Nick McKay;
## These have since been updated, with up-to-date versions available in dev mode from https://github.com/nickmckay/compositeR
##
## In this script I use the original versions used to make the figures shown in the text
##
## You will need to set your working directory to the folder containing the auxiliary data and functions
##
## ----------------------------------------------------------------------------

# =============================================================================
# set display options
# =============================================================================

options(scipen = 10, digits = 4)

# =============================================================================
# load required packages
# =============================================================================
# -------------------------------------------------------------
# Install packages like so: install.packages("forcats")
# -------------------------------------------------------------

library(forcats)
library(maps)
library(rgeos)
library(rgdal)
library(gridExtra)
library(data.table)
library(lipdR) # remotes::install_github("nickmckay/lipdR") - more instructions for installation https://github.com/nickmckay/lipdR
library(data.table)
library(magrittr)
library(tidyverse)
library(foreach)
library(doParallel)
library(geoChronR) # remotes::install_github("nickmckay/geoChronR") instructions for installation https://nickmckay.github.io/GeoChronR/
library(patchwork)
library(pammtools)
library(compositeR) # remotes::install_github("nickmckay/compositeR")

# =============================================================================
# Read in v1 of the Iso2k database
# =============================================================================

# -------------------------------------------------------------
# download to a location of your choice (edit path passed to 'destfile')
# -------------------------------------------------------------

download.file("https://lipdverse.org/iso2k/1_0_0/iso2k1_0_0.RData",
              destfile = "iso2k1_0_0.RData", method = 'curl')

load("iso2k1_0_0.RData")

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
              destfile = "Full_ensemble_median_and_95pct_range.txt", method = 'curl')

# -------------------------------------------------------------
# read in the downloaded GMST reconstruction
# -------------------------------------------------------------

p2k <- read.table("Full_ensemble_median_and_95pct_range.txt",
                               sep = "\t", header = TRUE)

# =============================================================================
# load source files
# =============================================================================

# -------------------------------------------------------------
# you will need to download the folder 'functions' and then edit the filepaths here so R can read them in
# -------------------------------------------------------------

# source("functions/compositeR/bin.R")
# source("functions/compositeR/composite.R")
# source("functions/compositeR/scale.R")
# source("functions/compositeR/standardize.R")
# source("functions/compositeR/uq.R")

# =============================================================================
# Extract variables required for database filtering, and make filters
# =============================================================================

# See Konecky et al (2020) for a complete description of the Iso2k database, including variable descriptions
# https://essd.copernicus.org/articles/12/2261/2020/

# -------------------------------------------------------------
# inspect all names in the sTS object
# -------------------------------------------------------------

allnames <- sapply(sTS,names) %>%
  unlist() %>%
  unique() %>%
  sort()

# -------------------------------------------------------------
# these are the variables we need for filtering
# -------------------------------------------------------------

primary <- pullTsVariable(sTS,variable = "paleoData_iso2kPrimaryTimeseries")
isoInt <- pullTsVariable(sTS,variable = "isotopeInterpretation1_variableGroup")
archive <- pullTsVariable(sTS, variable = "archiveType")
variable <- pullTsVariable(sTS, variable = "paleoData_variableName")

# ---------------------
# if you want to see what the options are for a particular variable
# ---------------------

table(isoInt) %>%
  as.data.frame

# -------------------------------------------------------------
# make initial filters
# -------------------------------------------------------------

effMoist <- which(primary == "TRUE" & isoInt == "EffectiveMoisture")
temp <- which(primary == "TRUE" & isoInt == "Temperature")
precIso <- which(primary == "TRUE" & isoInt == "P_isotope")

# =============================================================================
# Filter timeseries by isotope interpretation group
# =============================================================================

effMoistTS <- sTS[effMoist]
tempTS <- sTS[temp]
precIsoTS <- sTS[precIso]

# ---------------------
# check results
# ---------------------

pullTsVariable(effMoistTS,variable = "paleoData_inferredMaterial") %>%
  unique() %>%
  sort()

# =============================================================================
# Check data density and filter where necessary
# =============================================================================

# -------------------------------------------------------------
# check where there are sufficient entries in both the data and time columns, and filter for TS with more than ten non-NA entries.
# -------------------------------------------------------------

for(TSname in c("effMoistTS", "tempTS", "precIsoTS")) {

  thisTS <- get(TSname)

  checkLength <- map_dbl(thisTS,function(x) sum(!is.na(x$paleoData_values) & !is.na(x$year)))
  checkLength2 <- map_dbl(thisTS,function(x) length(x$paleoData_values))
  thisTS <- thisTS[which(checkLength > 10 & checkLength2 >10)]

  assign(TSname, thisTS)

  rm(thisTS)

}

# =============================================================================
# Set parameters for the composites
# =============================================================================

# -------------------------------------------------------------
# binning parameters
# -------------------------------------------------------------

binvec <-  seq(1, to = 2011, by = 30)
binYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

# -------------------------------------------------------------
# multi-core processing & number of ensemble members
# -------------------------------------------------------------

nens <- 100

registerDoParallel(3)

# =============================================================================
# Calculate composites
# =============================================================================

## Effective moisture

ensOutEffMoist <- foreach(i = 1:nens) %dopar% {
    tc <- compositeEnsembles(effMoistTS,
                             binvec,
                             binFun = simpleBinTs,
                             ageVar = "year",
                             alignInterpDirection = FALSE,
                             spread = TRUE,
                             duration = 90,
                             minN = 3,
                             searchRange = c(0,2011),
                             normalizeVariance = FALSE)

    return(list(composite = tc$composite,count = tc$count, contributed = tc$contributed))
}

## Temperature

ensOutTemp <- foreach(i = 1:nens) %dopar% {
  tc <- compositeEnsembles(tempTS,
                           binvec,
                           binFun = simpleBinTs,
                           ageVar = "year",
                           alignInterpDirection = FALSE,
                           spread = TRUE,
                           duration = 90,
                           minN = 3,
                           searchRange = c(0,2011),
                           normalizeVariance = FALSE)

  return(list(composite = tc$composite,count = tc$count, contributed = tc$contributed))
}

## P_isotope
ensOutPrecIso <- foreach(i = 1:nens) %dopar% {
  tc <- compositeEnsembles(precIsoTS,
                           binvec,
                           binFun = simpleBinTs,
                           ageVar = "year",
                           alignInterpDirection = FALSE,
                           spread = TRUE,
                           duration = 90,
                           minN = 3,
                           searchRange = c(0,2011),
                           normalizeVariance = FALSE)

  return(list(composite = tc$composite,count = tc$count, contributed = tc$contributed))
}

# =============================================================================
# Save composites if you like (edit filepaths as necessary)
# =============================================================================

saveRDS(ensOutEffMoist, "data/EffMoist_compens_nens100_30yr.rds")
saveRDS(ensOutTemp, "data/Temp_compens_nens100_30yr.rds")
saveRDS(ensOutPrecIso, "data/PrecIso_compens_nens100_30yr.rds")

# =============================================================================
# Extract composite timeseries from the outputs
# =============================================================================

effMoistComp_ens <-  as.matrix(purrr::map_dfc(ensOutEffMoist,magrittr::extract,"composite"))
tempComp_ens <-  as.matrix(purrr::map_dfc(ensOutTemp,magrittr::extract,"composite"))
precIsoComp_ens <-  as.matrix(purrr::map_dfc(ensOutPrecIso,magrittr::extract,"composite"))

# =============================================================================
# Determine quantiles for each composite ensemble
# =============================================================================

for(ensname in c("effMoistComp_ens", "tempComp_ens", "precIsoComp_ens")) {

  this_comp <- get(ensname)

  stats <- matrix(NA, nrow(this_comp), 5)
  for (i in 1:(nrow(this_comp))) {

    stats[i, ] <- quantile(this_comp[i, ], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

  }

  # ---------------------
  # wrangle to wreadiness for plotting
  # ---------------------

  stats <- as.data.frame(stats)
  colnames(stats) <- c("lwr2", "lwr1", "med", "upr1", "upr2")
  stats$binStart <- binYears-(unique(diff(binYears))*0.5)
  stats$binEnd <- binYears+(unique(diff(binYears))*0.5)

  # ---------------------
  # make a new object
  # ---------------------

  assign(paste0(gsub("_ens", "", ensname), "_stats"), stats)

  rm(stats)

}

# =============================================================================
# Shift (centre) so the mean of each median is zero
# =============================================================================

for(enssumname in c("effMoistComp_stats", "tempComp_stats", "precIsoComp_stats")) {

  this_comp_sum <- get(enssumname)

  median_mean <- mean(this_comp_sum$med)

  this_comp_sum <- this_comp_sum %>%
    mutate(lwr2 = lwr2-median_mean, lwr1 = lwr1-median_mean, med = med-median_mean, upr1 = upr1-median_mean, upr2 = upr2-median_mean)

  assign(enssumname, this_comp_sum)

  rm(this_comp_sum, median_mean)

}

# =============================================================================
# Plot composites on top of the PAGES 2k global mean temperature reconstruction
# =============================================================================

# -------------------------------------------------------------
# Note that ggplot does NOT like different y axes (because one of the creators, Hadley Wickham, does not 'believe' in dual axes) -
# see replies to https://stackoverflow.com/questions/3099219/ggplot-with-2-y-axes-on-each-side-and-different-scales
#
# So, we have to do a bit of mucking around to make this work nicely
# -------------------------------------------------------------

# ---------------------
# Effective moisture
# ---------------------

effMoistPlot <- ggplot(effMoistComp_stats, aes(x = binStart)) +
  # add PAGES 2k temperature reconstruction
  geom_ribbon(data = p2k, aes(x = Year, ymin = X31.year.filtered.full.ensemble.2.5th.percentile+0.3,
                              ymax = X31.year.filtered.full.ensemble.97.5th.percentile*1+0.3), fill = "grey") +
  # add Iso2k composite
  geom_stepribbon(aes(ymin = lwr2, ymax = upr2),
                  fill = "seagreen2", alpha = 0.7) +
  geom_stepribbon(aes(ymin = lwr1, ymax = upr1),
                  fill = "seagreen4", alpha = 0.7) +
  geom_step(aes(y = med), size = 1) +
  # sort out axis display and scaling
  scale_x_continuous(name = "Year (CE)",oob = scales::squish, breaks = seq(0, 2000, 200),
                     minor_breaks = seq(0, 2000, 100), limits = c(0, 2000), expand = c(0, 0))+
  scale_y_continuous(name = "Δ18O (permil)",limits = c(-0.7, 0.7), oob = scales::squish, sec.axis = sec_axis(~.-0.3, name = "GMST anomaly"))+
  # other aesthetics
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 10, colour = "black"),
    plot.margin = unit(c(0, 1, 0.5, 1), "cm"),
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_line(colour = "grey75"),
    panel.grid.minor.x = element_line(colour = "grey92")
  ) +
  ggtitle("Effective Moisture records")

# ---------------------
# Temperature
# ---------------------

TempPlot <- ggplot(tempComp_stats, aes(x = binStart)) +
  # add PAGES 2k temperature reconstruction
  geom_ribbon(data = p2k, aes(x = Year, ymin = X31.year.filtered.full.ensemble.2.5th.percentile*-0.2-0.06,
                              ymax = X31.year.filtered.full.ensemble.97.5th.percentile*-0.2-0.06), fill = "grey") +
  # add Iso2k composite
  geom_stepribbon(aes(ymin = lwr2, ymax = upr2),
                  fill = "indianred1", alpha = 0.7) +
  geom_stepribbon(aes(ymin = lwr1, ymax = upr1),
                  fill = "indianred4", alpha = 0.7) +
  geom_step(aes(y = med), size = 1) +
  # sort out axis display and scaling
  scale_x_continuous(name = "Year (CE)",oob = scales::squish, breaks = seq(0, 2000, 200),
                     minor_breaks = seq(0, 2000, 100), limits = c(0, 2000), expand = c(0, 0))+
  scale_y_reverse(name = "Δ18O (permil)",limits = c(0.15, -0.15), oob = scales::squish, sec.axis = sec_axis(~./-0.2-0.3, name = "GMST anomaly"))+
  # other aesthetics
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 10, colour = "black"),
    plot.margin = unit(c(0, 1, 0.5, 1), "cm"),
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_line(colour = "grey75"),
    panel.grid.minor.x = element_line(colour = "grey92")
  ) +
  ggtitle("Temperature records")

# ---------------------
# P_isotope
# ---------------------

precIsoPlot <- ggplot(precIsoComp_stats, aes(x = binStart)) +
  # add PAGES 2k temperature reconstruction
  geom_ribbon(data = p2k, aes(x = Year, ymin = X31.year.filtered.full.ensemble.2.5th.percentile+0.3,
                              ymax = X31.year.filtered.full.ensemble.97.5th.percentile+0.3), fill = "grey") +
  # add Iso2k composite
  geom_stepribbon(aes(ymin = lwr2, ymax = upr2),
                  fill = "mediumpurple1", alpha = 0.7) +
  geom_stepribbon(aes(ymin = lwr1, ymax = upr1),
                  fill = "mediumpurple4", alpha = 0.7) +
  geom_step(aes(y = med), size = 1) +
  # sort out axis display and scaling
  scale_x_continuous(name = "Year (CE)",oob = scales::squish, breaks = seq(0, 2000, 200),
                     minor_breaks = seq(0, 2000, 100), limits = c(0, 2000), expand = c(0, 0))+
  scale_y_continuous(name = "Δ18O (permil)",limits = c(-0.7, 0.7), oob = scales::squish, sec.axis = sec_axis(~.-0.3, name = "GMST anomaly"))+
  # other aesthetics
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 10, colour = "black"),
    plot.margin = unit(c(0, 1, 0.5, 1), "cm"),
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_line(colour = "grey75"),
    panel.grid.minor.x = element_line(colour = "grey92")
  ) +
  ggtitle("δprecip records")

# -------------------------------------------------------------
# Show all composites together
# -------------------------------------------------------------

precIsoPlot / effMoistPlot / TempPlot

# =============================================================================
# Calculate correlation of individual contributing records with the composite medians
# =============================================================================

int_names <- c("effMoist", "temp", "precIso")

# -------------------------------------------------------------
# first, make data frames with summary data for each component record
# -------------------------------------------------------------

startYear <- 1
endYear <- 2011

for(int_name in int_names) {

  this_ts <- get(paste0(int_name, "TS"))

  this_df <- matrix(NA, length(this_ts), 7)
  colnames(this_df) <- c("dsn", "lat", "lon", "archive", "infmat", "res", "dur")

  for(i in 1:(length(this_ts))){

    vals <- na.omit(data.frame(year = this_ts[[i]]$year, val = this_ts[[i]]$paleoData_values)) %>%
      filter(year >= startYear & year <= endYear)

    this_df[i, 1] <- this_ts[[i]]$dataSetName
    this_df[i, 2] <- this_ts[[i]]$geo_latitude
    this_df[i, 3] <- this_ts[[i]]$geo_longitude
    this_df[i, 4] <- this_ts[[i]]$archiveType
    this_df[i, 5] <- this_ts[[i]]$paleoData_inferredMaterialGroup
    this_df[i, 6] <- (max(vals$year) - min(vals$year)) / length(vals$year)
    this_df[i, 7] <- max(vals$year) - min(vals$year)
  }

  assign(paste0(int_name, "_df"), as.data.frame(this_df) %>%
           mutate(across(c("lat", "lon", "res", "dur"), as.numeric)))

  rm(this_ts, this_df)

}

# -------------------------------------------------------------
# calculate correlation of each record with its composite, in the interval where it contributed
# -------------------------------------------------------------

for(int_name in int_names) {

  this_ts <- get(paste0(int_name, "TS"))
  this_df <- get(paste0(int_name, "_df"))
  this_comp <- get(paste0(int_name, "Comp_stats"))

  this_df$cor <- NA
  this_df$pval <- NA

  for(i in 1:length(this_ts)){

    vals <- na.omit(data.frame(year = this_ts[[i]]$year, val = this_ts[[i]]$paleoData_values)) %>%
      filter(year >= startYear & year <= endYear)

    # bin the data using the same binning vector as for the composites

    binned_vals <- as.vector(tapply(vals$val, cut(vals$year, breaks = binvec), mean))

    # combine the composite median and the binned record into one data frame, and find correlation

    comDat <- na.omit(data.frame(rec = binned_vals, comp = this_comp$med))

    if (!nrow(comDat) <= 2){

      mod <- cor.test(comDat$rec, comDat$comp)
      this_df$cor[i] <- mod$estimate
      this_df$pval[i] <- mod$p.value
    } else {
      this_df$cor[i] <- NA
      this_df$pval[i] <- NA
    }

  }

  assign(paste0(int_name, "_df"), this_df)

  rm(this_ts, this_df, this_comp, vals, binned_vals, comDat)
}

# -------------------------------------------------------------
# give a correlation of zero where there were insufficient data points (compositeR deals with this via spread&bin)
# -------------------------------------------------------------

for(int_name in int_names){

  this_df <- get(paste0(int_name, "_df"))

  this_df$cor[which(is.na(this_df$cor))] <- 0
  this_df$pval[which(is.na(this_df$pval))] <- 1

  assign(paste0(int_name, "_df"), this_df)

  rm(this_df)
}

# -------------------------------------------------------------
# add some columns for nicer plotting
# -------------------------------------------------------------

for(int_name in int_names) {

  this_df <- get(paste0(int_name, "_df"))

  this_df <- this_df %>%
    mutate(cor_text = ifelse(cor < 0, "Negative", "Positive")) %>%
    mutate(sig = ifelse(pval <= 0.05, "Significant correlation", "No significant correlation")) %>%
    mutate(cor_text = replace(cor_text, sig == "No significant correlation", "No significant correlation"))

  assign(paste0(int_name, "_df"), this_df)

  rm(this_df)

}

# =============================================================================
# Prepare to create maps
# =============================================================================

# -------------------------------------------------------------
# at this point, head over to the 'create_world_basemap.R' script and run the whole thing
# -------------------------------------------------------------

proj_string_data <- "+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

for (int_name in int_names) {

  trans <- get(paste0(int_name, "_df"))
  coordinates(trans) <- c("lon", "lat")
  proj4string(trans) <- CRS("+proj=longlat +datum=WGS84")
  trans <- spTransform(trans, CRS(proj_string_data))
  assign(paste0(int_name, "_rob"), as.data.frame(trans))

}

# -------------------------------------------------------------
# set up some aesthetics
# -------------------------------------------------------------

shapes <- c("Coral" = 24, "Sclerosponge" = 24, "MolluskShells" = 25, "Speleothem" = 24, "MarineSediment" = 21, "LakeSediment" = 21, "Wood" = 22, "GlacierIce" = 23, "GroundIce" = 23, "TerrestrialSediment" = 25)

outlines <- c("Significant correlation" = "black", "No significant correlation" = "grey50")

# =============================================================================
# Maps showing the correlations of each component record with the relevant composite
# =============================================================================

source("utilities/R/create_world_basemap.R")

precIso_map <- world_basemap_rob +
  # add points: shaped by archive type, coloured by correlation, outlined by significance, sized by record duration
  geom_point(data = arrange(precIso_rob, abs(cor)), aes(x = lon, y = lat, fill = cor, colour = sig, size = dur, shape = archive)) +
  # define colours & fills
  scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1,1)) +
  scale_colour_manual(values = outlines) +
  scale_shape_manual(values = shapes) +
  scale_size_binned(breaks = c(0, 500, 1000, 1500, 2000)) +
  labs(size = "Record\nduration\n(years)", fill = "Correlation", colour = "Significant\ncorrelation?", shape = "Archive") +
  # tinker with this as desired, to show whichever legends you choose
  guides(fill = "none", colour = "none", size = "none", shape = "none")

effMoist_map <- world_basemap_rob +
  # add points: shaped by archive type, coloured by correlation, outlined by significance, sized by record duration
  geom_point(data = arrange(effMoist_rob, abs(cor)), aes(x = lon, y = lat, fill = cor, colour = sig, size = dur, shape = archive)) +
  # define colours & fills
  scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1,1)) +
  scale_colour_manual(values = outlines) +
  scale_shape_manual(values = shapes) +
  scale_size_binned(breaks = c(0, 500, 1000, 1500, 2000)) +
  labs(size = "Record\nduration\n(years)", fill = "Correlation", colour = "Significant\ncorrelation?", shape = "Archive") +
  # tinker with this as desired, to show whichever legends you choose
  guides(fill = "none", colour = "none", size = "none", shape = "none")

temp_map <- world_basemap_rob +
  # add points: shaped by archive type, coloured by correlation, outlined by significance, sized by record duration
  geom_point(data = arrange(temp_rob, abs(cor)), aes(x = lon, y = lat, fill = cor, colour = sig, size = dur, shape = archive)) +
  # define colours & fills
  scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1,1)) +
  scale_colour_manual(values = outlines) +
  scale_shape_manual(values = shapes) +
  scale_size_binned(breaks = c(0, 500, 1000, 1500, 2000)) +
  labs(size = "Record\nduration\n(years)", fill = "Correlation", colour = "Significant\ncorrelation?", shape = "Archive") +
  # tinker with this as desired, to show whichever legends you choose
  guides(shape = "none", colour = "none") +
  theme(legend.position = "bottom", legend.box = "vertical")

# ---------------------
# show maps together
# ---------------------

precIso_map + effMoist_map + temp_map + plot_layout(ncol = 1)
