
## ---------------------------------------------------------------------------
##
## Script name: Figure 4a & 4d
##
## Purpose of script: Recreate Figure 4a & 4d from Konecky et al (2023) Nature Geoscience
##
## Script authors: Georgina Falster & Jonathan Tyler
##
## Date create/updated: 2023-08-23
##
## Email: georgina.falster@anu.edu.au / georgina.falster@gmail.com

## ----------------------------------------------------------------------------
## 1. Perform PCA on Iso2k 'effective moisture' records (for years 1850-2004)
## 2. Compare with dSLP index of Pacific Walker Circulation variability
## 3. Calculate correlation of Iso2k 'effective moisture' PC1 with HadSLP
## 4. Plot correlation of Iso2k 'effective moisture' PC1 with HadSLP, along with PC1 loadings
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
library(ncdf4)
library(vegan)
library(lubridate)
library(smoothr)
library(tidyverse)
library(patchwork)

# =============================================================================
# Load required functions
# =============================================================================

# -------------------------------------------------------------
# you will need the folder 'functions' and then edit the filepaths here so R can read them in
# -------------------------------------------------------------

source("utilities/dineof.r")
source("utilities/iso2kPCA_fun.r")
source("utilities/project_field.R")

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

effMoistTS <- primary.isoTS[which(interpGroup == "EffectiveMoisture")]

# =============================================================================
# Create summary table
# =============================================================================

data_summary <- matrix(NA, length(effMoistTS), 12)
colnames(data_summary) <- c("firstdate", "lastdate", "recordlength", "resolution", "recordName",
                            "lat", "lon", "archive", "name", "infMat", "interp", "direction")

for(i in 1:length(effMoistTS)){
  thisIsoVec <- as.numeric(effMoistTS[[i]]$paleoData_values)
  thisYearVec <- as.numeric(effMoistTS[[i]]$year)
  if(length(thisIsoVec) == length(thisYearVec))
  {theseVec <- na.omit(cbind(thisIsoVec, thisYearVec))}
  iso <- theseVec[,1]
  year <- theseVec[,2]
  data_summary[i,1] <- min(year)
  data_summary[i,2] <- max(year)
  data_summary[i,3] <- max(year) - min(year)
  data_summary[i,4] <- (max(year) - min(year)) / length(iso)
  data_summary[i,5] <- effMoistTS[[i]]$paleoData_TSid
  data_summary[i,6] <- effMoistTS[[i]]$geo_latitude
  data_summary[i,7] <- effMoistTS[[i]]$geo_longitude
  data_summary[i,8] <- effMoistTS[[i]]$archiveType
  data_summary[i,9] <- effMoistTS[[i]]$dataSetName
  data_summary[i,10] <- effMoistTS[[i]]$paleoData_inferredMaterial
  data_summary[i,11] <- effMoistTS[[i]]$isotopeInterpretation1_variable
  data_summary[i,12] <- if
  (is.null(effMoistTS[[i]]$isotopeInterpretation1_direction)) {
    NA
  } else{
    effMoistTS[[i]]$isotopeInterpretation1_direction
  }
}
data_summary <- as.data.frame(data_summary) %>%
  mutate_at(c("lat", "lon", "recordlength"), as.numeric)

# =============================================================================
# Bin data to fixed range
# =============================================================================

# there is legacy code here which isn't actually relevant any more

XMIN = 1850
XMAX = 2004
RES = 3
YEAR <- seq(XMIN, XMAX, RES) #bin edges
binYear <- seq((XMIN + RES/2), (XMAX - RES/2), RES) #bin middles

isoInterpolate <- matrix(NA, length(YEAR), length(effMoistTS))
isoBin <- matrix(NA, length(YEAR)-1, length(effMoistTS))
isoOriginal <- vector("list")
isoOriginalDate <- vector("list")
isoIntName <- matrix(NA, 1, length(effMoistTS))
isoIntTSid <- matrix(NA, 1, length(effMoistTS))


for(i in 1:length(effMoistTS)){
  thisIsoVec <- as.numeric(effMoistTS[[i]]$paleoData_values)
  thisYearVec <- as.numeric(effMoistTS[[i]]$year)
  thisRecordName <- effMoistTS[[i]]$paleoData_TSid
  TSid <- effMoistTS[[i]]$paleoData_TSid

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
      # isoIntTSid[i] <- TSid
    }
  }
}

## screen datasets for just those with >80% data coverage

# First assign isoInt to be the binned data

isoInt <- isoBin
YEAR <- binYear

isoInt.keepCol <- matrix(NA, length(effMoistTS), 1)
isoOriginal.keep <- vector("list")
isoOriginalDate.keep <- vector("list")

for(i in 1:ncol(isoInt)){
  if( (length(na.omit(isoInt[,i])) / length(isoInt[,i])) >= .80
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

# =============================================================================
# Find site means if necessary
# =============================================================================

# -------------------------------------------------------------
# this is actually not necessary here, but would be useful for sites with multiple records. here we just go through the motions...
# -------------------------------------------------------------

uniqueNames <- unique(isoIntName.keep)
isoSiteMeans <- matrix(NA, nrow(isoInt.keep), length(uniqueNames))

for(i in 1:length(uniqueNames)){
  thisSite <- which(isoIntName.keep == uniqueNames[i])
  thisSiteDat <- data.frame(isoInt.keep[,thisSite])
  thisSiteMean <- apply(thisSiteDat, 1, mean)
  isoSiteMeans[,i] <- thisSiteMean
}

# =============================================================================
# Perform Principal Component Analysis
# =============================================================================

iso2kPCA(data=isoSiteMeans, recordName=uniqueNames, data_summary=data_summary)

# -------------------------------------------------------------
# at this point, if you like, you can e.g. inspect the broken stick plots. For example:
# -------------------------------------------------------------

screeplot(pca_out, bstick=TRUE, type="lines", main="")

# -------------------------------------------------------------
# tidy the outputs
# -------------------------------------------------------------

effMoist_load <- PCA_results %>%
  mutate_at(c("PC1", "PC2", "PC3", "PC4"), as.numeric) %>%
  mutate(PC1_adj = PC1) %>% # FLIP IF NECESSARY HERE
  dplyr::rename(uID = recordName)

effMoist_ts <- PCA_timeseries %>%
  as.data.frame() %>%
  mutate(PC1_adj = PC1) %>% # FLIP IF NECESSARY HERE (make sure you do the same as in the above)
  mutate(binmid = binYear)

# =============================================================================
# Read in and arrange the dSLP index for Pacific Walker Circulation variability
# =============================================================================

pwc_index <- read.csv("data/pwc_raw.csv", header = TRUE) %>%
  rename(gradient = Value, year = yy, month = mm) %>%
  filter(year >= 1850 & year <= 2004)

month_means <- pwc_index %>%
  group_by(month) %>%
  summarise(month_mean = mean(gradient)) %>%
  ungroup

pwc_index$pwc_monthly <- NA

## ----------------------------------------------------------------------------
# calculate index (deviation from long-term climatology)
## ----------------------------------------------------------------------------

for(i in 1:(nrow(pwc_index))) {

  pwc_index$pwc_monthly[i] <- pwc_index$gradient[i] - month_means$month_mean[which(month_means$month == pwc_index$month[i])]
}

## ----------------------------------------------------------------------------
# calculate annual index
## ----------------------------------------------------------------------------

pwc_annual <- pwc_index %>%
  group_by(year) %>%
  summarise(pwc_annual = mean(pwc_monthly)) %>%
  ungroup()

## ----------------------------------------------------------------------------
# bin the index, using same parameters as when preparing the Iso2k records above
## ----------------------------------------------------------------------------

XMIN = 1850
XMAX = 2004
RES = 3
YEAR <- seq(XMIN, XMAX, RES) #bin edges
binYear <- seq((XMIN + RES/2), (XMAX - RES/2), RES) #bin middles

pwc_binned <- as.vector(tapply(pwc_annual$pwc_annual, cut(pwc_annual$year, YEAR), mean))

# =============================================================================
# Show dSLP and Iso2k 'effective moisture' PC1 on the same plot (Fig. 4d)
# =============================================================================

# -------------------------------------------------------------
# some preparation for easier plotting
# -------------------------------------------------------------

pca_full_with_pwc <- data.frame(year = effMoist_ts$binmid, PWC_scaled = scale(pwc_binned), EffectiveMoisture_PC1 = effMoist_ts$PC1_adj)

pca_dSLP_long <- pivot_longer(pca_full_with_pwc, cols = -year, names_to = "ts", values_to = "val")

# -------------------------------------------------------------
# make the timeseries plot
# -------------------------------------------------------------

effmoist_pc1_with_dslp <- ggplot(data = pca_dSLP_long) +
  # add the Iso2k 'effective moisture' PC1 and dSLP timeseries
  geom_line(aes(x = year, y = val, group = ts, colour = ts), linewidth = 1) +
  scale_colour_manual(values = c("PWC_scaled" = "grey50", "EffectiveMoisture_PC1" = "black")) +
  scale_x_continuous(breaks = seq(from = 1850, to = 2005, by = 10),
                     minor_breaks = seq(from = 1850, to = 2005, by = 5),
                     expand = expansion(add = c(2,0))) +
  xlab("Year CE") + ylab("Scaled value") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        text = element_text(size = 14),
        axis.text.x = element_text(size = 12))

effmoist_pc1_with_dslp

# =============================================================================
# Read in gridded SLP data from HadSLP
# =============================================================================

url_grid <- "ftp://ftp.cdc.noaa.gov/Datasets.other/hadslp2/slp.mnmean.real.nc"


download.file(url_grid, "data/hadslp2.slp.mnmean.nc", method = "auto",
              quiet = FALSE, mode="wb", cacheOK = TRUE)

hadslp <- nc_open("data/hadslp2.slp.mnmean.nc")

## ----------------------------------------------------------------------------
# get coordinates and time index
## ----------------------------------------------------------------------------

lon <- ncvar_get(hadslp, "lon")

lat <- ncvar_get(hadslp, "lat")

time <- ncvar_get(hadslp, "time")
time_units <- ncatt_get(hadslp, "time", "units")

## ----------------------------------------------------------------------------
# extract SLP variable and close file
## ----------------------------------------------------------------------------

slp_array <- ncvar_get(hadslp, "slp")
long_name <- ncatt_get(hadslp, "slp", "long_name")
mslp_units <- ncatt_get(hadslp, "slp", "units")
fill_value <- ncatt_get(hadslp, "slp", "missing_value")

nc_close(hadslp)

## ----------------------------------------------------------------------------
# tidy up extracted components
## ----------------------------------------------------------------------------

# convert fill values to NA

slp_array[slp_array == fill_value$value] <- NA

# convert the time index into something useable; it's packaged as 'days since a certain time'

unit_split <- strsplit(time_units$value, " ")
origin <- ymd(unit_split[[1]][3])

date_index <- as.Date(origin + as.difftime(time, units = "days"))
year_index <- year(date_index)

rm(unit_split)

# =============================================================================
# Bin the SLP data to match the PCA, and calculate correlations
# =============================================================================

ind_cells <- expand.grid(lon, lat)
colnames(ind_cells) <- c("lon", "lat")

ind_cells$cor_em_pc1 <- NA
ind_cells$cor_pwc <- NA

for(i in 1:(nrow(ind_cells))) {

  lon_ind <- which(lon %in% ind_cells$lon[i])
  lat_ind <- which(lat %in% ind_cells$lat[i])

  this_cell_vec <- slp_array[lon_ind,lat_ind, ]

  if(all(is.na(this_cell_vec))) {
    next
  } else{

    thisRec <- data.frame(val = this_cell_vec, date = date_index) %>%
      mutate(year = year(date)) %>%
      group_by(year) %>%
      summarise(mean_slp = mean(val)) %>%
      ungroup()

    thisRec_binned <- as.vector(tapply(thisRec$mean_slp, cut(thisRec$year, YEAR), mean)) #bin using same edges as in the PCA

    ind_cells$cor_em_pc1[i] <- cor(thisRec_binned, effMoist_ts$PC1_adj) # 'PC1_adj' is the variable where PC1 has been flipped (if neccessary) to positively correlate with dSLP
    ind_cells$cor_pwc[i] <- cor(thisRec_binned, pwc_binned)
  }
}

rm(thisRec, lon_ind, lat_ind, this_cell_vec, thisRec_binned)

# =============================================================================
# Some initial mapping preparations
# =============================================================================

# -------------------------------------------------------------
# at this point, head over to the 'create_world_basemap_PCs.R' script and run the whole thing
# -------------------------------------------------------------

proj_string_data <- "+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

trans <- effMoist_load
coordinates(trans) <- c("lon", "lat")
proj4string(trans) <- CRS("+proj=longlat +datum=WGS84")
trans <- spTransform(trans, CRS(proj_string_data))
effMoist_rob <- as.data.frame(trans)
rm(trans)

# -------------------------------------------------------------
# Now project the gridded SLP data. Note that this is *absolutely not* the best way to do this(!) but it is what I did for this paper
# -------------------------------------------------------------

project_field(clim_dat = ind_cells, cell_size = 5, proj_new = proj_string_data)
slp_proj <- proj_main

rm(proj_main)

# =============================================================================
# Now make the SLP~dSLP correlation contours
# =============================================================================

WGS84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

## ----------------------------------------------------------------------------
# convert ind_cells to spatial points
## ----------------------------------------------------------------------------

pwc_coords <- cbind(ind_cells$lon, ind_cells$lat)

sp_pwc <- SpatialPointsDataFrame(coords = pwc_coords, data = data.frame(ind_cells$cor_pwc))

## ----------------------------------------------------------------------------
# make an empty raster and fill it in
## ----------------------------------------------------------------------------

cell_size <- 5

ncols <- ((max(ind_cells$lon) - min(ind_cells$lon))/cell_size) + 1; nrows <- ((max(ind_cells$lat) - min(ind_cells$lat))/cell_size) + 1

grid_pwc <- raster::raster(nrows = nrows, ncols = ncols,
                   xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat),
                   res = cell_size, crs = "+proj=longlat +datum=WGS84")

grid_pwc <- raster::rasterize(sp_pwc, grid_pwc, field = "ind_cells.cor_pwc", fun = mean)

smoothed_grid_pwc <- raster::disaggregate(grid_pwc, 10, method = 'bilinear')

## ----------------------------------------------------------------------------
# now following steps from here: https://stackoverflow.com/questions/14379828/how-does-one-turn-contour-lines-into-filled-contours
## ----------------------------------------------------------------------------

rc <- raster::cut(smoothed_grid_pwc, breaks= seq(-1, 1, 0.2))

rc_pos <- raster::cut(smoothed_grid_pwc, breaks= seq(0, 1, 0.2))
rc_neg <- raster::cut(smoothed_grid_pwc, breaks= seq(-1, 0, 0.2))

pols <- raster::rasterToPolygons(rc, dissolve=TRUE)
pols_pos <- raster::rasterToPolygons(rc_pos, dissolve=TRUE)
pols_neg <- raster::rasterToPolygons(rc_neg, dissolve=TRUE)

smoothed_pols <- smoothr::smooth(pols, method = "chaikin")

pols_pos <- smoothr::smooth(pols_pos, method = "chaikin")
pols_neg <- smoothr::smooth(pols_neg, method = "chaikin")

# -------------------------------------------------------------
# now project the positive contours
# -------------------------------------------------------------

shift <- 180

split.line = SpatialLines(list(Lines(list(Line(cbind(180-shift,c(-90,90)))), ID="line")),
                          proj4string=WGS84)

# intersecting line with contour polygons
line.gInt <- gIntersection(split.line, pols_pos)

# create a very thin polygon (buffer) out of the intersecting "split line"
bf <- gBuffer(line.gInt, byid=TRUE, width=0.000001)

# split contour polygons using intersecting thin polygon (buffer)
contour.split <- gDifference(pols_pos, bf, byid=TRUE)

# give the PROJ.4 string for the Robinson projection
rob_proj <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# transform split contour polygons in a data table that ggplot can use
Contour.DT.pos <- data.table(map_data(as(contour.split, "SpatialPolygonsDataFrame")))
# Shift coordinates
Contour.DT.pos[, long.new := long + shift]
Contour.DT.pos[, long.new := ifelse(long.new > 180, long.new-360, long.new)]
# project coordinates
Contour.DT.pos[, c("X","Y") := data.table(project(cbind(long.new, lat), proj=rob_proj))]

# -------------------------------------------------------------
# repeat for negative contours
# -------------------------------------------------------------

split.line = SpatialLines(list(Lines(list(Line(cbind(180-shift,c(-90,90)))), ID="line")),
                          proj4string=WGS84)

# intersecting line with contour polygons
line.gInt <- gIntersection(split.line, pols_neg)

# create a very thin polygon (buffer) out of the intersecting "split line"
bf <- gBuffer(line.gInt, byid=TRUE, width=0.000001)

# split contour polygons using intersecting thin polygon (buffer)
contour.split <- gDifference(pols_neg, bf, byid=TRUE)


# give the PROJ.4 string for the Robinson projection (if you want Eckert IV projection, use "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
rob_proj <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# transform split contour polygons in a data table that ggplot can use
Contour.DT.neg <- data.table(map_data(as(contour.split, "SpatialPolygonsDataFrame")))
# Shift coordinates
Contour.DT.neg[, long.new := long + shift]
Contour.DT.neg[, long.new := ifelse(long.new > 180, long.new-360, long.new)]
# project coordinates
Contour.DT.neg[, c("X","Y") := data.table(project(cbind(long.new, lat), proj=rob_proj))]

rm(WGS84, pwc_coords, sp_pwc, cell_size, ncols, grid_pwc, rc, pols, pols_pos, pols_neg, shift, split.line, line.gInt, bf, contour.split, rob_proj)

# =============================================================================
# Define vector for colourbar (matches the MATLAB red-blue gradient)
# =============================================================================

matlab_cols <- c("#0000FF", "#1A1AFF", "#3333FF", "#4D4DFF", "#6666FF", "#8080FF", "#9999FF", "#B3B3FF", "#CDCDFF", "#E6E6FF", "#FFFFFF",
                 "#FFFFFF", "#FFE6E6", "#FFCDCD", "#FFB3B3", "#FF9999", "#FF8080", "#FF6666", "#FF4D4D", "#FF3333", "#FF1A1A", "#FF0000")

matlab_ramp <- colorRampPalette(matlab_cols)

# =============================================================================
# Make map showing correlation of Iso2k 'effective moisture' PC1 with SLP (Fig. 4a)
# =============================================================================

effMoist_SLP_map <- world_basemap_rob +
  # SLP correlation grid
  geom_rect(data = na.omit(slp_proj), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cor_em_pc1), colour = "#00000000") +
  # overlay continent outlines
  geom_polygon(data = Country.DT,aes(x = X, y = Y, group = group),colour = "gray30",fill = "#00000000",size = 0.5) +
  # overlay PWC-SLP correlation contours
  geom_polygon(data = Contour.DT.pos,aes(x = X, y = Y, group = group),colour = "black",fill = "#00000000",size = 0.5) +
  geom_polygon(data = Contour.DT.neg,aes(x = X, y = Y, group = group),colour = "black",fill = "#00000000",size = 0.5, linetype = "69") + # 6 on 9 off
  # bounding box
  geom_path(data = grid.DT[(long %in% c(-180,180) & region == "NS")
                           |(long %in% c(-180,180) & lat %in% c(-90,90) & region == "EW")],
            aes(x = X, y = Y, group = group),
            linetype = "solid", colour = "black", size = .7) +
  # show Iso2k 'effective moisture' PC1 loadings (EOF1)
  geom_point(data = effMoist_rob, aes(lon, lat, fill = PC1_adj), shape = 21, size = 4) +
  scale_fill_stepsn(colours = matlab_ramp(21), breaks = seq(-1, 1, 0.1), limits = c(-1, 1)) +
  ggtitle("1850-2005 (Iso2k)") +
  labs(fill = NULL) +
  theme(
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)) +
  guides(fill= guide_colorbar(barwidth=unit(0.4, "npc"))) +
  theme(legend.position = "bottom")

effMoist_SLP_map
