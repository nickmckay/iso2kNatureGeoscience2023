## ---------------------------------------------------------------------------
##
## Script name: Fig. 2 (basemap)
##
## Purpose of script: create a global basemap, on a Pacific-centred Robinson projection
##
## Script author: Valentin Ștefan, ed. Georgina Maja Falster
##
## Date create/updated: 2023-02-23
##
## Email: georgina.falster@anu.edu.au (institutional) / georgina.falster@gmail.com (permanent)
##
## Citation: Falster et al. (2023) Nature 

## ----------------------------------------------------------------------------
##
## This script adapted from:
## https://github.com/valentinitnelav/valentinitnelav.github.io/blob/master/gallery/Pacific%20centered%20world%20map%20with%20ggplot.R#L137
## 
## ----------------------------------------------------------------------------
# =============================================================================
# Load coastline polygons
# =============================================================================

coastline <- readOGR("spatialData/GSHHS_c_L1.shp")
antarctica <- readOGR("spatialData/GSHHS_l_L6.shp")

wmap <- rbind(coastline, antarctica)
class(wmap) # is a SpatialPolygonsDataFrame object

# =============================================================================
# Split world map by "split line"
# =============================================================================

# inspired from:
# https://stat.ethz.ch/pipermail/r-sig-geo/2015-July/023168.html

# shift central/prime meridian towards west - positive values only
shift <- 180 #OG 180 + 30

# create "split line" to split country polygons
WGS84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
split.line = SpatialLines(list(Lines(list(Line(cbind(180-shift,c(-90,90)))), ID="line")), 
                          proj4string=WGS84)

# NOTE - in case of TopologyException' errors when intersecting line with country polygons,
# apply the gBuffer solution suggested at:
# http://gis.stackexchange.com/questions/163445/r-solution-for-topologyexception-input-geom-1-is-invalid-self-intersection-er
# wmap <- gBuffer(wmap, byid=TRUE, width=0)

# intersecting line with country polygons
line.gInt <- gIntersection(split.line, wmap)

# create a very thin polygon (buffer) out of the intersecting "split line"
bf <- gBuffer(line.gInt, byid=TRUE, width=0.000001)  

# split country polygons using intersecting thin polygon (buffer)
wmap.split <- gDifference(wmap, bf, byid=TRUE)
# plot(wmap.split) # check map
class(wmap.split) # is a SpatialPolygons object

# =============================================================================
# Create graticules
# =============================================================================

# create a bounding box - world extent
b.box <- as(raster::extent(-180, 180, -90, 90), "SpatialPolygons")

# assign CRS to box
proj4string(b.box) <- WGS84

# create graticules/grid lines from box
grid <- gridlines(b.box, 
                  easts  = seq(from=-180, to=180, by=10),
                  norths = seq(from=-90, to=90, by=10))

# create labels for graticules
grid.lbl <- labels(grid, side = 1:2)

# transform labels from SpatialPointsDataFrame to a data table that ggplot can use
grid.lbl.DT <- data.table(grid.lbl@coords, grid.lbl@data)

# filter so we only have labels every 30 degrees (whilst still having graticules at 10 degrees)
#grid.lbl.DT <- grid.lbl.DT[coords.x1 %in% seq(from=-180, to=180, by=30)] #for longitude labels every 30 degrees
grid.lbl.DT <- grid.lbl.DT[coords.x1 %in% seq(from=-180, to=180, by=60)] #for longitude labels every 60 degrees
grid.lbl.DT <- grid.lbl.DT[coords.x2 %in% seq(from=-90, to=90, by=30)]

# prepare labels with regular expression:
# - delete any unwanted labels i.e. this at the top or bottom or sides
grid.lbl.DT[, labels := gsub(pattern="180\\*degree|90\\*degree\\*N|90\\*degree\\*S", replacement="", x=labels)]

grid.lbl.DT$labels[which(grid.lbl.DT$coords.x1 == -180 & grid.lbl.DT$coords.x2 == -90 & grid.lbl.DT$pos == 1)] <- "180*degree"
grid.lbl.DT$labels[which(grid.lbl.DT$coords.x1 == 180 & grid.lbl.DT$coords.x2 == -90 & grid.lbl.DT$pos == 1)] <- "0*degree"

# - replace pattern "*degree" with "°" (* needs to be escaped with \\)
grid.lbl.DT[, lbl := gsub(pattern="\\*degree", replacement="°", x=labels)]

# - delete any remaining "*"
grid.lbl.DT[, lbl := gsub(pattern="*\\*", replacement="", x=lbl)]

# adjust coordinates of labels so that they fit inside the globe **EDIT I don't want this - I actually want them outside the globe!
# grid.lbl.DT[, long := ifelse(coords.x1 %in% c(-180,180), coords.x1*175/180, coords.x1)]
# grid.lbl.DT[, lat  := ifelse(coords.x2 %in% c(-90,90), coords.x2*82/90, coords.x2)]
grid.lbl.DT[, long := ifelse(coords.x1 %in% c(-180,180), coords.x1, coords.x1)]
grid.lbl.DT[, lat  := ifelse(coords.x2 %in% c(-90,90), coords.x2, coords.x2)]

# =============================================================================
# Prepare data for ggplot, shift & project coordinates
# =============================================================================

# give the PROJ.4 string for the Robinson projection (if you want Eckert IV projection, use "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 

# transform graticules from SpatialLines to a data table that ggplot can use
grid.DT <- data.table(map_data(SpatialLinesDataFrame(sl=grid, 
                                                     data=data.frame(1:length(grid)), 
                                                     match.ID = FALSE)))
# project coordinates
# assign matrix of projected coordinates as two columns in data table
grid.DT[, c("X","Y") := data.table(project(cbind(long, lat), proj=PROJ))]

# project coordinates of labels
grid.lbl.DT[, long.new := long + shift]
grid.lbl.DT[, long.new := ifelse(long.new > 180, long.new-360, long.new)]
grid.lbl.DT[, long.new := ifelse(long.new == 0 & pos == 2, -180, long.new)]
grid.lbl.DT[, long.new := ifelse(long.new == 0 & lbl == "0°", -180, long.new)]
grid.lbl.DT[, c("X","Y") := data.table(project(cbind(long.new, lat), proj=PROJ))]

# transform split country polygons in a data table that ggplot can use
Country.DT <- data.table(map_data(as(wmap.split, "SpatialPolygonsDataFrame")))
# Shift coordinates
Country.DT[, long.new := long + shift]
Country.DT[, long.new := ifelse(long.new > 180, long.new-360, long.new)]
# project coordinates 
Country.DT[, c("X","Y") := data.table(project(cbind(long.new, lat), proj=PROJ))]

# =============================================================================
# Plot
# =============================================================================

world_basemap_rob <- ggplot() + 
  # add projected countries
  geom_polygon(data = Country.DT, 
               aes(x = X, y = Y, group = group), 
               colour = "gray30", 
               fill = "gray90", 
               size = 0.6) +
  # add graticules
  #geom_path(data = grid.DT, 
  #          aes(x = X, y = Y, group = group), 
  #          linetype = "dotted", colour = "grey50", size = .25) +
  # add a bounding box (select graticules at edges)
  geom_path(data = grid.DT[(long %in% c(-180,180) & region == "NS")
                           |(long %in% c(-180,180) & lat %in% c(-90,90) & region == "EW")], 
            aes(x = X, y = Y, group = group), 
            linetype = "solid", colour = "black", size = .3) +
  # add graticule labels (longitudes)
  geom_text(data = filter(grid.lbl.DT, lat == -90),
            aes(x = X, y = Y-960000, label = lbl, angle = 30), 
            colour = "grey30", size = 6) +
  # add graticule labels (latitudes)
  geom_text(data = filter(grid.lbl.DT, long == -180 & pos == 2),
            aes(x = X-1450000, y = Y, label = lbl), 
            colour = "grey30", size = 6) +
  # ensures that one unit on the x-axis is the same length as one unit on the y-axis
  coord_equal() + # same as coord_fixed(ratio = 1)
  # set empty theme
  theme_void()

rm(b.box, bf, grid, grid.lbl, line.gInt, split.line, WGS84, wmap, wmap.split, shift, PROJ, antarctica, coastline)
