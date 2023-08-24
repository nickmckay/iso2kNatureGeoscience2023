
project_field <- function(clim_dat, cell_size, proj_orig = "+proj=longlat +datum=WGS84", proj_new) {
  
  shift_from_centre <- cell_size/2
  
  max_lat <- 90-shift_from_centre
  min_lat = -90+shift_from_centre
  
  clim_dat$xmax <- clim_dat$lon + shift_from_centre
  clim_dat$xmin <- clim_dat$lon - shift_from_centre
  clim_dat$ymax <- clim_dat$lat + shift_from_centre
  clim_dat$ymin <- clim_dat$lat - shift_from_centre

  clim_dat$lon_centre <- clim_dat$lon
  clim_dat$lat_centre <- clim_dat$lat
  
  clim_dat$lon_centre1 <- clim_dat$lon
  clim_dat$lat_centre1 <- clim_dat$lat
  
  clim_dat$lon_centre2 <- clim_dat$lon
  clim_dat$lat_centre2 <- clim_dat$lat

  
  clim_dat <- filter(clim_dat, clim_dat$lat <= max_lat & clim_dat$lat >= min_lat)
  
  # we only need the sliver if 360 longitude is not on an edge! here it is, so we're fine. 
  
  init_main <- na.omit(filter(clim_dat, clim_dat$lon < 360-shift_from_centre & clim_dat$lon > 0+shift_from_centre))
  #init_sliver <- na.omit(filter(clim_dat, clim_dat$lon > 360))
  
  trans <- init_main
  coordinates(trans) <- c("xmax", "lat_centre1")
  proj4string(trans) <- CRS(proj_orig)
  trans <- spTransform(trans, CRS(proj_new))
  temp <- as.data.frame(trans)
  
  trans <- temp
  coordinates(trans) <- c("xmin", "lat_centre2")
  proj4string(trans) <- CRS(proj_orig)
  trans <- spTransform(trans, CRS(proj_new))
  temp <- as.data.frame(trans)
  
  trans <- temp
  coordinates(trans) <- c("lon_centre1", "ymax")
  proj4string(trans) <- CRS(proj_orig)
  trans <- spTransform(trans, CRS(proj_new))
  temp <- as.data.frame(trans)
  
  trans <- temp
  coordinates(trans) <- c("lon_centre2", "ymin")
  proj4string(trans) <- CRS(proj_orig)
  trans <- spTransform(trans, CRS(proj_new))
  temp <- as.data.frame(trans)
  
  trans <- temp
  coordinates(trans) <- c("lon_centre", "lat_centre")
  proj4string(trans) <- CRS(proj_orig)
  trans <- spTransform(trans, CRS(proj_new))
  proj_main <- as.data.frame(trans)
  
  proj_main$lat_centre1 <- NULL
  proj_main$lat_centre2 <- NULL
  proj_main$lon_centre1 <- NULL
  proj_main$lon_centre2 <- NULL
  
  # trans <- init_sliver
  # coordinates(trans) <- c("top_left_x", "top_left_y")
  # proj4string(trans) <- CRS(proj_orig)
  # trans <- spTransform(trans, CRS(proj_new))
  # temp <- as.data.frame(trans)
  # 
  # trans <- temp
  # coordinates(trans) <- c("top_right_x", "top_right_y")
  # proj4string(trans) <- CRS(proj_orig)
  # trans <- spTransform(trans, CRS(proj_new))
  # temp <- as.data.frame(trans)
  # 
  # trans <- temp
  # coordinates(trans) <- c("bot_left_x", "bot_left_y")
  # proj4string(trans) <- CRS(proj_orig)
  # trans <- spTransform(trans, CRS(proj_new))
  # temp <- as.data.frame(trans)
  # 
  # trans <- temp
  # coordinates(trans) <- c("bot_right_x", "bot_right_y")
  # proj4string(trans) <- CRS(proj_orig)
  # trans <- spTransform(trans, CRS(proj_new))
  # proj_sliver <- as.data.frame(trans)
  
  proj_main <<- proj_main
  # proj_sliver <<- proj_sliver
  
}

