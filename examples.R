# examples.R

source("~/Documents/trafficVolume/functions.R")

# VARIOGRAM 

data(meuse)
coordinates(meuse) <- ~x+y
meuse.proj4string <- "+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs "
proj4string(meuse) <- meuse.proj4string
data(meuse.grid)
gridded(meuse.grid) <- ~x+y
meuse.raster <- raster(meuse.grid)
# my functions
meuse.raster <- createRaster(meuse, 100)
meuse.avgspdf <- createAvgSpdf(meuse, meuse.raster, c("zinc"))
#
meuse.variogram <- createVariogram(c("zinc~1"), meuse, meuse.raster)
meuse.variogramaf <- autofitVariogram(zinc~1, meuse)

createVariogram <- function(equation, spdf, rast) {
  
  #############
  
  # PURPOSE
  # creates variogram through modded autoFitVariogram function
  
  # INPUTS
  # equation = equation for the variogram model (ex: AADT ~ 1)
  # spdf     = SpatialPointsDataFrame
  # width    =  width between pairs of points
  
  # OUTPUTS
  # variogram = variogram for kriging
  
  #############
  
  opts <- list(orig.behavior = FALSE)
  
  # generate the cutoff or range
  var_cutoff <- getCutoff(spdf)
  
  # generate the width or lag size
  var_width <- (res(rast)[1]) / 2
  
  variogram <- afvmod(as.formula(equation), input_data = spdf, width = var_width, cutoff = var_cutoff, 
                      verbose = TRUE) # , miscFitOptions = opts)
  
  return (variogram)
}

getCutoff <- function(spdf) {
  spdf_bbox <- bbox(spdf)
  x1 <- spdf_bbox[1]
  x2 <- spdf_bbox[2]
  y1 <- spdf_bbox[3]
  y2 <- spdf_bbox[4]
  
  spdf_bbox_dist <- sqrt((x1-x2)^2 + (y1-y2)^2)
  spdf_range <- spdf_bbox_dist * .35
  
  return(spdf_range) 
}

getWidth <- function(rast) {
  width <- res(rast)[1] / 2
  return (width)
}