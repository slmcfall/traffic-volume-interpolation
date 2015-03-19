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
meuse.cutoff <- getBufferCutoff(meuse)
meuse.width <- getBufferWidth(meuse.raster)
meuse.variogram <- createVariogram(c("zinc~1"), meuse, 3000, 5000)
meuse.variogramaf <- autofitVariogram(zinc~1, meuse)

