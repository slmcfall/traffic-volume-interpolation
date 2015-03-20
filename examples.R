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
meuse.variogram <- createVariogram(c("zinc~1"), meuse, meuse.width, meuse.cutoff)
meuse.variogramaf <- autofitVariogram(zinc~1, meuse)

# BUFFER STUFF

# use range of each adjacent state to buffer State A

adjWA <- getAdjStates("WA")

range_df <- getStateRanges(adjWA, projection, c("AADT"), c("AADT~1"))

# do.call(function(ranges,...) print(ranges), range_df )

# iterate over dimensions of range_df
# dim(range_df)[1]

singleStateAbbrv <- selectState(as.character(range_df$adjStates[2]))
singleStateRange <- range_df$ranges[2]

WA <- selectState("WA")
adjState_buff <- gBuffer(WA, width = singleStateRange) #, byid = TRUE)

# clip buffered State A using adjacent state

clipped_range_buffer <- gIntersection(adjState_buff, singleStateAbbrv )

# clip points outside of State A according to clipped buffers

ID_points <- createSpdf("/home/sean/traffic/WA_example/ID", c("AADT"), projection)
ID_points_clipped <- gIntersection(clipped_range_buffer, ID_points)
WA_points <- createSpdf("/home/sean/traffic/WA_example/WA", c("AADT"), projection)

# merge clipped points 

WA_points_all <- spRbind(WA_points, ID_points_clipped)

