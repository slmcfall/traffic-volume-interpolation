# main.R

# define directory with data
data.directory <- "/home/sean/traffic/AADT/"
data.suffix <- "_AADT"

# call functions and libraries
source("~/Documents/trafficVolume/functions.R")
source("~/Documents/trafficVolume/bufferMethod.R")

# define projection, albers equal area
projection <- "+proj=aea +lat_1=38 +lat_2=41 +lat_0=34 +lon_0=-114 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs"

# polygon for the continental US
cUS <- readOGR(dsn = "/home/sean/Documents/trafficVolume/cUS.shp", layer = "cUS")

# state to be buffered
stateAbbreviation <- "WA"


#####

adjacent.states <- getAdjStates(state.abbrv = stateAbbreviation)

state.ranges <- getStateRanges(adj.states = adjacent.states, 
                               proj.str   = projection, 
                               col.names  = c("AADT"), 
                               equation   = c("AADT ~ 1"), 
                               pnt.suffix = data.suffix, 
                               pnt.dir    = data.directory, 
                               rst.res    = 1000
                               )

buffered.spdf <- getBufferedSpdf(state.abbrv = stateAbbreviation,
                                 pnt.dir     = data.directory,
                                 col.names   = c("AADT"),
                                 proj.str    = projection,
                                 ranges.df   = state.ranges,
                                 pnt.suffix  = data.suffix,
                                 adj.states  = adjacent.states)

WA.raster <- createRaster(spdf = buffered.spdf, resolution = 1000 )
#WA.avg.spdf <- createAvgSpdf(spdf = buffered.spdf, rast = WA.raster, col_names = c("AADT"))
WA.cutoff <- getBufferCutoff(spdf = buffered.spdf)
WA.width <- getBufferWidth(rast = WA.raster)
WA.variogram <- createVariogram(equation = c("AADT~1"), spdf = buffered.spdf, varWidth = WA.width, varCutoff = WA.cutoff)

# the data I automatically created with just AADT, AADT column not showing up for some reason....
# the problem is buffered.spdf is creating a SpatialPoints class, not a SpatialPointsDataFrame