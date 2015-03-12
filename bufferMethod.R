#
# 
#
#
#
#
#
#
#
#

# get kriging functions, and all applicable libraries
source("~/Documents/trafficVolume/functions.R")
projection <- "+proj=aea +lat_1=38 +lat_2=41 +lat_0=34 +lon_0=-114 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs"

# bring in polygon for whole US
cUS <- readOGR(dsn = "/home/sean/traffic/borders/cUS.shp", layer = "cUS")

# select the polygon of interest, say State A

selectState <- function(stateName) {
  state <- cUS[cUS$STATE_ABBR == stateName,]
  return(state)
}

# get all states adjacent to State A

adjTable <- gTouches(cUS, byid = TRUE)

stateFactor <- cUS@data$STATE_ABBR
stateVector <- as.vector(stateFactor)

colnames(adjTable) <- stateVector
rownames(adjTable) <- stateVector

getAdjStates <- function(stateAbbrv) {
  adjStates <- which(adjTable[stateAbbrv,] == "TRUE")
  adjStates <- names(adjStates)
  
  return(adjStates)
}

# calculate range for each adjacent state

getStateRanges <- function(adjStates, proj4str, column_names, equation) {

  proj <- proj4str
  ranges <- c()
  
  for (stateAbbrv in adjStates) {
    state_abbrv <- stateAbbrv
    vct_path <- paste("/home/sean/traffic/WA_example/", state_abbrv, sep = "")
    column_names <- column_names
    
    spdf <- createSpdf(vct_path, column_names, proj)
    raster <- createRaster(spdf = spdf, resolution = 1000)
    avg_spdf <- createAvgSpdf(spdf, raster, column_names)
    # HOW should the width and cutoff be generated?
    variogram <- createVariogram(equation, spdf, 300, 5000)
    
    range_state <- round(getRange(variogram))
    print(range_state)
    ranges <- c(ranges, range_state)
  }
  
  return (ranges)
  
}

# use range of each adjacent state to buffer State A

# clip buffered State A using adjacent state

# merge State A, and clipped buffers

# clip points outside of State A according to clipped buffers

# merge clipped points 