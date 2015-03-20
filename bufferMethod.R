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

getAdjStates <- function(state.abbrv) {
  adjStates <- which(adjTable[state.abbrv,] == "TRUE")
  adjStates <- names(adjStates)
  
  return(adjStates)
}

# calculate range for each adjacent state

getStateRanges <- function(adj.states, proj.str, col.names, equation, pnt.suffix, pnt.dir, rst.res) {

  ranges <- c()
  
  for (state.abbrv in adj.states) {
    vct.path <- paste(pnt.dir, state.abbrv, pnt.suffix, sep = "")
    
    spdf <- createSpdf(vct.path, col.names, proj.str)
    raster <- createRaster(spdf = spdf, resolution = rst.res)
    avg_spdf <- createAvgSpdf(spdf, raster, col.names)
    
    buffer.width <- getBufferWidth(raster)
    buffer.cutoff <- getBufferCutoff(spdf)
    
    variogram <- createVariogram(equation, spdf, buffer.width, buffer.cutoff)
    
    range.state <- round(getRange(variogram))
    ranges <- c(ranges, range.state)
  }
  
  ranges.df <- data.frame(adj.states, ranges)
  return (ranges.df)
}

getBufferedSpdf <- function(state.abbrv, pnt.dir, col.names, proj.str, ranges.df, pnt.suffix, adj.states) {
  
  originalState_border <- selectState(state.abbrv)
  
  vct.path <- paste(pnt.dir, state.abbrv, pnt.suffix, sep = "")
  originalState_points <- createSpdf(vct.path, col.names, projection)

  for (state in seq(1:dim(ranges.df)[1])){
    
    stateAbbrv <- as.character(ranges.df$adj.states[state])
    
    # select border of adjacent state
    singleState_border <- selectState(stateAbbrv)
    # select range of adjacent state
    singleState_range <- (ranges.df$ranges[state])
    
    # create buffer of original state with adjacent state range
    originalState_buffer <- gBuffer(originalState_border, width = singleState_range)
    
    # pull out points for the adjacent state
    vct_path <- paste("/home/sean/traffic/WA_example/", stateAbbrv, sep = "")
    column_names <- c("AADT")
    
    adjState_points <- createSpdf(vct_path, column_names, projection)
    
    # clip out the points that intersect the buffer
    adjState_clipped_points <- gIntersection(originalState_buffer, adjState_points)
    
    # merge the points to the original state
    originalState_points <- spRbind(originalState_points, adjState_clipped_points)
  } 
  return(originalState_points)
}

























