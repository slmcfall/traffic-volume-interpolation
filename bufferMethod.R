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

  ranges <- c()
  
  for (stateAbbrv in adjStates) {
    state_abbrv <- stateAbbrv
    vct_path <- paste("/home/sean/traffic/WA_example/", state_abbrv, sep = "")
    column_names <- column_names
    
    spdf <- createSpdf(vct_path, column_names, proj4str)
    raster <- createRaster(spdf = spdf, resolution = 1000)
    avg_spdf <- createAvgSpdf(spdf, raster, column_names)
    # HOW should the width and cutoff be generated?
    buffer_width <- getBufferWidth(raster)
    buffer_cutoff <- getBufferCutoff(spdf)
    
    variogram <- createVariogram(equation, spdf, buffer_width, buffer_cutoff)
    
    range_state <- round(getRange(variogram))
    ranges <- c(ranges, range_state)
  }
  
  ranges_df <- data.frame(adjStates, ranges)
  return (ranges_df)
}

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



getBufferedSpdf <- function(state.abbrv, directory, column.names, projection, ranges.df) {
  
  originalState_border <- selectState(state.abbrv)
  
  vct_path <- paste(directory, state.abbrv, sep = "")
  originalState_points <- createSpdf(vct_path, column.names, projection)

  for (state in seq(1:dim(ranges.df)[1])){
    
    stateAbbrv <- as.character(ranges.df$adjState[state])
    
    # select border of adjacent state
    singleState_border <- selectState(stateAbbrv)
    # select range of adjacent state
    singleState_range <- (ranges.df$ranges[state])
    
    # create buffer of original state with adjacent state range
    originalState_buffer <- gBuffer(originalState, width = singleState_range)
    
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

























