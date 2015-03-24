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
  
  # data frame to be combined with original state
  add.df <- data.frame(x=numeric(0), y=numeric(0), AADT=numeric(0))

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
    adjState_clipped_points <- gIntersects(originalState_buffer, adjState_points, byid = TRUE)
    
    # change the column of the intersection data frame to something constant
    colnames(adjState_clipped_points) <- c("Boolean")
    
    # create subset of intersection data frame where overlap occurs
    adjState.overlap.true <- subset(as.data.frame(adjState_clipped_points), Boolean == TRUE)
    
    # capture the subset in a numeric vector
    adjState.overlap.ids <- rownames(adjState.overlap.true)
    
    # adjacent state, make it a data frame
    adjState.df <- as.data.frame(adjState_points)
    
    # for every overlap id and value, add to data frame??
    # need to change out defined AADT value
    for (id in adjState.overlap.ids) {
        id <- as.numeric(id)
      
        id.value <- adjState_points$AADT[id]
        id.coords <- as.numeric(adjState_points@coords[id,])
        
        id.df.row <- c(id.coords, id.value)
        # probably faster to create one smaller data frame to add in one step to addendum df
        id.df.add <- data.frame(x = id.df.row[1], y = id.df.row[2], AADT = id.df.row[3])
        
        add.df <- rbind(add.df, id.df.add)
    }
  } 
  
  # add in the adjacent state points to the original state
  originalState.df <- as.data.frame(originalState_points)
  adjacentStates.df <- add.df
  
  bufferedState.df <- rbind(originalState.df, adjacentStates.df)
  
  # setup parameters for creation of new spatial points data frame
  bufferedState.coords <- as.matrix(bufferedState.df[,1:2])
  bufferedState.AADT <- bufferedState.df["AADT"]
  
  bufferedState.spdf <- SpatialPointsDataFrame(bufferedState.coords, bufferedState.AADT)
  proj4string(bufferedState.spdf) <- proj.str
  names(bufferedState.spdf) <- col.names
  
  return(bufferedState.spdf)
}

# gIntersect attempt

# create the intersection of original border and adj state points
OR_cp_intersect <- gIntersects(WA_border, OR_points, byid = TRUE)
# change the column of the intersection data frame to non-integer
colnames(OR_cp_intersect) <- c("Boolean")

# create subset of intersection data frame where overlap occurs
OR_olp <- subset(as.data.frame(OR_cp_intersect), Boolean == TRUE)
# capture the subset in a numeric vecto
OR_ids <- rownames(OR_olp)

# pull in the data frame with the coordinates, values
OR_pnts_df <- as.data.frame(OR_points)

# capture the values, the coordinates
OR.value <- OR_points$AADT[1]
OR.coords <- as.numeric(OR_points@coords[1,])

# create a row to add the the addendum data frame 
OR.df.row <- c(OR.coords, OR.value)
# add to addendum data frame
OR.df.add <- data.frame(x = OR.df.row[1], y = OR.df.row[2], AADT = OR.df.row[3])

# add addendum data frame to original state
OR.df <- rbind(OR_pnts_df, OR.df.add)

# setup parameters for creation of new spatial points data frame
OR.coords <- as.matrix(OR.df[,1:2])
OR.AADT <- OR.df["AADT"]

# create spatial points data frame
OR.spdf <- SpatialPointsDataFrame(OR.coords, OR.AADT)
proj4string(OR.spdf) <- projection
names(OR.spdf) <- "AADT"



###
###
###




# buffer method from scratch

OR_border <- selectState("OR")
OR_range <- 5969

WA_border <- selectState("WA")
WA_state_buffer <- gBuffer(WA_border, width = OR_range)

vct_path <- paste("/home/sean/traffic/WA_example/", "OR", sep = "")
column_names <- c("AADT")
OR_points <- createSpdf(vct_path, column_names, projection)

vct_path <- paste("/home/sean/traffic/WA_example/", "WA", sep = "")
WA_points <- createSpdf(vct_path, column_names, projection)










for (id in OR_ids) {
  id <- as.integer(id)
  val <- (OR_points$AADT[id])
  
  coords <- as.numeric(OR_points@coords[id,])
}



OR_clipped_points <- gIntersects(WA_border, OR_points, byid = TRUE)
colnames(OR_clipped_points) <- c("Intersects")
OR_cp_df <- as.data.frame(OR_clipped_points)
OR_cp_true <- subset(OR_cp_df, Intersects==TRUE)

# correct

OR_clipped_points <- gIntersection(WA_border, OR_points, byid = TRUE)
tmp <- strsplit(row.names(OR_clipped_points), " ")
OR.id <- (sapply(tmp, "[[", 2))

OR.vals <- c()
for (id in OR.id) {
  val <- OR_points$AADT[as.integer(id)]
  OR.vals <- c(OR.vals, val)
}

WA_aadt <- WA_points$AADT
OR.vals <- OR.vals
as.data.frame(c(WA_aadt, OR.vals))

OR.clipped.point.value <- OR_points$AADT[as.integer(OR.id[1])]
OR.clipped.point.coord <- coordinates(OR_points[as.integer(OR.id[1]),])
























