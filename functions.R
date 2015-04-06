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

library(geoR)
library(rgdal)
library(rgeos)
library(gstat)
library(sp)
library(raster)
library(maptools)
library(rasterVis)  
library(rgeos)      
library(gtools)     
library(colorRamps) 
library(gridExtra)
library(automap)
library(tools)
require(lattice)
library(plyr)

createSpdf <- function(vct_path, col_names, proj4) {
  
  #############
  
  # PURPOSE
  # Creates SpatialPointsDataFrame w/ specified attributes
  
  # INPUTS
  # vct_path  = vector name, including file ending
  # vct_dir   = vector directory
  # proj      = proj4string
  # col_names = column names, form of...a vector!
  
  # OUTPUTS
  # vector_spdf      = spatial points data frame that contains desired column/s
  
  #############
  
  # get basename and directory of file
  vct_name <- basename(vct_path)
  vct_dir <- dirname(vct_path)
  
  # bring vector file into memory
  vector <- readOGR(dsn = vct_dir, layer = vct_name)
  
  # get proj4string
  # proj <- "+proj=aea +lat_1=38 +lat_2=41 +lat_0=34 +lon_0=-114 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs"
  proj <- proj4
  
  # convert targeted column into df
  vector_df <- as.data.frame(vector)
  vector_df <- vector_df[col_names]
  
  # get shpfile coordinates pair
  vector_coords <- coordinates(vector)[,1:2]
  
  vector_spdf <- SpatialPointsDataFrame(vector_coords, vector_df)
  # set spdf projection, column names
  proj4string(vector_spdf) <- proj
  names(vector_spdf) <- col_names
  
  return (vector_spdf)
  
}

createRaster <- function(spdf, resolution) {
  
  #############
  
  # PURPOSE
  # creates raster from SpatialPointsDataFrame at specified resolution (integer)
  
  # INPUTS
  # spdf       = SpatialPointsDataFrame
  # resolution = as one number! not a vector
  
  # OUTPUTS
  # rast = spdf converted into raster w/ desired resolution
  
  #############
  
  # bounding box of spatial points object
  spdf_bb <- bbox(spdf)
  # create extent object from bbox
  spdf_ex <- extent(spdf_bb)
  # create raster from extent object
  rast <- raster(spdf_ex)
  
  # coarse resolution for testing
  res(rast) <- c(resolution,resolution)
  # define projection
  proj4string(rast) <- proj4string(spdf)
  # set all raster cells to value of -9999, so no N/A problems
  values(rast) <- rep(-9999, ncell(rast))
  
  return (rast)
  
}

createAvgSpdf <- function(spdf, rast, col_names) {
  
  ############
  
  # PURPOSE
  # creates SpatialPointsDataFrame with values averaged
  
  # INPUTS
  # spdf      = SpatialPointsDataFrame
  # rast      = Extent raster
  # col_names = column names, as vector
  
  # OUTPUTS
  # spdf_avg  = spdf with values averaged
  
  #############
  
  # create raster with data points averaged
  raster_avg <- rasterize(spdf, rast, field = col_names, fun = mean, background = NA)
  
  # convert to spdf
  spdf_avg <- rasterToPoints(raster_avg, spatial = TRUE)
  names(spdf_avg) <- col_names
  
  return (spdf_avg)
  
}

createVariogram <- function(equation, spdf, varWidth, varCutoff) {
  
  #############
  
  # PURPOSE
  # creates variogram through modded autoFitVariogram function
  
  # INPUTS
  # equation = equation for the variogram model (ex: AADT ~ 1)
  # spdf     = SpatialPointsDataFrame
  # width    =  width between pairs of points
  # cutoff   = range of the variogram, essentially
  
  # OUTPUTS
  # variogram = variogram for kriging
  
  #############
  
  # original behavior MUST be ignored else cutoff and width parameters will have no effect! 
  opts <- list(orig.behavior = FALSE)
  
  variogram <- afvmod(as.formula(equation), input_data = spdf, width = varWidth, cutoff = varCutoff, 
                      verbose = TRUE, miscFitOptions = opts)
  
  return (variogram)
}

getBufferCutoff <- function(spdf) {
  spdf_bbox <- bbox(spdf)
  x1 <- spdf_bbox[1]
  x2 <- spdf_bbox[2]
  y1 <- spdf_bbox[3]
  y2 <- spdf_bbox[4]
  
  spdf_bbox_dist <- sqrt((x1-x2)^2 + (y1-y2)^2)
  spdf_range <- spdf_bbox_dist * .01
  
  return(spdf_range) 
}

getBufferWidth <- function(rast) {
  width <- res(rast)[1] / 2
  return (width)
}

createKrigeLayer <- function(spdf, grid, raster, equation, variogram, rst_name) {
  
  spdf <- remove.duplicates(spdf)
  spdf$AADT <- as.numeric(spdf$AADT)
  
  #
  # automap autokrige 
  #
  
  krige <- krige(formula = as.formula(equation), spdf, grid,
                       model = variogram$var_model)
  
  return(krige)
  
  #
  # vector/raster output workaround 
  #
  
  # prediction layer
  krige_pred <- krige$var1.pred
  krige_pred_rst <- raster
  names(krige_pred_rst) <- "krige_pred"
  clearValues(krige_pred_rst) 
  values(krige_pred_rst) <- as.numeric(krige_pred)
  
  # error layer
  krige_var <- krige$var1.var
  krige_var_rst <- raster
  names(krige_var_rst) <- "krige_var"
  clearValues(krige_var_rst) 
  values(krige_var_rst) <- as.numeric(krige_var)
  
  writeRaster(krige_pred_rst, paste(rst_name,"pred_autofit.tif", sep=""), overwrite = TRUE)
  
  writeRaster(krige_var_rst, paste(rst_name,"var_autofit.tif", sep=""), overwrite = TRUE)
}

writeKrigeLayer <- function(krige.layer, extent.raster, name.prefix) {
  
  # prediction layer
  krige_pred <- krige.layer$var1.pred
  krige_pred_rst <- extent.raster
  names(krige_pred_rst) <- "krige_pred"
  clearValues(krige_pred_rst) 
  values(krige_pred_rst) <- as.numeric(krige_pred)
  
  # error layer
  krige_var <- krige.layer$var1.var
  krige_var_rst <- extent.raster
  names(krige_var_rst) <- "krige_var"
  clearValues(krige_var_rst) 
  values(krige_var_rst) <- as.numeric(krige_var)
  
  writeRaster(krige_pred_rst, paste(name.prefix,"pred_autofit.tif", sep=""), overwrite = TRUE)
  
  writeRaster(krige_var_rst, paste(name.prefix,"var_autofit.tif", sep=""), overwrite = TRUE)
}

getKrigeRasters<- function(krige.layer, extent.raster) {
  
  # prediction layer
  krige_pred <- krige.layer$var1.pred
  krige_pred_rst <- extent.raster
  names(krige_pred_rst) <- "krige_pred"
  clearValues(krige_pred_rst) 
  values(krige_pred_rst) <- as.numeric(krige_pred)
  
  # error layer
  krige_var <- krige.layer$var1.var
  krige_var_rst <- extent.raster
  names(krige_var_rst) <- "krige_var"
  clearValues(krige_var_rst) 
  values(krige_var_rst) <- as.numeric(krige_var)
  
  return(list(krige_pred_rst, krige_var_rst))
}

doKrige <- function(shapefile, columnname, resolution, equation, width, cutoff, outputname) {
  
  spdf <- createSpdf(shapefile, columnname)
  
  raster <- createRaster(spdf, resolution)
  
  grid <- as(raster, 'SpatialGridDataFrame')
  
  spdf_avg <- createAvgSpdf(spdf, raster, columnname)
  
  variog <- createVariogram(equation, spdf, width, cutoff)
  
  kr <- createKrigeLayer(spdf, grid, raster, equation, variog, outputname)
  
  saveObjects(variogram, kr, "RI")
}

saveObjects <- function(variogram, krigeoutput, output_string) {
  
  object_list <- list(variogram, krigeoutput)
  
  save(object_list, output_string)
  
}

getRange <- function(variogram) {
  range <- variogram$var_model$range[2]
  return (range)
}


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

getBufferedSpdf <- function(state.abbrv, pnt.dir, col.names, proj.str, 
                            ranges.df, pnt.suffix, adj.states, data.directory) {
  
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
    vct_path <- paste(data.directory, "/", stateAbbrv, "_AADT", sep = "")
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

# Function that returns Root Mean Squared Error
rmse <- function(error)
{
  sqrt(mean(error^2, na.rm = TRUE))
}

# Function that returns Mean Absolute Error
mae <- function(error)
{
  mean(abs(error), na.rm = TRUE)
}

writeRasters <- function(krige.layer, extent.raster, output.name) {
  
  # prediction layer
  krige_pred <- krige.layer$var1.pred
  krige_pred_rst <- extent.raster
  names(krige_pred_rst) <- "krige_pred"
  clearValues(krige_pred_rst) 
  values(krige_pred_rst) <- as.numeric(krige_pred)
  
  # error layer
  krige_var <- krige.layer$var1.var
  krige_var_rst <- extent.raster
  names(krige_var_rst) <- "krige_var"
  clearValues(krige_var_rst) 
  values(krige_var_rst) <- as.numeric(krige_var)
  
  writeRaster(krige_pred_rst, paste(output.name, "_pred.tif", sep=""), overwrite = TRUE)
  writeRaster(krige_var_rst, paste(output.name, "_var.tif",  sep=""), overwrite = TRUE)
  
  return(list(krige_pred_rst, krige_var_rst))
}

# need to fix inclusion of column.names
getErrorValues <- function(validation.spdfs, column.names) {
  
  observed.spdf <- tst.validation.dfs[[1]]
  predicted.spdf <- tst.validation.dfs[[2]]
  
  obs.vals <- as.data.frame(observed.spdf)
  obs.vals <- obs.vals[column.names]
  obs.vals <- data.matrix(obs.vals)
  
  pred.vals <- predicted.spdf$AADT
  
  resid.vals <- obs.vals - pred.vals
  
  mae.val <- mae(resid.vals)
  rmse.val <- rmse(resid.vals)
  
  error.list <- list(mae.val, rmse.val)
  names(error.list) <- c("MAE", "RMSE")
  
  return(error.list)
  
}
# this needs a better function name
createValidationData <- function(spdf, testingProportion) {
  n <- nrow(spdf)
  ns <- n - round(n * testingProportion)
  nv <- n - ns
  
  index.training <- sample(nrow(spdf), size = ns, replace = FALSE)
  index.testing <- setdiff(1:nrow(spdf), index.training)
  
  spdf.training <- spdf[index.training,]
  spdf.testing <- spdf[index.testing,]
  
  return(list(spdf.training, spdf.testing))
}

# this needs a better function name
createValidationDataFrames <- function(spdf, projection, krige.surface) {
  
  # get points where values are compared
  spdf.points <- SpatialPoints(coordinates(spdf))
  proj4string(spdf.points) <- projection
  
  # extract values from kriged raster surface
  extracted.values <- extract(krige.surface, spdf.points)
  
  observed.values.df <- spdf
  predicted.values.df <- SpatialPointsDataFrame(coords = as.data.frame(spdf@coords),
                                                data   = as.data.frame(extracted.values))
  names(predicted.values.df) <- c("AADT") 
    
  output.list <- list(observed.values.df, predicted.values.df)
  names(output.list) <- c("Observed", "Predicted")
  
  #
  # rename column names
  #
  
  return(output.list)
}

logTransformSpdf <- function(spdf, column.names, projection) {
  
  # spdf = spatialpointsdataframe
  # column.name = vector with attribute names
  # projection = desired proj4string
  
  log.coordinates <- as.data.frame(spdf@coords)
  log.values <- as.data.frame(spdf)
  log.values <- log(log.values[column.names])
  
  log.spdf <- SpatialPointsDataFrame(coords = log.coordinates, data = log.values)
  names(log.spdf) <- c(column.names)
  proj4string(log.spdf) <- projection
  
  return(log.spdf)
}

# next three functions, credit to Michael Koohafkan

annotatedplot <- function(krigeobj) {
  annotatedplot <- xyplot(gamma ~ dist, data = krigeobj$exp_var, panel = automap:::autokrige.vgm.panel,
                          #labels = as.character(krigeobj$exp_var$np * 0), 
                          labels = "",
                          shift = 0.03, model = krigeobj$var_model,
                          direction = c(krigeobj$exp_var$dir.hor[1], krigeobj$exp_var$dir.ver[1]),
                          ylim = c(min(0, 1.04 * min(krigeobj$exp_var$gamma)), 1.04 * max(krigeobj$exp_var$gamma)),
                          xlim = c(0, 1.04 * max(krigeobj$exp_var$dist)), xlab = "Distance", ylab = "Semi-variance",
                          main = "Experimental variogram and fitted variogram model", mode = "direct") 
  return(annotatedplot)
}

gg_variogram <- function(krigeobj) {
  ggplot(krigeobj$exp_var, aes(x=dist, y=gamma)) + geom_point(color='blue') + geom_text(aes(label=np)) 
  
}

afvmod <- function(formula, input_data, model = c("Sph", "Exp", "Gau", "Ste"),
                   kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10), fix.values = c(NA,NA,NA),
                   verbose = FALSE, GLS.model = NA, start_vals = c(NA,NA,NA), 
                   miscFitOptions = list(), fit.method = 7, ...) {
  # This function automatically fits a variogram to input_data
  
  # Check for anisotropy parameters
  if('alpha' %in% names(list(...))) warning('Anisotropic variogram model fitting not supported, see the documentation of autofitVariogram for more details.')
  
  # Take the misc fit options and overwrite the defaults by the user specified ones
  miscFitOptionsDefaults = list(merge.small.bins = TRUE, min.np.bin = 5, 
                                orig.behavior = TRUE, 
                                num.bins = NA, equal.width.bins = FALSE,
                                equal.np.bins = FALSE)
  miscFitOptions = modifyList(miscFitOptionsDefaults, miscFitOptions)
  
  # Create boundaries
  longlat = !is.projected(input_data)
  if(is.na(longlat)) longlat = FALSE
  diagonal = spDists(t(bbox(input_data)), longlat = longlat)[1,2] * 0.35 # 0.35 times the length of the central axis through the area
  
  ### BEGIN MODIFICATIONS ###
  
  if(miscFitOptions[["orig.behavior"]]){
    if(verbose) cat ("Boundaries as defined by original autofitVariogram...\n\n")
    # compute boundaries the old way
    boundaries = c(2,4,6,9,12,15,25,35,50,65,80,100) * diagonal/100         # Boundaries for the bins in km
    # If you specifiy a variogram model in GLS.model the Generelised least squares sample variogram is constructed
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, boundaries = boundaries, ...)
    } else {
      if(verbose) cat("Calculating GLS sample variogram\n")
      g = gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls=1))
      experimental_variogram = variogram(g, boundaries = boundaries, ...)
    }
    # request by Jon Skoien
    if(miscFitOptions[["merge.small.bins"]]) {
      # bin the old way
      if(verbose) cat("Checking if any bins have less than 5 points, merging bins when necessary...\n\n")
      while(TRUE) {
        if(length(experimental_variogram$np[experimental_variogram$np < miscFitOptions[["min.np.bin"]]]) == 0 | length(boundaries) == 1) break
        boundaries = boundaries[2:length(boundaries)]                 
        if(!is(GLS.model, "variogramModel")) {
          experimental_variogram = variogram(formula, input_data, boundaries = boundaries, ...)
        } else {
          experimental_variogram = variogram(g, boundaries = boundaries, ...)
        }
      } 
    }
    ### equal-width bins (approximately, variogram does its own binning too) ###
  } else if(miscFitOptions[["equal.width.bins"]]){
    if(verbose) cat("Using equal-width bins...\n")
    if('width' %in% names(list(...))) stop('Cannot pass width when equal.width.bins = TRUE. Supply "init.width" in miscFitOptions instead.')
    # replace diagonal with cutoff, if provided
    if('cutoff' %in% names(list(...))) diagonal <- list(...)[['cutoff']]
    # user must supply either width or num.bins
    if(!'init.width' %in% names(miscFitOptions)){
      if(is.na(miscFitOptions[['num.bins']])) stop('when equal.width.bins = TRUE, user must supply either init.width or num.bins as well.')
      width <- diagonal/miscFitOptions[['num.bins']]
      if(verbose) cat("initial width not provided. Calculating using num.bins.\n")
    } else {
      width <- miscFitOptions[['init.width']]
      if(verbose) cat("initial width provided.\n")
    }
    # get the experimental variogram
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, width = width, ...)
    } else {
      experimental_variogram = variogram(g, width = width, ...)
    }
    # merge small bins if requested
    if(miscFitOptions[['merge.small.bins']]){
      if(verbose) cat("Checking if any bins have less than ", miscFitOptions[["min.np.bin"]], " points, merging bins when necessary...\n")
      iter <- 0
      maxiter <- 1000
      while(TRUE){
        if(!any(experimental_variogram$np < miscFitOptions[["min.np.bin"]])) break  		
        # increase width by 10% and try again
        width <- width*1.1
        if(!is(GLS.model, "variogramModel")) {
          experimental_variogram = variogram(formula, input_data, width = width, ...)
        } else {
          experimental_variogram = variogram(g, width = width, ...)
        }
        iter <- iter + 1
        if(iter > maxiter){
          cat('maximum number of interations reached. Try decreasing min.np.bin or init.width.\n\n')
          break
        }
      }
    }
    ### equal observation count bins ###
  } else if(miscFitOptions[["equal.np.bins"]]){
    if(verbose) cat("Using bins of equal observation counts...\n")
    if('boundaries' %in% names(list(...))) stop('Cannot pass boundaries when equal.np.bins is TRUE. Pass num.bins or min.np.bin instead.')
    # replace diagonal with cutoff, if provided
    if('cutoff' %in% names(list(...))) diagonal <- list(...)[['cutoff']]
    # get a sorted list of distances 
    dists <- sort(spDists(input_data))
    # apply the cutoff
    dists <- dists[dists < diagonal & dists > 0]
    # split the data into bins based on number of observations
    if(is.na(miscFitOptions[['num.bins']])){
      # compute number of bins based on the minimum number of observations per bin
      miscFitOptions[['num.bins']] <- floor(0.5*length(dists)/miscFitOptions[['min.np.bin']])
      if(verbose) cat("num.bins not supplied. Setting num.bins =", miscFitOptions[['num.bins']], 'based on min.np.bin.\n')
    }
    cat("checking bins, decreasing num.bins if necessary... \n")
    while(TRUE){
      # compute interval based on the number of bins
      interval <- length(dists)/miscFitOptions[['num.bins']]
      # define boundaries
      boundaries <- rep(NA, miscFitOptions[['num.bins']])
      for(i in 1:miscFitOptions[['num.bins']]){
        boundaries[i] <- dists[round(i*interval)]
      }
      if(length(boundaries == length(unique(boundaries)))) break
      # reduce number of bins
      miscFitOptions[['num.bins']] <- miscFitOptions[['num.bins']] - 1 
    }
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, boundaries = boundaries, ...)
    } else {
      if(verbose) cat("Calculating GLS sample variogram\n")
      g = gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls=1))
      experimental_variogram = variogram(g, boundaries = boundaries, ...)
    }
  } else {
    # default behavior of variogram
    cat("No binning action specified in miscFitOptions.\n\n")
    if(!is(GLS.model, "variogramModel")) {
      experimental_variogram = variogram(formula, input_data, ...)
    } else {
      if(verbose) cat("Calculating GLS sample variogram\n")
      g = gstat(NULL, "bla", formula, input_data, model = GLS.model, set = list(gls=1))
      experimental_variogram = variogram(g, ...)
    }
  }
  ### END MODIFICATIONS ###
  
  # set initial values
  if(is.na(start_vals[1])) {  # Nugget
    initial_nugget = min(experimental_variogram$gamma)
  } else {
    initial_nugget = start_vals[1]
  }
  if(is.na(start_vals[2])) { # Range
    initial_range = 0.1 * diagonal   # 0.10 times the length of the central axis through the area
  } else {
    initial_range = start_vals[2]
  }
  if(is.na(start_vals[3])) { # Sill
    initial_sill = mean(c(max(experimental_variogram$gamma), median(experimental_variogram$gamma)))
  } else {
    initial_sill = start_vals[3]
  }
  
  # Determine what should be automatically fitted and what should be fixed
  # Nugget
  if(!is.na(fix.values[1]))
  {
    fit_nugget = FALSE
    initial_nugget = fix.values[1]
  } else
    fit_nugget = TRUE
  
  # Range
  if(!is.na(fix.values[2]))
  {
    fit_range = FALSE
    initial_range = fix.values[2]
  } else
    fit_range = TRUE
  
  # Partial sill
  if(!is.na(fix.values[3]))
  {
    fit_sill = FALSE
    initial_sill = fix.values[3]
  } else
    fit_sill = TRUE
  
  getModel = function(psill, model, range, kappa, nugget, fit_range, fit_sill, fit_nugget, fit.method, verbose)
  {
    if(verbose) debug.level = 1 else debug.level = 0
    if(model == "Pow") {
      warning("Using the power model is at your own risk, read the docs of autofitVariogram for more details.")
      if(is.na(start_vals[1])) nugget = 0
      if(is.na(start_vals[2])) range = 1    # If a power mode, range == 1 is a better start value
      if(is.na(start_vals[3])) sill = 1
    }
    if(fit.method==5){
      #fit.variogram.reml(formula, locations, data, model, debug.level = 1, set, degree = 0)
      #	obj = try(fit.variogram.reml(experimental_variogram,
      #					        model=vgm(psill=psill, model=model, range=range,
      #							nugget=nugget,kappa = kappa),
      #			  fit.ranges = c(fit_range), fit.sills = c(fit_nugget, fit_sill),
      #			  debug.level = 0, fit.method=fit.method), TRUE)		
    }
    else if(fit.method==8){
      #fit.variogram.gls(formula, data, model, maxiter = 30,
      #eps = .01, trace = TRUE, ignoreInitial = TRUE, cutoff = Inf,
      #plot = FALSE
      #obj = try()
    }
    else{
      obj = try(fit.variogram(experimental_variogram,
                              model=vgm(psill=psill, model=model, range=range,
                                        nugget=nugget,kappa = kappa),
                              fit.ranges = c(fit_range), fit.sills = c(fit_nugget, fit_sill),
                              debug.level = 0, fit.method=fit.method), TRUE)
    }
    if("try-error" %in% class(obj)) {
      #print(traceback())
      warning("An error has occured during variogram fitting. Used:\n", 
              "\tnugget:\t", nugget, 
              "\n\tmodel:\t", model, 
              "\n\tpsill:\t", psill,
              "\n\trange:\t", range,
              "\n\tkappa:\t",ifelse(kappa == 0, NA, kappa),
              "\n  as initial guess. This particular variogram fit is not taken into account. \nGstat error:\n", obj)
      return(NULL)
    } else return(obj)
  }
  
  
  # Automatically testing different models, the one with the smallest sums-of-squares is chosen
  test_models = model
  SSerr_list = c()
  vgm_list = list()
  counter = 1
  
  for(m in test_models) {
    if(m != "Mat" && m != "Ste") {        # If not Matern and not Stein
      model_fit = getModel(initial_sill - initial_nugget, m, initial_range, kappa = 0, initial_nugget, fit_range, fit_sill, fit_nugget, verbose = verbose, fit.method=fit.method)
      if(!is.null(model_fit)) {       # skip models that failed
        vgm_list[[counter]] = model_fit
        SSerr_list = c(SSerr_list, attr(model_fit, "SSErr"))}
      counter = counter + 1
    } else {                 # Else loop also over kappa values
      for(k in kappa) {
        model_fit = getModel(initial_sill - initial_nugget, m, initial_range, k, initial_nugget, fit_range, fit_sill, fit_nugget, verbose = verbose, fit.method=fit.method)
        if(!is.null(model_fit)) {
          vgm_list[[counter]] = model_fit
          SSerr_list = c(SSerr_list, attr(model_fit, "SSErr"))}
        counter = counter + 1
      }
    }
  }
  
  # Check for negative values in sill or range coming from fit.variogram
  # and NULL values in vgm_list, and remove those with a warning
  strange_entries = sapply(vgm_list, function(v) any(c(v$psill, v$range) < 0) | is.null(v))
  if(any(strange_entries)) {
    if(verbose) {
      print(vgm_list[strange_entries])
      cat("^^^ ABOVE MODELS WERE REMOVED ^^^\n\n")
    }
    warning("Some models where removed for being either NULL or having a negative sill/range/nugget, \n\tset verbose == TRUE for more information")
    SSerr_list = SSerr_list[!strange_entries]
    vgm_list = vgm_list[!strange_entries]
  }
  
  if(verbose) {
    cat("Selected:\n")
    print(vgm_list[[which.min(SSerr_list)]])
    cat("\nTested models, best first:\n")
    tested = data.frame("Tested models" = sapply(vgm_list, function(x) as.character(x[2,1])), 
                        kappa = sapply(vgm_list, function(x) as.character(x[2,4])), 
                        "SSerror" = SSerr_list)
    tested = tested[order(tested$SSerror),]
    #####
    ##### add in line to output this data frame into a .csv
    #####
    print(tested)
  }
  
  result = list(exp_var = experimental_variogram, var_model = vgm_list[[which.min(SSerr_list)]], sserr = min(SSerr_list))
  class(result) = c("autofitVariogram","list")    
  
  return(result)
}

#
# example 
#

#file <- paste(getwd(), "/sec_pnts_inter", sep="")
#doKrige(file, c("AADT"), 1000, c("AADT~1"), 300, 5000, "RI_sec_" )

# GA is acting a fool
# ME is acting a fool
# IA is huge
# NC is huge
# NH is acting a fool











