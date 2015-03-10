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

# function sequence
setwd("/home/sean/Documents/trafficVolume/")
source("autofitVariogram_mk.R")
setwd("/media/sean/Data/Research/trafficVolume/Buffers/RI/Sec_Buffers/")

getProj <- function(vct_path) {
  
  #############
   
  # PURPOSE
  # Extracts proj4string from a shapefiles .prj file
  
  # INPUTS
  # vct_path  = vector name, including file ending
  
  # OUTPUTS
  # proj4string = proj4string
  
  #############
  
  # run system command
  # outputs to temp txt file
  system(paste("sudo bash -c gdalsrsinfo", vct_path, "> proj4str.txt", sep = " "))
  
  # get second line of txt file
  proj4string <- readLines("proj4str.txt", n = 2)
  
  # save second line as variable
  proj4string <- substr(proj4string[2], 10, 200)
  
  # remove first and last quotes from the string
  proj4string <- substr(proj4string, 2, nchar(proj4string)-1)
  
  return(proj4string)
}

createSpdf <- function(vct_path, col_names) {
  
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
  proj <- "+proj=aea +lat_1=38 +lat_2=41 +lat_0=34 +lon_0=-114 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs"

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
  
  #############
  
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

createVariogram <- function(equation, spdf, width, cutoff) {
  
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
  
  opts <- list(orig.behavior = FALSE)
  
  # equation <- lapply(equation, as.formula)
  
  print(equation)
  
  variogram <- afvmod(as.formula(equation), input_data = spdf, width = width, cutoff = cutoff, 
                      verbose = TRUE) # , miscFitOptions = opts)
  
  return (variogram)
  
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

saveObjects <- function(variogram, krigeoutput, output_string) {
  
  object_list <- list(variogram, krigeoutput)
  
  save(object_list, output_string)
  
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

#
# example 
#

file <- paste(getwd(), "/sec_pnts_inter", sep="")
doKrige(file, c("AADT"), 1000, c("AADT~1"), 300, 5000, "RI_sec_" )

spdf <- createSpdf(paste(getwd(),"/sec_pnts_all",sep=""), c("AADT"))

raster <- createRaster(spdf, 100)

grid <- as(raster, 'SpatialGridDataFrame')

spdf_avg <- createAvgSpdf(spdf, raster, c("AADT"))

variogram <- afvmod(AADT~1, input_data = spdf )
variogram2 <- afvmod(AADT~1, input_data = spdf, width = 300, cutoff = 5000)


#
# for loop to create range 
#
states <- c()
ranges <- c()

ready_files = list.files(path = "~/traffic/data", pattern = "*shp*", full.names = TRUE)

for (file in ready_files) {
  file <- file_path_sans_ext(file)
  state_abrv <- substring(basename(file), 1 , 2)
  states <- c(states, state_abrv)
  
   
  spdf <- createSpdf(file, c("AADT")
  raster <- createRaster(spdf, 500)
  grid <- as(raster, 'SpatialGridDataFrame')
  spdf_avg <- createAvgSpdf(spdf, raster, c("AADT"))
  variogram <- createVariogram(c("AADT ~ 1"), spdf, 300, 5000) )
  
  range <- variogram$var_model$range[2]
  ranges <- c(ranges, range)
  #a <- c(a, file)
}

#
# single case use, how to get range by itself
#

AL <- "/home/sean/traffic/data/AL_ready"

projection <- "+proj=aea +lat_1=38 +lat_2=41 +lat_0=34 +lon_0=-114 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs"
spdf <- createSpdf(AL, c("AADT"))
raster <- createRaster(spdf, 500)
grid <- as(raster, 'SpatialGridDataFrame')
spdf_avg <- createAvgSpdf(spdf, raster, c("AADT"))
variogram <- createVariogram(c("AADT ~ 1"), spdf, 300, 5000)












