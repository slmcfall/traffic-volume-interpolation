#
# This script averages data points in a raster cell to remove possibility of singular
# matrix for kriging applications.
# AUTHORS: Sean McFall, Benoit Parmentier 
# CREATED ON: 02/05/2015  
# MODIFIED ON: 02/05/2015            
# Version: 1
# PROJECT: Sean McFall's Traffic Volume Project 
#

# load libraries, setup script
source("libraries_and_wdirectories.R")

# set working directory
set_directory("SWEMCGA18")

# define albers equal area projection
aea <- "+proj=aea +lat_1=38 +lat_2=41 +lat_0=34 +lon_0=-114 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs +ellps=clrk66"

######### BEGIN SCRIPT #########

#                             #
# spdf manipulation, creation #
#                             #

# bring the shpfile into memory
## primary points in OR
state_sample <- readOGR(dsn = getwd(), layer = "pr_pnts_inter") #, p4s = NAD27)
## secondary points in ID
#state_sample <- readOGR(dsn = getwd(), layer = "sec_pnts_all") #, p4s = NAD27)

# convert AADT column to data frame
state_sample_df <- as.data.frame(state_sample$AADT)
# get shpfile coordinates pair
state_sample_coords <- coordinates(state_sample)[,1:2]

# make spatial points data frame
state_sample_spdf <- SpatialPointsDataFrame(state_sample_coords, state_sample_df)
# define projection
proj4string(state_sample_spdf) <- aea
# define the name as what it is
# Annual Average Daily Traffic, AADT
names(state_sample_spdf) <- c("AADT")  

# make shp file
# state_sample_shp <- writeOGR(state_sample_spdf, getwd(), "state_sample_pr_AADT", driver = "ESRI Shapefile", overwrite_layer = TRUE)

#                                                       #
# create raster skeleton with bounding box of shapefile #
#                                                       #

# bounding box of spatial points object
state_sample_bb <- bbox(state_sample_spdf)
# create extent object from bbox
state_sample_ex <- extent(state_sample_bb)
# create raster from extent object
state_sample_rast <- raster(state_sample_ex)

# coarse resolution for testing
res(state_sample_rast) <- c(1000,1000)
# define projection
proj4string(state_sample_rast) <- aea
# set all raster cells to value of -9999, so no N/A problems
values(state_sample_rast) <- rep(-9999, ncell(state_sample_rast))

# covert to spatial grid data frame
state_sample_grid <- as(state_sample_rast, 'SpatialGridDataFrame')

#
# create raster with data points averaged
#

state_sample_raster_avg <- rasterize(state_sample_spdf, state_sample_rast, 
                                     field = "AADT", fun = mean, background = NA)

#
# convert to spdf
#

state_sample_spdf_avg <- rasterToPoints(state_sample_raster_avg, spatial = TRUE)
names(state_sample_spdf_avg) <- c("AADT")

#
# create variogram 
#

state_variogram <- variogram(AADT ~ 1, state_sample_spdf, width = 300, cutoff = 10000)
state_vgm <- vgm(psill = 1 , model = "Exp", range = 10000/3, nugget = 1) 
state_fvariogram <- fit.variogram(state_variogram, model = state_vgm)

#
# vanilla gstat krige
#

state_krige <- krige(AADT ~ 1, state_sample_spdf_avg, state_sample_grid,
                     model = state_fvariogram)

#
# attempt autofitVariogram
#

# fails
## does not take on width, cutoff parameters
# opts <- list(orig.behavior = FALSE)
# state_autoVariogram <- autofitVariogram(AADT ~ 1, input_data = state_sample_spdf,
#                                          width = 300, cutoff = 100, verbose = TRUE,
#                                            miscFitOptions = opts)

# succeeds 
state_afvmod <- afvmod(AADT ~ 1, input_data = state_sample_spdf,
                       width = 300, cutoff = 5000, verbose = TRUE,
                       miscFitOptions = opts)

state_sample_spdf_avg <- remove.duplicates(state_sample_spdf)
state_sample_spdf_avg$AADT <- as.numeric(state_sample_spdf_avg$AADT)

#
# automap autokrige 
#

state_krige <- krige(AADT ~ 1, state_sample_spdf_avg, state_sample_grid,
                     model = state_afvmod$var_model)
  

#
# vector/raster output workaround
#

# prediction layer
krige_pred <- state_krige$var1.pred
krige_pred_rst <- state_sample_rast
names(krige_pred_rst) <- "krige_pred"
clearValues(krige_pred_rst) 
values(krige_pred_rst) <- as.numeric(krige_pred)

# error layer
krige_var <- state_krige$var1.var
krige_var_rst <- state_sample_rast
names(krige_var_rst) <- "krige_var"
clearValues(krige_var_rst) 
values(krige_var_rst) <- as.numeric(krige_var)

#
# write raster 
#


# structure
# str(state_krige)

# change name from krige_pred_rst...terribly redundant
writeRaster(krige_pred_rst, "state_krige_pred_autofit.tif", overwrite = TRUE)

writeRaster(krige_var_rst, "state_krige_var_autofit.tif", overwrite = TRUE)




