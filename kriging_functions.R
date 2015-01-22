# starting a session

library(geoR)
library(rgdal)
library(gstat)
library(sp)
library(raster)
library(maptools)
library(automap)
library(grDevices)

#                 #
# data collection #
#                 #

# go through all folders
# find the exact name if it's there
# record folder name 
shp_files <- list.files(pattern = "sec_pnts_inter.shp$", recursive = TRUE)

# find file of interest by state

#                  #
# data preparation #
#                  #

# setwd("E:/Research/trafficVolume/ID/Pri_Sec_combined2/")
setwd("/media/sean/Data/Research/trafficVolume/ID/Pri_Sec_combined2/")

# define albers equal area projection
aea <- "+proj=aea +lat_1=38 +lat_2=41 +lat_0=34 +lon_0=-114 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs +ellps=clrk66"

# bring the shpfile into memory
ID <- readOGR(dsn = getwd(), layer = "pr_pnts_inter") #, p4s = NAD27)
# define projection for the shpfile
# proj4string(ID) <- aea

# convert AADT column to data frame
OR_df <- as.data.frame(ID$AADT)

# get shpfile coordinates pair
OR_coords <- coordinates(ID)[,1:2]

#                                             #
# make the shapefile with only two dimensions #
#                                             #

## make spatial points data frame
OR_spdf <- SpatialPointsDataFrame(OR_coords, OR_df)
proj4string(OR_spdf) <- aea

# rename the column name
names(OR_spdf) <- c("AADT")

## make shp file
OR_shp <- writeOGR(OR_spdf, getwd(), "OR_pr_AADT", driver = "ESRI Shapefile", overwrite_layer = TRUE)

#                                              #
# create raster with bounding box of shapefile #
#                                              #

# bounding box of spatial points object
bb_OR <- bbox(OR_spdf)

# create extent object from bbox
ex_OR <- extent(bb_OR)

# split into four chunks

# spatial grid

OR_rast <- raster(ex_OR)
res(OR_rast) <- c(1000,1000)
proj4string(OR_rast) <- aea
#setValues(OR_rast) <- rep(1, ncell(OR_rast))
values(OR_rast) <- rep(1, ncell(OR_rast))

# covert to griddf
OR_grid <- as(OR_rast, 'SpatialGridDataFrame')
s_sgdf <- as(OR_rast, 'SpatialGridDataFrame')

#s_sgdf<-as(r_stack_covar,"SpatialGridDataFrame") #Conversion to spatial grid data frame, only convert the necessary layers!!
s_spdf<-as.data.frame(s_sgdf) #Note that this automatically removes all NA rows
s_spdf<-na.omit(s_spdf) #removes all rows that have na...
coords<- s_spdf[,c('s1','s2')]
coordinates(s_spdf)<-coords
proj4string(s_spdf)<-proj4string(s_sgdf) #Need to assign coordinates...
#           #
# autoKrige #
#           #

# fit variogram
variogram <- autofitVariogram(AADT ~ 1, input_data = OR_spdf)

#this removes duplicates, did not  affectthe results
data <- remove.duplicates(OR_spdf)#, zero = dist_val) #spatially sub sample...

kriging_result = autoKrige(AADT ~ 1, input_data = OR_spdf, new_data = OR_grid)
data$AADT <- as.numeric(data$AADT)
kriging_result = autoKrige(AADT ~ 1, input_data = data, new_data = s_spdf,data_variogram=data)


#data <- remove.duplicates(OR_spdf, zero = dist_val) #spatially sub sample...

##
###
#####
###
##


# spatial pixel data frame

## df to spdf


## will need to figure out how to convert to the correct format

#
# exploratory tools/normalcy
#

## test for normalcy

hist(IDdf$AADT)

## inclusion of log transformation

hist(log(AADT))

## lagged scatterplots

hscat(log(AADT) ~ 1, ID, (3:18) * 100 )

## variogram plot, cloud

# cloud
varIDcloud <- variogram(log(AADT) ~ 1, ID, cloud = TRUE)
# plot
varID <- variogram(log(AADT) ~1, ID)
plot(varID)

# fitted exponential curve
varID_fit <- fit.variogram(varID, vgm(1, "Exp", 3000, 1))


#
# variograms
#

# variogram comparisons

#
#
#
#
#




