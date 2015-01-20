# starting a session

library(geoR)
library(rgdal)
library(gstat)
library(sp)
library(raster)
library(maptools)

#                  #
# data preparation #
#                  #

setwd("E:/Research/trafficVolume/ID/Pri_Sec_combined2/")

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

## make shp file
OR_shp <- writeOGR(OR_spdf, getwd(), "OR_pr_AADT", driver = "ESRI Shapefile", overwrite_layer = TRUE)

#                                              #
# create raster with bounding box of shapefile #
#                                              #

# bounding box of spatial points object
bb_OR <- bbox(OR_spdf)

# create extent object from bbox
ex_OR <- extent(bb_OR)

OR_rast <- raster(ex_OR)
res(OR_rast) <- c(30,30)
proj4string(OR_rast) <- aea

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




