# starting a session

###Loading R library and packages                                                      

library(geoR)
library(rgdal)
library(gstat)
library(sp)
library(raster)
library(maptools)
library(rasterVis) #Raster visualization
library(rgeos) #GEOS binding for R
library(gtools) #general additional tools
library(colorRamps) #Palette/coloramp for display,contains matlab.like color palette
library(gridExtra)

###### Functions used in this script

#                  #
# data preparation #
#                  #

#setwd("E:/Research/trafficVolume/ID/Pri_Sec_combined2/")

in_dir <- "/home/parmentier/Data/kriging_traffic"
setwd(in_dir)
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
names(OR_spdf) <- c("AADT") #rename column in df
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

#res(OR_rast) <- c(30,30)
res(OR_rast) <- c(1000,1000) #coarse resolution for testing...

proj4string(OR_rast) <- aea
values(OR_rast) <- rep(1, ncell(OR_rast))

# covert to griddf
OR_grid <- as(OR_rast, 'SpatialGridDataFrame')

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

hist(ID$AADT)

## inclusion of log transformation

hist(log(ID$AADT))

## lagged scatterplots

hscat(log(AADT) ~ 1, ID, (3:18) * 100 )
hscat(AADT~1,OR_spdf,(3:18) * 100)
hscat(log(AADT)~1,OR_spdf,(3:18) * 100)

p_hscat <- hscat(AADT~1,OR_spdf,(1:100) * 100)
breaks_hscat <- seq(0,10000,200)
p_hscat <- hscat(AADT~1,OR_spdf,breaks_hscat)

tb_hscat <- hscat(AADT~1,OR_spdf,breaks_hscat,as.table=T)

## variogram plot, cloud
variogram(AADT ~ 1, ID)
vgm1 <- variogram(AADT~1,OR_spdf,width=(3:18) * 100) #need better control range and bin/lag distance here?
vgm1 <- variogram(AADT~1,OR_spdf,cutoff) #need better control range and bin/lag distance here? I think we can use cutoff

##other possiblities?
#-other variogram package translated in gstat
# use vgm or other functions from http://cran.r-project.org/web/packages/gstat/gstat.pdf ?

vg_fit <- fit.variogram(vgm1, model=vgm(1,"Sph",1800,1))
#vg_fit <- fit.variagram(vgm1, model=vgm(1,"Sph",1800,1),fit.ranges=FALSE)
plot(vgm1) #plot shown in the power point.
plot(vg_fit)

#Use krige or krige0 function from gstat package to set the variogram used in the modleing
#x - krige(log(zinc)~1, meuse, meuse.grid, model = m, beta = 5.9)

krige_mod <- krige(AADT~1, OR_spdf, s_spdf, model = vg_fit)#, beta = 5.9)
#we need to have a vg_fit that has a limited range to predict ....


# cloud
varIDcloud <- variogram(log(AADT) ~ 1, ID, cloud = TRUE)
# plot

varID <- variogram(log(AADT) ~1, ID)
varID <- variogram(AADT ~1, OR_spdf,vgm(1,"Sph",1800,1))

plot(varID)

# fitted exponential curve
varID_fit <- fit.variogram(varID, vgm(1, "Exp", 3000, 1))

# fit variogram
variogram <- autofitVariogram(AADT ~ 1, input_data = OR_spdf)

data <- remove.duplicates(OR_spdf)#, zero = dist_val) #spatially sub sample...

kriging_result = autoKrige(AADT ~ 1, input_data = OR_spdf, new_data = OR_grid)
data$AADT <- as.numeric(data$AADT)
kriging_result = autoKrige(AADT ~ 1, input_data = data, new_data = s_spdf,data_variogram=data)

#
# variograms
#

# variogram comparisons

#Read in ArcGIS variogram and compare to R gstat
#
#
#
#
#




