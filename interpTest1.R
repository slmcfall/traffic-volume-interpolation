####################################  INTERPOLATION OF Traffic   #######################################
############################  Script for kriging of traffic    ##############################
#This script tests various variogram models for kriging of traffic volume data.
#AUTHORS: Sean McFall, Benoit Parmentier 
#CREATED ON: 01/12/2015 ? (change here) 
#MODIFIED ON: 01/23/2015            
#Version: 1
#PROJECT: Sean McFall's WM?? (change here)                                     
#################################################################################################

### Loading R library and packages        

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

######################### NOW taking what has beeen learned above ##############
# variogram comparisons

OR_rast <- raster(ex_OR)
res(OR_rast) <- c(100,100)
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

#Read in ArcGIS variogram and compare to R gstat
#
#
#
#
#

#Use cutoff to cut variogram...
vgm1 <- variogram(AADT~1,OR_spdf,cutoff=10000) #need better control range and bin/lag distance here? I think we can use cutoff
plot(vgm1)
vgm1 <- variogram(AADT~1,OR_spdf,width=200,cutoff=10000) #need better control range and bin/lag distance here? I think we can use cutoff
plot(vgm1)
#20km is probably too long, set to 10km for now? could be also defined on distrubtion of distances from input...
vgm1 <- variogram(AADT~1,OR_spdf,width=200,cutoff=20000) #need better control range and bin/lag distance here? I think we can use cutoff
plot(vgm1)

vgm1 <- variogram(AADT~1,OR_spdf,width=200,cutoff=10000) #need better control range and bin/lag distance here? I think we can use cutoff
plot(vgm1) #this is the sample variogram

vgm_model <- vgm(psill=1,model="Exp",range=10000/3,nugget=1) #define theoretical model for the variogram...
#read the documentation to set the values above in some sensible way...

vg_fit <- fit.variogram(vgm1, model=vgm_model) #now fit
#vg_fit <- fit.variagram(vgm1, model=vgm(1,"Sph",1800,1),fit.ranges=FALSE)
plot(vg_fit)

#Use krige or krige0 function from gstat package to set the variogram used in the modleing
#x - krige(log(zinc)~1, meuse, meuse.grid, model = m, beta = 5.9)

data <- remove.duplicates(OR_spdf)#, zero = dist_val) #spatially sub sample...
krige_mod <- try(krige(AADT~1, OR_spdf, newdata= s_spdf, model = vg_fit))#, beta = 5.9)
#we need to have a vg_fit that has a limited range to predict ....

#TODO:
#chekc if cutoff option of fix.value can be used to set the range or cutoff variagrma in the Automap package...
#chekc the error for krige:
#"chfactor.c", line 131: singular matrix in function LDLfactor() Error in predict.gstat(g, newdata = newdata, block = block, nsim = nsim,  :    LDLfactor
#from the google search...:
#Yes. Also near-identical locations in combination with variograms 
#without a nugget effect may caus it.

#The GSTAT online manual suggests to remove "noaverage" or increase "zero"
#to overcome the problem. However, I have not figured out how to do this
#when running gstat under R.

#Either package gstat, but more likely package sp has a function zerodist 
#which will find point pairs with (nearly) zero distances. You may use 
#the indexes it returns to do something with these points (e.g. remove 
#one of them).

#Should I use only a subset of the orignial data or are there other
#alternatives to overcome the problem?

#Yes. An alternative would be to average observations, if they are scalar 
#data.

#Maybe once the grid is created average one observation by grid cell and use this new sample point dataset in the kriging? 