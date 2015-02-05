############################  INTERPOLATION OF Traffic   ############################
############################  Script for kriging of traffic    ######################
# This script tests various variogram models for kriging of traffic volume data.
# AUTHORS: Sean McFall, Benoit Parmentier 
# CREATED ON: 01/12/2015  
# MODIFIED ON: 01/29/2015            
# Version: 1
# PROJECT: Sean McFall's Traffic Volume Project                                  
#####################################################################################

library(geoR)
library(rgdal)
library(gstat)
library(sp)
library(raster)
library(maptools)
library(rasterVis)  # Raster visualization
library(rgeos)      # GEOS binding for R
library(gtools)     # general additional tools
library(colorRamps) # Palette/coloramp for display,contains matlab-like color palette
library(gridExtra)
library(automap)

### FUNCTIONS USED IN SCRIPT ###

#                  #
# data preparation #
#                  #

# for Sean's windows location
# setwd("E:/Research/trafficVolume/ID/Pri_Sec_combined2/")
# setwd("/media/sean/Data/Research/trafficVolume/ID/Pri_Sec_combined2/")
setwd("W:/slmcfall/traffic_from_desktop/ID/Pri_Sec_combined2/")

# define albers equal area projection
aea <- "+proj=aea +lat_1=38 +lat_2=41 +lat_0=34 +lon_0=-114 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs +ellps=clrk66"

# bring the shpfile into memory
state_sample <- readOGR(dsn = getwd(), layer = "pr_pnts_inter") #, p4s = NAD27)

# convert AADT column to data frame
state_sample_df <- as.data.frame(state_sample$AADT)

# get shpfile coordinates pair
state_sample_coords <- coordinates(state_sample)[,1:2]

#                                             #
# make the shapefile with only two dimensions #
#                                             #

# make spatial points data frame
state_sample_spdf <- SpatialPointsDataFrame(state_sample_coords, state_sample_df)

# define projection
proj4string(state_sample_spdf) <- aea

# define the name as what it is, Annual Average Daily Traffic, AADT
names(state_sample_spdf) <- c("AADT")  

# make shp file
# state_sample_shp <- writeOGR(state_sample_spdf, getwd(), "state_sample_pr_AADT", driver = "ESRI Shapefile", overwrite_layer = TRUE)

#                                              #
# create raster with bounding box of shapefile #
#                                              #

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
# set all raster cells to value of 1, so no N/A problems
values(state_sample_rast) <- rep(1, ncell(state_sample_rast))

# covert to spatial grid data frame
state_sample_grid <- as(state_sample_rast, 'SpatialGridDataFrame')

#                            #
# exploratory tools/normalcy #
#                            #

# test for normalcy
hist(state_sample$AADT)
# inclusion of log transformation
hist(log(state_sample$AADT))

# lagged scatterplots
hscat(log(AADT) ~ 1, state_sample, (3:18) * 100 )
hscat(log(AADT) ~ 1, state_sample_spdf, (3:18) * 100)
hscat(AADT ~ 1, state_sample_spdf, (3:18) * 100)

# p_hscat <- hscat(AADT~1,state_sample_spdf,(1:100) * 100)

# from 0 to 10000 in increments of 200
breaks_hscat <- seq(0,10000,200)
p_hscat <- hscat(AADT~1,state_sample_spdf,breaks_hscat)
# create h plots as table
tb_hscat <- hscat(AADT~1,state_sample_spdf,breaks_hscat,as.table=T)

# variogram plot, cloud
variogram(AADT ~ 1, state_sample)
vgm1 <- variogram(AADT~1,state_sample_spdf,width=(3:18) * 100) #need better control range and bin/lag distance here?

# other possiblities?
# other variogram package translated in gstat
# use vgm or other functions from http://cran.r-project.org/web/packages/gstat/gstat.pdf ?

vg_fit <- fit.variogram(vgm1, model=vgm(1,"Sph",1800,1))
# vg_fit <- fit.variogram(vgm1, model=vgm(1,"Sph",1800,1),fit.ranges=FALSE)
plot(vgm1) # plot shown in the power point
plot(vg_fit)

# Use krige or krige0 function from gstat package to set the variogram used in the modeling
# x - krige(log(zinc)~1, meuse, meuse.grid, model = m, beta = 5.9)

krige_mod <- krige(AADT~1, state_sample_spdf, state_sample_grid, model = vg_fit)#, beta = 5.9)
# we need to have a vg_fit that has a limited range to predict ....

# cloud
varstate_samplecloud <- variogram(log(AADT) ~ 1, state_sample, cloud = TRUE)
# plot, looks pretty good 
varstate_sample <- variogram(log(AADT) ~1, state_sample)

# ERROR THROWN
## varstate_sample w/ "Sph" fit
### Error in as(data, "data.frame") : 
### internal problem in as(): “variogramModel” is(object, "data.frame") is TRUE
### but the metadata asserts that the 'is' relation is FALSE
varstate_sample <- variogram(AADT ~ 1, state_sample_spdf, vgm(1, "Sph", 1800, 1))

# not plotting
plot(varstate_sample)

# fitted exponential curve
varstate_sample_fit <- fit.variogram(varstate_sample, vgm(1, "Exp", 3000, 1))
# fit variogram
variogram <- autofitVariogram(AADT ~ 1, input_data = state_sample_spdf)
# remove duplicate values, seems like two duplicates per state
data <- remove.duplicates(state_sample_spdf)  #, zero = dist_val)  # spatially sub sample...
# krige. horizontal line syndrome...
kriging_result = autoKrige(AADT ~ 1, input_data = state_sample_spdf, new_data = state_sample_grid)
data$AADT <- as.numeric(data$AADT)
kriging_result = autoKrige(AADT ~ 1, input_data = data, new_data = state_sample_grid ,data_variogram=data)

#            #
# using geoR #
#            #

# the first point of semivariance is a negative number...
variogram_geor <- variog(coords = state_sample_coords, data = state_sample_df)
variofit <- variofit(variogram_geor)


#            #
# variograms #
#            #

# variogram comparisons

# state_sample_rast <- raster(state_sample_ex)
# res(state_sample_rast) <- c(100,100)
# proj4string(state_sample_rast) <- aea
# setValues(state_sample_rast) <- rep(1, ncell(state_sample_rast))
values(state_sample_rast) <- rep(1, ncell(state_sample_rast))

# covert to griddf
state_sample_grid <- as(state_sample_rast, 'SpatialGridDataFrame')
# create spatial grid data frame
s_sgdf <- as(state_sample_rast, 'SpatialGridDataFrame')

# s_sgdf <- as(r_stack_covar, "SpatialGridDataFrame") 
# Conversion to spatial grid data frame, only convert the necessary layers!!
s_spdf <- as.data.frame(s_sgdf)  # Note that this automatically removes all NA rows
s_spdf <- na.omit(s_spdf)  # removes all rows that have NA
coords <- s_spdf[,c('s1','s2')]
coordinates(s_spdf) <- coords
proj4string(s_spdf) <- proj4string(s_sgdf)  # Need to assign coordinates...

# Read in ArcGIS variogram and compare to R gstat
#
#
#
#
#

# Use cutoff to cut variogram...

# 10km
vgm1 <- variogram(AADT~1,state_sample_spdf,cutoff=10000) 
plot(vgm1)
# 10km, 200m space between pairs of points
vgm1 <- variogram(AADT~1,state_sample_spdf,width=200,cutoff=10000) 
plot(vgm1)
# 20km, 200m
# 20km is probably too long, set to 10km for now? could be also defined on distrubtion of distances from input...
vgm1 <- variogram(AADT~1,state_sample_spdf,width=200,cutoff=20000) 
plot(vgm1)

vgm1 <- variogram(AADT~1,state_sample_spdf,width=200,cutoff=10000) 
plot(vgm1) #this is the sample variogram

# define theoretical model for the variogram
vgm_model <- vgm(psill = 1 , model = "Exp", range = 10000/3, nugget = 1) 
# read the documentation to set the values above in some sensible way...
# going to set the values based off of what ARCmap does for now
# definitely not a long term solution, as singular model fit still occurs...

vg_fit <- fit.variogram(vgm1, model=vgm_model) #now fit
# vg_fit <- fit.variagram(vgm1, model=vgm(1,"Sph",1800,1),fit.ranges=FALSE)
plot(vg_fit)

# Use krige or krige0 function from gstat package to set the variogram used in the modleing
# x - krige(log(zinc)~1, meuse, meuse.grid, model = m, beta = 5.9)

data <- remove.duplicates(state_sample_spdf)#, zero = dist_val) #spatially sub sample...
krige_mod <- try(krige(AADT~1, state_sample_spdf, newdata = s_spdf, model = vg_fit))#, beta = 5.9)
# we need to have a vg_fit that has a limited range to predict ....

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