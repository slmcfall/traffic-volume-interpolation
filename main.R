# main.R

# define directory with point data
data.directory <- "/home/sean/traffic/AADT/"
data.suffix <- "_AADT"

# call functions and libraries
source("~/Documents/trafficVolume/functions.R")

# define projection, albers equal area
aea <- "+proj=aea +lat_1=38 +lat_2=41 +lat_0=34 +lon_0=-114 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs"

# polygon for the continental US
cUS <- readOGR(dsn = "/home/sean/Documents/trafficVolume/cUS.shp", layer = "cUS")

# state to be buffered
stateAbbreviations <- list("OR", "NV", "CA")

analysis <- lineprof(main(state.list = stateAbbreviations, column.names = c("AADT"), equation = c("AADT ~ 1"), 
                          data.suffix = data.suffix, data.directory = data.directory, resolution = 1000, 
                          projection = aea))

#####


#
# additional possible analysis of the data
#

# main function (with lineprof, with try/catch stuff)

## for loop over state abbreviations
  # the functions
  # save objects
    # write the rasters
    # data frames of observed and predicted, as objects
    # variogram saved as object
    # error values saved as objects 

############# log transformed data

WAplot <- annotatedplot(WA.variogram)
WA.cutoff
WA.width

# output stuff

# residuals
WA.lm <- lm(WA.observed ~ WA.predicted)
WA.lm.res <- resid(WA.lm)
qqnorm(WA.lm.res)
qqline(WA.lm.res)
hist(WA.lm.res, breaks = 200)

# validation
WA.error <- WA.observed - WA.predicted
rmse(WA.error)
mae(WA.error)
adj.r.sq <- summary(WA.lm)$adj.r.squared
plot(WA.observed~WA.predicted)

# maps

WA.layers <- writeRasters(WA.krige2, WA.raster, "WA2")

min(WA.krige2$var1.pred)
max(WA.krige2$var1.pred)
mean(WA.krige2$var1.pred)
median(WA.krige2$var1.pred)
sd(WA.krige2$var1.pred)



