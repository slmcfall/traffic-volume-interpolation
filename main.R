# main.R

# define directory with data
data.directory <- "/home/sean/traffic/AADT/"
data.suffix <- "_AADT"

# call functions and libraries
source("~/Documents/trafficVolume/functions.R")
#source("~/Documents/trafficVolume/bufferMethod.R")

# define projection, albers equal area
projection <- "+proj=aea +lat_1=38 +lat_2=41 +lat_0=34 +lon_0=-114 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs"

# polygon for the continental US
cUS <- readOGR(dsn = "/home/sean/Documents/trafficVolume/cUS.shp", layer = "cUS")

# state to be buffered
stateAbbreviation <- "WA"

#####

adjacent.states <- getAdjStates(state.abbrv = stateAbbreviation)

state.ranges <- getStateRanges(adj.states = adjacent.states, proj.str = projection, 
                               col.names  = c("AADT"), equation = c("AADT ~ 1"), 
                               pnt.suffix = data.suffix, pnt.dir = data.directory, 
                               rst.res    = 1000
                               )

buffered.spdf <- getBufferedSpdf(state.abbrv = stateAbbreviation, pnt.dir = data.directory,
                                 col.names   = c("AADT"), proj.str = projection,
                                 ranges.df   = state.ranges, pnt.suffix = data.suffix,
                                 adj.states  = adjacent.states
                                 )

log.spdf <- logTransformSpdf(buffered.spdf, c("AADT"), projection)
validation.spdf.list <- createValidationData(log.spdf, .3)
training.spdf <- validation.spdf.list[[1]]
test.spdf <- validation.spdf.list[[2]]

tr.raster <- createRaster(spdf = training.spdf, resolution = 1000 )
tr.avg.spdf <- createAvgSpdf(spdf = training.spdf, rast = tr.raster, col_names = c("AADT"))
tr.cutoff <- getBufferCutoff(spdf = training.spdf)
tr.width <- getBufferWidth(rast = tr.raster)
tr.variogram <- createVariogram(equation = c("AADT~1"), spdf = traning.spdf, 
                                varWidth = tr.width, varCutoff = tr.cutoff)
# save variogram

tr.grid <- as(tr.raster, 'SpatialGridDataFrame')

tr.krige <- createKrigeLayer(spdf = traning.spdf, grid = WA.grid, 
                              raster = WA.raster, equation = c("AADT ~ 1"), 
                              variogram = WA.variogram, rst_name = "WA")

# write out rasters, save them as objects too
writeKrigeLayer(krige.layer = WA.krige2, extent.raster = tr.raster, name.prefix = stateAbbreviation)

tr.surfaces <- getKrigeRasters(krige.layer = WA.krige2, extent.raster = tr.raster)

tr.pred.rst <- tr.surfaces[[1]]

tst.validation.dfs <- createValidationDataFrames(spdf = test.spdf, projection = projection,
                                                 krige.surface = tr.pred.rst)
tr.validation.dfs <- createValidationDataFrames(spdf = training.spdf, projection = projection,
                                                krige.surface = tr.pred.rst)

tr.err.vals <- getErrorValues(validation.spdfs = tr.validation.dfs, column.names = c("AADT"))
tst.err.vals <- getErrorValues(validation.spdfs = tst.validation.dfs, column.names = c("AADT"))

output.list <- list(WA.krige2, tr.validation.dfs, tst.validation.dfs, tr.variogram, tr.err.vals, tst.err.vals)
names(output.list) <- c("Krige Object", "Training Validation Table", "Testing Validation Table",
                        "Variogram", "Training Error Values", "Testing Error Values")

saveRDS(output.list, paste(stateAbbreviation, "_output.rds", sep=""))


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



