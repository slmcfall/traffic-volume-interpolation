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
stateAbbreviations <- list("OR")

main(state.list = stateAbbreviations, column.names = c("AADT"), equation = c("AADT ~ 1"), 
     data.suffix = data.suffix, data.directory = data.directory, resolution = 1000, projection = aea)

#####

main <- function(state.list, column.names, equation, data.suffix, data.directory, resolution, projection) {
  
  for (state.name in state.list){

    adjacent.states <- getAdjStates(state.abbrv = state.name)
    
    state.ranges <- getStateRanges(adj.states = adjacent.states, proj.str = projection, 
                                   col.names  = column.names, equation = equation, 
                                   pnt.suffix = data.suffix, pnt.dir = data.directory, 
                                   rst.res    = resolution
                                   )
    
    buffered.spdf <- getBufferedSpdf(state.abbrv = state.name, pnt.dir = data.directory,
                                     col.names   = column.names, proj.str = projection,
                                     ranges.df   = state.ranges, pnt.suffix = data.suffix,
                                     adj.states  = adjacent.states, data.directory = data.directory
                                     )
    
    log.spdf <- logTransformSpdf(buffered.spdf, column.names, projection)
    validation.spdf.list <- createValidationData(log.spdf, .3)
    training.spdf <- validation.spdf.list[[1]]
    test.spdf <- validation.spdf.list[[2]]
    
    tr.raster <- createRaster(spdf = training.spdf, resolution = resolution )
    tr.avg.spdf <- createAvgSpdf(spdf = training.spdf, rast = tr.raster, col_names = column.names)
    tr.cutoff <- getBufferCutoff(spdf = training.spdf)
    tr.width <- getBufferWidth(rast = tr.raster)
    tr.variogram <- createVariogram(equation = equation, spdf = training.spdf, 
                                    varWidth = tr.width, varCutoff = tr.cutoff)
    # save variogram
    
    tr.grid <- as(tr.raster, 'SpatialGridDataFrame')
    
    tr.krige <- createKrigeLayer(spdf = training.spdf, grid = WA.grid, 
                                  raster = tr.raster, equation = equation, 
                                  variogram = tr.variogram, rst_name = state.name)
    
    # write out rasters, save them as objects too
    writeKrigeLayer(krige.layer = tr.krige, extent.raster = tr.raster, name.prefix = state.name)
    
    tr.surfaces <- getKrigeRasters(krige.layer = tr.krige, extent.raster = tr.raster)
    
    tr.pred.rst <- tr.surfaces[[1]]
    
    tst.validation.dfs <- createValidationDataFrames(spdf = test.spdf, projection = projection,
                                                     krige.surface = tr.pred.rst)
    tr.validation.dfs <- createValidationDataFrames(spdf = training.spdf, projection = projection,
                                                    krige.surface = tr.pred.rst)
    
    tr.err.vals <- getErrorValues(validation.spdfs = tr.validation.dfs, column.names = column.names)
    tst.err.vals <- getErrorValues(validation.spdfs = tst.validation.dfs, column.names = column.names)
    
    output.list <- list(tr.krige, tr.validation.dfs, tst.validation.dfs, tr.variogram, tr.err.vals, tst.err.vals)
    names(output.list) <- c("Krige Object", "Training Validation Table", "Testing Validation Table",
                            "Variogram", "Training Error Values", "Testing Error Values")
    
    saveRDS(output.list, paste(state.name, "_output.rds", sep=""))
  }
}

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



