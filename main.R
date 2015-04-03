# main.R

# define directory with data
data.directory <- "/home/sean/traffic/AADT/"
data.suffix <- "_AADT"

# call functions and libraries
source("~/Documents/trafficVolume/functions.R")
source("~/Documents/trafficVolume/bufferMethod.R")

# define projection, albers equal area
projection <- "+proj=aea +lat_1=38 +lat_2=41 +lat_0=34 +lon_0=-114 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs"

# polygon for the continental US
cUS <- readOGR(dsn = "/home/sean/Documents/trafficVolume/cUS.shp", layer = "cUS")

# state to be buffered
stateAbbreviation <- "WA"


#####

adjacent.states <- getAdjStates(state.abbrv = stateAbbreviation)

state.ranges <- getStateRanges(adj.states = adjacent.states, 
                               proj.str   = projection, 
                               col.names  = c("AADT"), 
                               equation   = c("AADT ~ 1"), 
                               pnt.suffix = data.suffix, 
                               pnt.dir    = data.directory, 
                               rst.res    = 1000
                               )

buffered.spdf <- getBufferedSpdf(state.abbrv = stateAbbreviation,
                                 pnt.dir     = data.directory,
                                 col.names   = c("AADT"),
                                 proj.str    = projection,
                                 ranges.df   = state.ranges,
                                 pnt.suffix  = data.suffix,
                                 adj.states  = adjacent.states
                                 )

log.spdf <- logTransformSpdf(buffered.spdf, c("AADT"), projection)
validation.spdf.list <- createValidationData(buffered.spdf, .3)
training.spdf <- validation.spdf.list[[1]]
test.spdf <- validation.spdf.list[[2]]

tr.raster <- createRaster(spdf = training.spdf, resolution = 1000 )
tr.avg.spdf <- createAvgSpdf(spdf = training.spdf, rast = tr.raster, col_names = c("AADT"))
tr.cutoff <- getBufferCutoff(spdf = training.spdf)
tr.width <- getBufferWidth(rast = tr.raster)
tr.variogram <- createVariogram(equation = c("AADT~1"), spdf = traning.spdf, varWidth = tr.width, varCutoff = tr.cutoff)
tr.grid <- as(tr.raster, 'SpatialGridDataFrame')

tr.krige <- createKrigeLayer(spdf = traning.spdf, grid = WA.grid, 
                              raster = WA.raster, equation = c("AADT ~ 1"), 
                              variogram = WA.variogram, rst_name = "WA")

writeKrigeLayer(krige.layer = tr.krige, raster = tr.raster, name.prefix = stateAbbreviation)

tr.surfaces <- getKrigeRasters(krige.layer = WA.krige2, extent.raster = tr.raster)



############# log transformed data

## get points for testing data
WA.log.sp.tst <- SpatialPoints(coordinates(tst.spdf))
proj4string(traning.spdf) <- projection
WA.log.sp.tr <- SpatialPoints(coordinates(traning.spdf))
proj4string(traning.spdf) <- projection

## pull out layers

# prediction layer
krige_pred <- WA.krige2$var1.pred
krige_pred_rst <- WA.raster
names(krige_pred_rst) <- "krige_pred"
clearValues(krige_pred_rst) 
values(krige_pred_rst) <- as.numeric(krige_pred)

# error layer
krige_var <- WA.krige2$var1.var
krige_var_rst <- WA.raster
names(krige_var_rst) <- "krige_var"
clearValues(krige_var_rst) 
values(krige_var_rst) <- as.numeric(krige_var)

## extract testing points calues from training raster
WA.log.tst.ext <- extract(krige_pred_rst, WA.log.sp.tst)
WA.log.tr.ext <- extract(krige_pred_rst, WA.log.sp.tr)

## get observed and predicted values
WA.spdf.obs <- tst.spdf
WA.spdf.pred <- SpatialPointsDataFrame(coords = as.data.frame(WA.log.sp.tst@coords), 
                                       data   = as.data.frame(WA.log.tst.ext))

WA.predicted <- WA.log.tst.ext
WA.observed <- tst.spdf$AADT
WA.residuals <- WA.observed - WA.predicted

WA.predicted.tr <- WA.log.tr.ext
WA.observed.tr <- traning.spdf$AADT
WA.residuals <- WA.observed.tr - WA.predicted.tr
mae(WA.residuals)
rmse(WA.residuals)

## residuals as spdf
#WA.spdf.res <- SpatialPointsDataFrame(coords = as.data.frame(WA.sp.tst@coords),
#                                     data   = as.data.frame(WA.residuals))

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


################ END log transformed 

# manual creation of washington state buffered points
# WA.spdf <- createSpdf(vct_path = "/home/sean/traffic/WA_example/WA_buffer_points", col_names = c("AADT"), proj4 = projection)
WA.raster <- createRaster(spdf = traning.spdf, resolution = 1000 )
WA.avg.spdf <- createAvgSpdf(spdf = traning.spdf, rast = WA.raster, col_names = c("AADT"))
WA.cutoff <- getBufferCutoff(spdf = traning.spdf)
WA.width <- getBufferWidth(rast = WA.raster)
WA.variogram <- createVariogram(equation = c("AADT~1"), spdf = traning.spdf, varWidth = WA.width, varCutoff = WA.cutoff)
WA.grid <- as(WA.raster, 'SpatialGridDataFrame')

WA.krige2 <- createKrigeLayer(spdf = traning.spdf, grid = WA.grid, 
                 raster = WA.raster, equation = c("AADT ~ 1"), 
                 variogram = WA.variogram, rst_name = "WA")

############# NON TRANSFORMED

# get points for testing data
WA.sp.tst <- SpatialPoints(coordinates(WA.spdf.tst))
proj4string(WA.sp.tst) <- projection

# extract testing points values from training raster
WA.tst.extr <- extract(krige_pred_rst, WA.sp.tst)

# observed and predicted as spatial point data frames
WA.spdf.obs <- WA.spdf.tst
WA.spdf.pred <- SpatialPointsDataFrame(coords = as.data.frame(WA.sp.tst@coords), 
                                       data   = as.data.frame(WA.tst.extr))

# get the values only
WA.predicted <- WA.tst.extr
WA.observed <- WA.spdf.tst$AADT
WA.residuals <- WA.observed - WA.predicted

# residuals as a spdf
WA.spdf.res <- SpatialPointsDataFrame(coords = as.data.frame(WA.sp.tst@coords),
                                      data   = as.data.frame(WA.residuals))

# output data

# variogram stuff
WAplot <- plot(WA.variogram)
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
rmse(WA.error)
mae(WA.error)
adj.r.sq <- summary(WA.lm)$adj.r.squared
plot(WA.observed~WA.predicted)

# maps

WA.layers <- writeRasters(WA.krige, WA.raster, "WA")

min(WA.krige$var1.pred)
max(WA.krige$var1.pred)
mean(WA.krige$var1.pred)
median(WA.krige$var1.pred)
sd(WA.krige$var1.pred)



