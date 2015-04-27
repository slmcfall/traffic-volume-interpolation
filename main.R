# main.R

# call functions and libraries
# source("/Users/seanmcfall/Documents/traffic-volume-interpolation/functions.R")

# define projection, albers equal area
aea <- "+proj=aea +lat_1=38 +lat_2=41 +lat_0=34 +lon_0=-114 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs"

# define directory with input AADT point data
data.directory <- "/home/sean/Documents/traffic-volume-interpolation/AADT/"
data.suffix <- "_AADT"

# polygon for the continental US
## need to have selectState function have cUS called from elsewhere....?


stateAbbreviations <- list("WA", "OR")

lapply(stateAbbreviations, main)

########################

az <- readRDS("AZ_output.rds")
ca <- readRDS("CA_output.rds")
or <- readRDS("OR_output.rds")
wa <- readRDS("WA_output.rds")
nv <- readRDS("NV_output.rds")

annotatedplot(krigeobj = wa[[4]], title = "Washington Variogram")
annotatedplot(krigeobj = or[[4]], title = "Oregon Variogram")
annotatedplot(krigeobj = nv[[4]], title = "Nevada Variogram")

### box plot log data

# bring in the raw AADT data as spdf
az.spdf <- createSpdf(vct_path = "/Users/seanmcfall/Documents/traffic-volume-interpolation/AADT/AZ_AADT.shp", col_names = c("AADT"), proj = aea)
#....stupid R studio won't copypasta....

# log transform spdf's
az.log <- logTransformSpdf(spdf = az.spdf, column.names = c("AADT"), projection = aea)

# display as boxplots
box.list <- list(az.log$AADT, az.spdf$AADT)

boxplot(box.list, outline = FALSE, main = "Log Transformed AADT by State")
mtext("STATE", side = 1, line = 3, font = 4)
mtext("AADT", side = 2, line = 3, font = 4)

### box plot observed and expected

getObsList <- function(rds) {
  
  # capture the data from the lists
  rds.obs <- rds[2][[1]]$Observed
  
  # get them as values
  rds.obs.vals <- rds.obs$AADT
  
  return (rds.obs.vals)
}

getPredList <- function(rds) {
  
  # capture data from rds list
  rds.pred <- rds[2][[1]]$Predicted
  
  # get as values
  rds.pred.vals <- rds.pred$AADT
  
  return(rds.pred.vals)
}

box.names = c("        Arizona", "", 
              "        California", "",
              "        Nevada", "",
              "        Oregon", "",
              "        Washington", "")

boxplot(list(getObsList(az), getPredList(az), getObsList(ca), getPredList(ca), getObsList(nv), getPredList(nv), getObsList(or), getPredList(or), getObsList(wa), getPredList(wa)), 
        names = box.names, at =c(1,2,4,5,7,8,10,11,13,14), outline = FALSE, col = c("orange", "blue"),
        par(mar = c(5, 5, 2, 15)+ 0.1),
        xlab = 'STATE', ylab = 'log(AADT)', main = 'Observed vs. Predicted Distribution ')

a <- getObsExpList(az)


foo[[1]]

# box plot 
boxplot(az.list)






