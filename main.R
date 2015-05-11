# main.R

# call functions and libraries
source("/home/sean/Documents/trafficVolume/functions.R")

# define projection, albers equal area
aea <- "+proj=aea +lat_1=38 +lat_2=41 +lat_0=34 +lon_0=-114 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs"

# define directory with input AADT point data
data.directory <- "/home/sean/traffic/AADT/"
data.suffix <- "_AADT"

# polygon for the continental US
## need to have selectState function have cUS called from elsewhere....?
# IL didn't have enough memory on the 4gb computer
stateAbbreviations <- list("NH", "NJ", "NM", "NY", "ND")

# NC 
# Error in getErrorValues(validation.spdfs = tr.validation.dfs, column.names = column.names) : 
# unused argument (validation.spdfs = tr.validation.dfs) 

main("NH")

lapply(stateAbbreviations, main)

########################
# setwd("/home/sean/traffic/forExample/)

az <- readRDS("AZ_output.rds")
ca <- readRDS("CA_output.rds")
or <- readRDS("OR_output.rds")
wa <- readRDS("WA_output.rds")
nv <- readRDS("NV_output.rds")



