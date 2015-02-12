#
#
#
#
#
#
#
#
#
#



getProj <- function(vct_name) {
  
  #############
  
  # INPUTS
  # vct_name  = vector name, including file ending
  
  # OUTPUTS
  # proj4string = proj4string
  
  #############
  
  # run system command
  # outputs to temp txt file
  system(paste("gdalsrsinfo", vct_name, "> proj4str.txt", sep = " "))
  
  # get second line of txt file
  proj4string <- readLines("proj4str.txt", n = 2)
  
  # save second line as variable
  proj4string <- substr(proj4string[2], 10, 200)
  
  return(proj4string)
}

createSpdf <- function(vct_name, vct_dir, proj, col_names) {
  
  #############
  
  # INPUTS
  # vct_name  = vector name, including file ending
  # vct_dir   = vector directory
  # proj      = proj4string
  # col_names = column names, form of a vector! (an ice dragon!)
  
  # OUTPUTS
  # spdf      = spatial points data frame that contains desired column/s
  
  #############
  
  # bring vector file into memory
  vector <- readOGR(dsn = vct_dir, layer = vct_name) #, p4s = NAD27)
  
  # get proj4string
  
  # run system command to pull .prj info
  
  # access as string

  # convert targeted column into df
  vector_df <- as.data.frame(vector$col_names)
  # get shpfile coordinates pair
  vector_coords <- coordinates(vector)[,1:2]
  
  vector_spdf <- SpatialPointsDataFrame(vector_coords, vector_df)
  # set spdf projection, column names
  proj4string(vector_spdf) <- proj
  names(vector_spdf) <- col_names
  
}