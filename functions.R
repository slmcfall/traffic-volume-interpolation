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

getProj <- function(vct_path) {
  
  #############
  
  # INPUTS
  # vct_path  = vector name, including file ending
  
  # OUTPUTS
  # proj4string = proj4string
  
  #############
  
  # run system command
  # outputs to temp txt file
  system(paste("gdalsrsinfo", vct_path, "> proj4str.txt", sep = " "))
  
  # get second line of txt file
  proj4string <- readLines("proj4str.txt", n = 2)
  
  # save second line as variable
  proj4string <- substr(proj4string[2], 10, 200)
  
  # remove first and last quotes from the string
  proj4string <- substr(proj4string, 2, nchar(proj4string)-1)
  
  return(proj4string)
}

createSpdf <- function(vct_path, col_names) {
  
  #############
  
  # INPUTS
  # vct_path  = vector name, including file ending
  # vct_dir   = vector directory
  # proj      = proj4string
  # col_names = column names, form of...a vector!
  
  # OUTPUTS
  # vector_spdf      = spatial points data frame that contains desired column/s
  
  #############
  
  # get basename and directory of file
  vct_name <- basename(vct_path)
  vct_dir <- dirname(vct_path)
  
  # bring vector file into memory
  vector <- readOGR(dsn = vct_dir, layer = vct_name)
  
  # get proj4string
  proj <- getProj(paste(vct_path,".prj", sep = ""))

  # convert targeted column into df
  vector_df <- as.data.frame(vector)
  vector_df <- vector_df[col_names]
  
  # get shpfile coordinates pair
  vector_coords <- coordinates(vector)[,1:2]
  
  vector_spdf <- SpatialPointsDataFrame(vector_coords, vector_df)
  # set spdf projection, column names
  proj4string(vector_spdf) <- proj
  names(vector_spdf) <- col_names
  
  return (vector_spdf)
  
}

