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

createSpdf <- function(vct_name, vct_dir, proj, col_names) {
  
  #############
  
  # INPUTS
  # vct_name  = vector name, including file ending
  # vct_dir   = vector directory
  # proj      = proj4string
  # col_names = column name
  
  # OUTPUTS
  # spdf      = spatial points data frame that contains desired column/s
  
  #############
  
  # bring vector file into memory
  vector <- readOGR(dsn = vct_dir, layer = vct_name) #, p4s = NAD27)

  # convert targeted column into df
  vector_df <- as.data.frame(vector$col_names)
  # get shpfile coordinates pair
  vector_coords <- coordinates(vector)[,1:2]
  
  vector_spdf <- SpatialPointsDataFrame(vector_coords, vector_df)
  # set spdf projection, column names
  proj4string(vector_spdf) <- proj
  names(vector_spdf) <- col_names
  
  
}