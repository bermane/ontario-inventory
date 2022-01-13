# This function aggregates polygons from mean shift algorithm with an 
# orginial NA value into a single large polygon

vect_in <- 'D:/ms_10_10_100.shp'

agg_na <- function(
  vect_in, # mean shift generated forest stand shape file path
  vect_out = stringr::str_replace(vect_in, '.shp', '_agg_na.shp') # output shape file with aggregated NA polygons
){
  
  # load polygon dataset
  p <- terra::vect(vect_in)
  
  # subset by polygons that only have one pixel (NA) and polygons that have more
  p_na <- p[p$nbPixels==1,]
  p_real <- p[p$nbPixels>1,]
  
  # dissolve polygons that only have 1 pixels
  p2 <- terra::aggregate(p_na, by='nbPixels')
  
  # add back into single file
  p3 <- rbind(p_real, p2)
  
  # write to file
  terra::writeVector(p3, vect_out)
}