# This function loads and extracts raster attributes over polygons

polygons <- 'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp'
vars <- c('cc' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif',
           'avg' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_avg.tif',
           'max' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_max.tif',
           'qav' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_qav.tif',
           'ske' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_ske.tif',
           'kur' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_kur.tif',
           'cv' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_cv.tif',
           'lor' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_lor.tif',
           'ba' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_ba_ha.tif',
           'qmdbh' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_qmdbh.tif',
           'dens' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_dens.tif',
           'agb' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_AGB_ha.tif')
out_path <- 'D:/ontario_inventory/imputation/test_polygon_attributes/fri_lidar_vars.shp'
overwrite <- T

build_polygon_attributes <- function(
  polygons, # shapefile of polygon dataset over which to load attributes
  vars, # a named character vector of variables to extract and the path to the raster files (could make this a multiband?)
  out_path, # shapefile to save polygon dataset with loaded attributes 
  overwrite = F # overwrite output shapefile
){
  
  # load polygons
  poly <- terra::vect(polygons)
  
  # convert polygons to data frame
  dat <- terra::as.data.frame(poly)
  
  # loop through LiDAR datasets and add to main dataframe
  for(i in 1:length(vars)){
    
    # load raster variable
    ras <- terra::rast(vars[i])
    
    # project poly to crs of raster
    poly_ras <- terra::project(poly, ras)
    
    # extract median values within each polygon
    ras_med <- terra::extract(ras, poly_ras, fun = function(x){median(x, na.rm = T)})
    
    # add new column into dat
    dat <- tibble::add_column(dat, ras = ras_med[,2])
    
    # change column name
    colnames(dat)[NCOL(dat)] <- names(vars[i])
    
  }
  
  # repopulate shape file with new data
  terra::values(poly) <- dat
  
  # write polygons with new attributes
  terra::writeVector(poly, filename = out_path, overwrite = overwrite)
  
}

