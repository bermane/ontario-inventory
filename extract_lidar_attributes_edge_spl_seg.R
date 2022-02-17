# this code extracts lidar attributes within spl segmented polygons
# it only uses data if more than 95% of a given pixel is inside the polygon
# the idea being edge pixels may negatively impact performance

# load packages
library(terra)
library(tidyverse)
library(exactextractr)
library(sf)

# load spl segmented polygons
poly <- vect('D:/ontario_inventory/segmentation/ms_10_10_100_for_only_agg_na.shp')

# convert to df
dat <- as.data.frame(poly)

# load LiDAR datasets we need to extract over polygons
# create named vector with variable names and data links

lidar <- c('cc' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif',
           'avg' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_avg.tif',
           'max' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_max.tif',
           'p95' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_p95.tif',
           'qav' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_qav.tif',
           'ske' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_ske.tif',
           'kur' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_kur.tif',
           'cv' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_cv.tif',
           'lor' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_lor.tif',
           'ba' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_ba_ha.tif',
           'qmdbh' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_qmdbh.tif',
           'dens' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_dens.tif',
           'agb' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_AGB_ha.tif',
           'top_height' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_top_height.tif',
           'v' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_V_ha.tif',
           'v_merch' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_Vmerch_ha.tif')

# load LiDAR rasters as raster stack
lidar_ras <- rast(lidar)

# project poly to crs of raster
poly_ras <- project(poly, lidar_ras)
  
# convert to sf
poly_ras <- st_as_sf(poly_ras)

#extract median values
vec <- exact_extract(lidar_ras, poly_ras, function(values, coverage_fraction){
  values <- values[coverage_fraction == 1,]
  apply(values, 2, function(x) median(x, na.rm = T))
})

# transpose matrix and make data.frame
vec <- t(vec) %>% as.data.frame

# change column names
colnames(vec) <- names(lidar)

# add new column into dat
dat <- cbind(dat, vec)

dat_lidar <- dat

# save extracted dataframe for fast rebooting
save(dat_lidar, file = 'D:/ontario_inventory/imputation/seg_df_ms_10_10_100_agg_na.RData')
