# this code extracts lidar attributes within FRI polygons
# it uses the median value taken as the fraction of 
# pixel covered by each polygon

# load packages
library(terra)
library(tidyverse)
library(exactextractr)
library(sf)

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

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
vec <- exact_extract(lidar_ras, poly_ras, 'median')

# change column names
colnames(vec) <- names(lidar)

# add new column into dat
dat <- cbind(dat, vec)

################################
### ADD SENTINAL REFLECTANCE ###
################################

# load sentinel BOA mosaic
sent_ras <- rast('D:/ontario_inventory/romeo/Sentinel/Mosaic/BOA/S2_BOA_20_Mosaic.tif')

# project poly to crs of raster
poly_ras <- project(poly, sent_ras)

# convert to sf
poly_ras <- st_as_sf(poly_ras)

#extract median values
vec <- exact_extract(sent_ras, poly_ras, 'median')

# change column names
colnames(vec) <- names(sent_ras)

# add new columns into dat
dat <- cbind(dat, vec)

####################################
### ADD ADDITIONAL LIDAR METRICS ###
####################################

# load PZABOVE data
ras <- rast('D:/ontario_inventory/romeo/SPL metrics/PZABOVE_MOSAIC/RMF_PZABOVE_MOSAIC.tif')

# set band names
names(ras) <- c('PZABOVEMEAN', 'PZABOVE2', 'PZABOVE5')

# project poly to crs of raster
poly_ras <- project(poly, ras)

# convert to sf
poly_ras <- st_as_sf(poly_ras)

#extract median values
vec <- exact_extract(ras, poly_ras, 'median')

# change column names
colnames(vec) <- names(ras)

# add new columns into dat
dat <- cbind(dat, vec)

# I already extracted the rumple and depth_q25 into the individual metrics
# folder so we can skip these two and move directly to individual

# # load RUMPLE data
# ras <- rast('D:/ontario_inventory/romeo/SPL metrics/RUMPLE_MOSAIC/RMF_RUMPLE_MOSAIC.tif')
# 
# # set band names
# names(ras) <- 'RUMPLE'
# 
# # project poly to crs of raster
# poly_ras <- project(poly, ras)
# 
# # convert to sf
# poly_ras <- st_as_sf(poly_ras)
# 
# #extract median values
# vec <- exact_extract(ras, poly_ras, 'median')
# 
# # add new columns into dat
# dat <- cbind(dat, vec)
# 
# # load SAD data
# ras <- rast('D:/ontario_inventory/romeo/SPL metrics/SAD_MOSAIC/RMF_SAD_MOSAIC.tif')
# 
# # set band names
# names(ras) <- c("depth_mean", "depth_max", 
#                 "depth_min", "depth_sd", 
#                 "depth_q10", "depth_q25", 
#                 "depth_q50", "depth_q75", 
#                 "depth_q95", "depth_q99")
# 
# # project poly to crs of raster
# poly_ras <- project(poly, ras)
# 
# # convert to sf
# poly_ras <- st_as_sf(poly_ras)
# 
# #extract median values
# vec <- exact_extract(ras, poly_ras, 'median')
# 
# # add new columns into dat
# dat <- cbind(dat, vec)

# load additional Z metrics data
ras_files <- list.files('D:/ontario_inventory/romeo/SPL metrics/Z_METRICS_MOSAIC/individual',
                        full.names = T)
ras_names <- list.files('D:/ontario_inventory/romeo/SPL metrics/Z_METRICS_MOSAIC/individual')
ras_names <- sapply(ras_names, function(x){
  str_split_fixed(x, pattern = '_', n = 5) %>% .[5] %>%
    str_replace(pattern = '.tif', replacement = '')
})

# load rasters
ras <- rast(ras_files)

# set band names
names(ras) <- ras_names

# project poly to crs of raster
poly_ras <- project(poly, ras)

# convert to sf
poly_ras <- st_as_sf(poly_ras)

#extract median values
vec <- exact_extract(ras, poly_ras, 'median')

# change column names
colnames(vec) <- names(ras)

# add new columns into dat
dat <- cbind(dat, vec)

# save extracted dataframe
save(dat, file = 'D:/ontario_inventory/dat/dat_fri_extr.RData')


