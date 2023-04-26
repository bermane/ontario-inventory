# this code extracts lidar and aux attributes within FRI polygons
# it uses the median value taken as the fraction of 
# pixel covered by each polygon

# load packages
library(terra)
library(tidyverse)
library(exactextractr)
library(sf)

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/FSF/FRI/FSF_opi_polygon_CSRS_NAD83_17.shp')

# convert to df
dat <- as.data.frame(poly)

# load LiDAR datasets we need to extract over polygons
# create named vector with variable names and data links

# all layers must have same extent

lidar <- c('cc' = 'D:/ontario_inventory/FSF/ALS/cov_2m.img',
           'avg' = 'D:/ontario_inventory/FSF/ALS/zmean.img',
           'max' = 'D:/ontario_inventory/FSF/ALS/zmax.img',
           'p95' = 'D:/ontario_inventory/FSF/ALS/zq95.img',
           'ske' = 'D:/ontario_inventory/FSF/ALS/zskew.img',
           'sd' = 'D:/ontario_inventory/FSF/ALS/zsd.img',
           'kur' = 'D:/ontario_inventory/FSF/ALS/zkurt.img',
           'cv' = 'D:/ontario_inventory/FSF/ALS/cv.img',
           'rumple' = 'D:/ontario_inventory/FSF/ALS/rumpleIndex.img',
           'zentropy' = 'D:/ontario_inventory/FSF/ALS/zentropy.img', 
           'zpcum8' = 'D:/ontario_inventory/FSF/ALS/zpcum8.img',
           'lor' = 'D:/ontario_inventory/FSF/EFI/LoreyHeight.img',
           'ba' = 'D:/ontario_inventory/FSF/EFI/BasalArea.img',
           'qmdbh' = 'D:/ontario_inventory/FSF/EFI/QMD.img',
           'stems' = 'D:/ontario_inventory/FSF/EFI/Stems.img',
           'agb' = 'D:/ontario_inventory/FSF/EFI/Biomass.img',
           'top_height' = 'D:/ontario_inventory/FSF/EFI/TopHt.img',
           'gtv' = 'D:/ontario_inventory/FSF/EFI/GTV.img',
           'gmvwl' = 'D:/ontario_inventory/FSF/EFI/GMV_WL.img')

# load LiDAR rasters as raster stack
lidar_ras <- rast(lidar)

# change layer names
names(lidar_ras) <- names(lidar)

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

###################################################
### ADD SENTINEL RED EDGE 2 SURFACE REFLECTANCE ###
###################################################

s2 <- c('b6' = 'D:/ontario_inventory/FSF/sentinel/boa/b6_fsf_2018.tif')

# load as raster
s2_ras <- rast(s2)

# change layer names
names(s2_ras) <- names(s2)

# project poly to crs of raster
poly_ras <- project(poly, s2_ras)

# convert to sf
poly_ras <- st_as_sf(poly_ras)

#extract median values
vec <- exact_extract(s2_ras, poly_ras, 'median')

# add new column into dat
dat %<>% add_column(b6 = vec)

###########################
### SAVE EXTRACTED DATA ###
###########################

# save extracted dataframe
save(dat, file = 'D:/ontario_inventory/dat/dat_fri_fsf.RData')
