# this code mosaics lidar metrics derived using LAStools

# load packages
library(terra)
library(tidyverse)
library(gdalUtilities)

##############################
###MOSAIC INDIVIDUAL BLOCKS###
##############################

# create list of lidar metric extensions
ext <- c('avg', 'cov', 'kur', 'max', 'p05', 'p10',
         'p20', 'p30', 'p40', 'p50', 'p60', 'p70',
         'p80', 'p90', 'p95', 'p99', 'qav', 'ske',
         'std')

# loop through files and mosiac
for(ex in ext){
  
  # find raster files
  ras <- list.files(path = 'E:/Block_2M/LAZ_Metrics',
                    pattern = str_c('*', ex, '.bil'), full.names = T)
  
  # build vrt
  vrt <- gdalbuildvrt(ras, output.vrt = str_c('D:/ontario_inventory/nipissing/als/BLOCK2M_Metrics/BLOCK2M_', ex, '.vrt'))
  
  r <- writeRaster(rast(vrt), filename = str_c('D:/ontario_inventory/nipissing/als/BLOCK2M_Metrics/BLOCK2M_', ex, '.tif'))
  
}

############################
###MOSAIC BLOCKS TOGETHER###
############################

# create list of lidar metric extensions
ext <- c('avg', 'cov', 'kur', 'max', 'p05', 'p10',
         'p20', 'p30', 'p40', 'p50', 'p60', 'p70',
         'p80', 'p90', 'p95', 'p99', 'qav', 'ske',
         'std')

# loop through files and mosiac
for(ex in ext){
  
  # find raster files
  ras <- list.files(path = 'D:/ontario_inventory/nipissing/als/BLOCK_Metrics',
                    pattern = str_c('*', ex, '.tif'), full.names = T, recursive = T)
  
  # build vrt
  vrt <- gdalbuildvrt(ras, output.vrt = str_c('D:/ontario_inventory/nipissing/als/NIPISSING_Metrics/NIPISSING_', ex, '.vrt'))
  
  r <- writeRaster(rast(vrt), filename = str_c('D:/ontario_inventory/nipissing/als/NIPISSING_Metrics/NIPISSING_', ex, '.tif'))
  
}

##################
###CALCULATE CV###
##################

sd <- rast('D:/ontario_inventory/nipissing/als/NIPISSING_Metrics/NIPISSING_std.tif')
avg <- rast('D:/ontario_inventory/nipissing/als/NIPISSING_Metrics/NIPISSING_avg.tif')

cv <- sd/avg

writeRaster(cv, filename = 'D:/ontario_inventory/nipissing/als/NIPISSING_Metrics/NIPISSING_cv.tif')


