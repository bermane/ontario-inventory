# load packages
library(terra)
library(tidyverse)

# load shapefile used to crop all layers
mask <- vect('D:/Image Segmentation Student Project/boundary file/RMF_Sample_student.shp')

# load files from lakes and rivers
name <- list.files(path = 'D:/ontario_inventory/romeo/RMF_EFI_layers/Lakes and rivers',
                   pattern = glob2rx('*.shp'))
loc <- list.files(path = 'D:/ontario_inventory/romeo/RMF_EFI_layers/Lakes and rivers',
                  pattern = glob2rx('*.shp'), full.names = T)

# loop through files, crop and put in output dir
for(i in 1:length(name)){
  
  # load vector
  vec <- vect(loc[i])
  
  # project mask to match input file
  mask1 <- project(mask, vec)
  
  # crop input file and write
  out <- crop(vec, mask1)
  
  # write to disk
  writeVector(out, filename = str_c('D:/Image Segmentation Student Project/Lakes and rivers/clip_',
                                    name[i]))
}

# load files from roads
name <- list.files(path = 'D:/ontario_inventory/romeo/RMF_EFI_layers/Roads',
                   pattern = glob2rx('*.shp'))
loc <- list.files(path = 'D:/ontario_inventory/romeo/RMF_EFI_layers/Roads',
                  pattern = glob2rx('*.shp'), full.names = T)

# loop through files, crop and put in output dir
for(i in 1:length(name)){
  
  # load vector
  vec <- vect(loc[i])
  
  # project mask to match input file
  mask1 <- project(mask, vec)
  
  # crop input file and write
  out <- crop(vec, mask1)
  
  # write to disk
  writeVector(out, filename = str_c('D:/Image Segmentation Student Project/Roads/clip_',
                                    name[i]))
}

# load files from Polygons Inventory
name <- list.files(path = 'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory',
                   pattern = glob2rx('*.shp'))
loc <- list.files(path = 'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory',
                  pattern = glob2rx('*.shp'), full.names = T)

# loop through files, crop and put in output dir
for(i in 1:length(name)){
  
  # load vector
  vec <- vect(loc[i])
  
  # project mask to match input file
  mask1 <- project(mask, vec)
  
  # crop input file and write
  out <- crop(vec, mask1)
  
  # write to disk
  writeVector(out, filename = str_c('D:/Image Segmentation Student Project/Polygons Inventory/clip_',
                                    name[i]))
}

# load files from Fire extent 2012
name <- list.files(path = 'D:/ontario_inventory/romeo/RMF_EFI_layers/Fire extent 2012',
                   pattern = glob2rx('*.shp'))
loc <- list.files(path = 'D:/ontario_inventory/romeo/RMF_EFI_layers/Fire extent 2012',
                  pattern = glob2rx('*.shp'), full.names = T)

# loop through files, crop and put in output dir
for(i in 1:length(name)){
  
  # load vector
  vec <- vect(loc[i])
  
  # project mask to match input file
  mask1 <- project(mask, vec)
  
  # crop input file and write
  out <- crop(vec, mask1)
  
  # write to disk
  writeVector(out, filename = str_c('D:/Image Segmentation Student Project/Fire extent 2012/clip_',
                                    name[i]))
}

# load files from ABA layers SPL 2018
name <- list.files(path = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018',
                   pattern = glob2rx('*FOR.tif'))
loc <- list.files(path = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018',
                  pattern = glob2rx('*FOR.tif'), full.names = T)

# loop through files, crop and put in output dir
for(i in 1:length(name)){
  
  # load raster
  ras <- rast(loc[i])
  
  # project mask to match input file
  mask1 <- project(mask, ras)
  
  # crop input file and write
  out <- crop(ras, mask1)
  
  # write to disk
  writeRaster(out, filename = str_c('D:/Image Segmentation Student Project/ABA layers SPL 2018/clip_',
                                    name[i]))
}

# load files from SPL100 metrics
name <- list.files(path = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics',
                   pattern = glob2rx('*.tif'))
loc <- list.files(path = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics',
                  pattern = glob2rx('*.tif'), full.names = T)

# loop through files, crop and put in output dir
for(i in 1:length(name)){
  
  # load raster
  ras <- rast(loc[i])
  
  # project mask to match input file
  mask1 <- project(mask, ras)
  
  # crop input file and write
  out <- crop(ras, mask1)
  
  # write to disk
  writeRaster(out, filename = str_c('D:/Image Segmentation Student Project/SPL100 metrics/clip_',
                                    name[i]))
}


