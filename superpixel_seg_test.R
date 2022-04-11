library(OpenImageR)
library(terra)
library(scales)
library(tidyverse)
library(magrittr)
library(SuperpixelImageSegmentation)

# create function to rescale values from 0 to 255
scale_255 <- function(x){
  
  # calculate 1st and 99th percentile of input raster
  x <- values(x, mat=F)
  
  # rescale raster values
  x <- scales::rescale(x, to = c(0,255))
  
  return(x)
}

# load lidar raster for segmentation
spl <- rast('D:/ontario_inventory/segmentation/spl_stack_mask_wat_ucl_polytype.tif')

# rescale rasters from 0 to 100
values(spl[[1]]) <- scale_255(spl[[1]])
values(spl[[2]]) <- scale_255(spl[[2]])
values(spl[[3]]) <- scale_255(spl[[3]])

# crop to test area
area <- vect('D:/ontario_inventory/romeo/RMF_Sample.shp') %>%
  project(., spl)

spl %<>% terra::crop(., area)

# write spl to disk
writeRaster(spl, 'D:/ontario_inventory/superpixel_test/spl.tif')

#turn spl into image
im <- readImage('D:/ontario_inventory/superpixel_test/spl.tif')

# convert to superpixel
res_slico <- superpixels(input_image = im,
                         method = 'slico',
                         superpixel = 40000,
                         return_slic_data = T,
                         return_labels = T,
                         write_slic = '',
                         verbose = T)

plot_slico = OpenImageR::NormalizeObject(res_slico$slic_data)
plot_slico = grDevices::as.raster(plot_slico)
graphics::plot(plot_slico)

# try superpixel segmentation
im <- OpenImageR::readImage('D:/ontario_inventory/superpixel_test/spl.tif')

# run segmentation
init <- Image_Segmentation$new()
spx <- init$spixel_segmentation(input_image = im, 
                                superpixel = 40000, 
                                AP_data = TRUE,
                                use_median = TRUE, 
                                sim_wL = 1, 
                                sim_wA = 1, 
                                sim_wB = 1,
                                sim_color_radius = 10, 
                                verbose = TRUE)

str(spx)

spl2 <- spl
values(spl2) <- spx$AP_image_data

spl2

plot(spl2[[1]])
plot(spl2[[2]])
