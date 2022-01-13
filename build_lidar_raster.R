# This function builds a 3-band LiDAR raster, which will be input into
# the mean shift image segmentation algorithm

lidar_bands <- c('D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_lor.tif',
                 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif',
                 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_cv.tif')
window <- 5
mask <- c('D:/ontario_inventory/romeo/RMF_EFI_layers/Roads/RMF_roads_hpb_mask.shp',
          'D:/ontario_inventory/romeo/RMF_EFI_layers/Lakes and rivers/RMF_waterbodies.shp')
filename <- 'D:/ontario_inventory/segmentation/spl_stack.tif'
overwrite = T

build_lidar_raster <- function(
  bands, # 3 part character vector with paths to raster bands. must be same projection
  window = 5, # smoothing window for lidar bands, defaults to 5
  mask, # character vector with paths to vector files all features are masked
  filename, # path to write out 3-band raster from lidar
  overwrite = F # overwrite output raster, defaults to F
                           ){
  
  # stack lidar bands
  spl <- terra::rast(lidar_bands)
  
  # apply smoothing function on 5 cell square
  spl[[1]] <- terra::focal(spl[[1]], w=5, fun="mean")
  spl[[2]] <- terra::focal(spl[[2]], w=5, fun="mean")
  spl[[3]] <- terra::focal(spl[[3]], w=5, fun="mean")
  
  # load mask layers and project to match spl
  mask_lyrs <- lapply(mask, FUN = function(x) {
    x <- terra::vect(x)
    terra::project(x, spl)
    })

  # apply mask to lidar bands
  for(i in 1:length(mask_lyrs)){
    spl <- terra::mask(spl, mask_lyrs[[i]], inverse = T)
  }
  
  # if any band is missing values set all to NA
  spl[is.na(spl[[1]])] <- NA
  spl[is.na(spl[[2]])] <- NA
  spl[is.na(spl[[3]])] <- NA
  
  # create function to rescale values from 0 to 100 using 1 and 99 percentile
  scale_100 <- function(x){
    
    # calculate 1st and 99th percentile of input raster
    perc <- terra::values(x, mat=F)
    perc <- quantile(perc, probs=c(0.01, 0.99), na.rm=T)
    
    # rescale raster using 1st and 99th %
    x <- (x-perc[1])/(perc[2] - perc[1]) * 100
    
    #reset values below 0 and above 100
    x[x < 0] <- 0
    x[x > 100] <- 100
    
    return(x)
  }
  
  # rescale rasters from 0 to 100
  spl[[1]] <- scale_100(spl[[1]])
  spl[[2]] <- scale_100(spl[[2]])
  spl[[3]] <- scale_100(spl[[3]])
  
  # write raster to tif
  terra::writeRaster(spl, filename = filename, overwrite = overwrite)
  
}



