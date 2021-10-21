# This code applies image segmentation to SPL data over the Romeo Malette Forest in Ontario
# in order to create forest stand polygons

# load packages
library(terra)
library(meanShiftR)
library(tidyverse)

# set names of SPL rasters to stack
lor <- 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_lor.tif'
cc <- 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif'
cv <- 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_cv.tif'

# stack rasters
spl <- rast(c(lor, cc, cv))

# apply smoothing function on 5 cell square
spl[[1]] <- focal(spl[[1]], w=5, fun="mean")
spl[[2]] <- focal(spl[[2]], w=5, fun="mean")
spl[[3]] <- focal(spl[[3]], w=5, fun="mean")

# if any band is missing values set all to NA
spl[is.na(spl[[1]])] <- NA
spl[is.na(spl[[2]])] <- NA
spl[is.na(spl[[3]])] <- NA

# scale values
spl_scale <- scale(spl, center=F)

# set no data values to 100
spl_scale[is.na(spl_scale)] <- 100
spl_scale[is.nan(spl_scale)] <- 100

# create function to rescale values from 0 to 100 using 1 and 99 percentile
scale_100 <- function(x){
  
  # calculate 1st and 99th percentile of input raster
  perc <- values(x, mat=F) %>% quantile(., probs=c(0.01, 0.99), na.rm=T)
  
  # rescale raster using 1st and 99th %
  x <- (x-perc[1])/(perc[2] - perc[1]) * 100
  
  #reset values below 0 and above 100
  x[x < 0] <- 0
  x[x > 100] <- 100
  
  return(x)
}

# rescale rasters from 0 to 100
spl_scale_100 <- spl
spl_scale_100[[1]] <- scale_100(spl_scale_100[[1]])
spl_scale_100[[2]] <- scale_100(spl_scale_100[[2]])
spl_scale_100[[3]] <- scale_100(spl_scale_100[[3]])

# write raster to tif
writeRaster(spl_scale, filename='D:/ontario_inventory/test/spl_scale.tif', overwrite=T)
writeRaster(spl_scale_100, filename='D:/ontario_inventory/test/spl_scale_100.tif', overwrite=T)
writeRaster(spl, filename='D:/ontario_inventory/test/spl.tif', overwrite=T)

# create function to extract polygon stats from qgis meanshift test runs
extract_ms_stats <- function(file, name){
  
  # load polygon and remove polygons with 1 pixel (eventually need to merge these maybe?)
  # need a spatially contiguous file
  p <- vect(file) %>% .[.$nbPixels>1,]
  
  # create dataframe with values wanted
  df <- data.frame(name = name,
                   min_pix = min(p$nbPixels),
                   mean_pix = mean(p$nbPixels) %>% round(1),
                   med_pix = median(p$nbPixels) %>% round(1),
                   max_pix = max(p$nbPixels))
  
  # plot histogram
  png(str_c('D:/ontario_inventory/test/plots/', name, '.png'))
  hist(p$nbPixels, 
       xlim=quantile(pd$nbPixels, probs = c(0.01, 0.99)), 
       breaks = 1000,
       main = name,
       xlab = 'Number of Pixels')
  dev.off()
  
  #return df
  return(df)
}

# compare meanshift test run specs
ms_df <- extract_ms_stats('D:/ontario_inventory/test/scale_100/scale_100_ms_10_75_50.shp',
                          'scale_100_ms_10_7-5_50')
ms_df <- rbind(ms_df,
               extract_ms_stats('D:/ontario_inventory/test/scale_100/scale_100_ms_10_5_100.shp',
                                'scale_100_ms_10_5_100'))

# INCLUDE TOTAL NUMBER OF POLYGONS???

# load interpreter derived polygons to extract statistics
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp') %>%
  as.data.frame

# create df of interpreter polygons
# from full interpreter dataset
int_df <- data.frame(name = 'interp_full',
           min_pix = (min(poly$AREA)/400) %>% round(1),
           mean_pix = (mean(poly$AREA)/400) %>% round(1),
           med_pix = (median(poly$AREA)/400) %>% round(1),
           max_pix = (max(poly$AREA)/400) %>% round(1))

# subset forested polygons only
poly_for <- poly[poly$POLYTYPE=='FOR',]

# from forested polygons only
int_df <- rbind(int_df,
                data.frame(name = 'interp_for',
                           min_pix = (min(poly_for$AREA)/400) %>% round(1),
                           mean_pix = (mean(poly_for$AREA)/400) %>% round(1),
                           med_pix = (median(poly_for$AREA)/400) %>% round(1),
                           max_pix = (max(poly_for$AREA)/400) %>% round(1)))

# subset all non water polygons
poly_land <- poly[poly$POLYTYPE!='WAT',]

# from non water polygons only
int_df <- rbind(int_df,
                data.frame(name = 'interp_land',
                           min_pix = (min(poly_land$AREA)/400) %>% round(1),
                           mean_pix = (mean(poly_land$AREA)/400) %>% round(1),
                           med_pix = (median(poly_land$AREA)/400) %>% round(1),
                           max_pix = (max(poly_land$AREA)/400) %>% round(1)))

# output histograms from interpreter polygons
# area hist all poly types
hist(poly$AREA/400, 
     xlim=(quantile(poly$AREA, probs = c(0.01, 0.99)))/400,
     ylim = c(0,15000),
     breaks=10000,
     main = 'All Interpreter Polygons',
     xlab = 'Number of Pixels')

# area hist only forest polys
hist(poly_for$AREA/400, 
     xlim=(quantile(poly_for$AREA, probs = c(0.01, 0.99)))/400,
     ylim = c(0,15000),
     breaks=10000,
     main = 'Forested Interpreter Polygons',
     xlab = 'Number of Pixels')

# area hist non water polys
hist(poly_land$AREA/400, 
     xlim=(quantile(poly_land$AREA, probs = c(0.01, 0.99)))/400,
     ylim = c(0,15000),
     breaks=10000,
     main = 'Non-Water Interpreter Polygons',
     xlab = 'Number of Pixels')

# meanshift_otb <- function(otb.path = "", raster.in = "", out.path = "", name ="", filter.meanshift.spatialr = "5",
#                            filter.meanshift.ranger = "0.1", filter.meanshift.thres = "0.1",
#                            filter.meanshift.maxiter = "100", filter.meanshift.minsize = "100",
#                            mode.vector.outmode = "ovw", mode.vector.inmask = "", mode.vector.neighbor = "false",
#                            mode.vector.stitch = "true", mode.vector.minsize = 1, mode.vector.simplify = 0.1,
#                            mode.vector.layername = "layer", mode.vector.fieldname = "DN", mode.vector.tilesize = 1024,
#                            mode.vector.startlabel = 1){
#   # Set configuration      
#   conf <- paste("-in",raster.in,"-filter meanshift","-filter.meanshift.spatialr",filter.meanshift.spatialr,
#                 "-filter.meanshift.ranger",filter.meanshift.ranger,"-filter.meanshift.thres",filter.meanshift.thres,
#                 "-filter.meanshift.maxiter",filter.meanshift.maxiter,"-filter.meanshift.minsize",filter.meanshift.minsize,
#                 "-mode vector","-mode.vector.out",paste(out.path,"/",name,".shp",sep=""),"-mode.vector.outmode",mode.vector.outmode,
#                 ifelse(missingArg(mode.vector.inmask),"",paste("-mode.vector.inmask",mode.vector.inmask)),
#                 "-mode.vector.neighbor", mode.vector.neighbor,
#                 "-mode.vector.stitch",mode.vector.stitch,
#                 "-mode.vector.minsize",mode.vector.minsize,
#                 "-mode.vector.simplify",mode.vector.simplify,
#                 "-mode.vector.layername",mode.vector.layername,
#                 "-mode.vector.fieldname",mode.vector.fieldname,
#                 "-mode.vector.tilesize",mode.vector.tilesize,
#                 "-mode.vector.startlabel",mode.vector.startlabel)
#   # apply function in command line
#   system(paste(otb.path,"/otbcli_Segmentation"," ",conf,sep=""))
#   # save configuration for further use
#   write.table(x = conf,file = paste(out.path,"/",name,"_conf.txt",sep=""),row.names = F, col.names = F)
# }
# 
# # usage, you can set any option listed above
# meanshift_otb(otb.path="C:/OTB/bin", 
#               raster.in='D:/ontario_inventory/test/spl.tif', 
#               out.path="D:/ontario_inventory/test/", 
#               name="ms_test")
# 
# # or, to read the output into R
# out_path = '/out/path'
# lyr_name = 'test'
# 
# # usage
# meanshift.segm(otb.path = "OTB-5.8.0-Darwin64/bin", raster.in = "path to/rater_in.tif", out.path=out_path, name = lyr_name)
# 
# shp <- readOGR(dsn=out_path, layer = lyr_name, driver = 'ESRI Shapefile')


# # set missing values to 0 just to try
# spl[is.na(spl)] <- 0
# spl[is.nan(spl)] <- 0
# 
# # create matrix of values
# spl_mat <- t(as.matrix(spl))
# 
# # try to apply mean shift
# ms <- meanShift(t(as.matrix(spl)))
# 
# # create as raster
# ms_ras <- terra::setValues(spl[[1]], as.numeric(ms$value[1,]))
