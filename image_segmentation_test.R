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

# load roads, watercourses and waterbodies
roads <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Roads/RMF_roads.shp') %>%
  project(., spl)
waterc <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Lakes and rivers/RMF_watercourses.shp') %>%
  project(., spl)
waterb <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Lakes and rivers/RMF_waterbodies.shp') %>%
  project(., spl)

# subset roads to only mask the main types used by interpreter
# RDTYPE = H, P, B
roads <- roads[roads$RDTYPE %in% c('H', 'P', 'B'),]

# mask road and water body pixels to NA
spl <- spl %>% 
  mask(., roads, inverse = T) %>% 
  mask(., waterb, inverse = T)

#mask(., waterc, inverse = T) %>% 

# if any band is missing values set all to NA
spl[is.na(spl[[1]])] <- NA
spl[is.na(spl[[2]])] <- NA
spl[is.na(spl[[3]])] <- NA

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
spl[[1]] <- scale_100(spl[[1]])
spl[[2]] <- scale_100(spl[[2]])
spl[[3]] <- scale_100(spl[[3]])

# write raster to tif
writeRaster(spl, filename='D:/ontario_inventory/segmentation/spl_stack.tif', overwrite=T)

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
                   max_pix = max(p$nbPixels),
                   num_poly = NROW(p))
  
  # plot density
  ggplot(data.frame(nbPixels = p$nbPixels), aes(x = nbPixels)) +
    geom_density() +
    xlim(c(0,1500)) +
    ylim(c(0, 0.015)) +
    geom_vline(aes(xintercept = median(nbPixels, na.rm = T)), 
               linetype = "dashed", size = 0.6) +
    theme_bw() +
    xlab('Number of Pixels') +
    ylab('Density') +
    ggtitle(name)
  
  # save plot
  ggsave(str_c('D:/ontario_inventory/test/plots/', name, '.png'))
  
  #return df
  return(df)
}

# compare meanshift test run specs
ms_df <- extract_ms_stats('D:/ontario_inventory/test/scale_100/scale_100_ms_10_75_50.shp',
                          'scale_100_ms_10_7-5_50')
ms_df <- rbind(ms_df,
               extract_ms_stats('D:/ontario_inventory/test/scale_100/scale_100_ms_10_5_100.shp',
                                'scale_100_ms_10_5_100'))
ms_df <- rbind(ms_df,
               extract_ms_stats('D:/ontario_inventory/test/scale_100/scale_100_ms_10_10_100.shp',
                                'scale_100_ms_10_10_100'))
ms_df <- rbind(ms_df,
               extract_ms_stats('D:/ontario_inventory/test/scale_100/scale_100_ms_10_10_50.shp',
                                'scale_100_ms_10_10_50'))
ms_df <- rbind(ms_df,
               extract_ms_stats('D:/ontario_inventory/test/scale_100/scale_100_ms_10_15_50.shp',
                                'scale_100_ms_10_15_50'))
ms_df <- rbind(ms_df,
               extract_ms_stats('D:/ontario_inventory/test/scale_100/scale_m_100_ms_10_15_50.shp',
                                'scale_m_100_ms_10_15_50'))
ms_df <- rbind(ms_df,
               extract_ms_stats('D:/ontario_inventory/test/scale_100/scale_m_100_ms_10_10_100.shp',
                                'scale_m_100_ms_10_10_100'))

#ms_df <- rbind(ms_df,
#               extract_ms_stats('D:/ontario_inventory/test/scale_100/scale_100_ms_10_10_50.shp',
#                                'scale_100_ms_10_15_50'))

# load interpreter derived polygons to extract statistics
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp') %>%
  as.data.frame

# create df of interpreter polygons
# from full interpreter dataset
int_df <- data.frame(name = 'interp_full',
           min_pix = (min(poly$AREA)/400) %>% round(1),
           mean_pix = (mean(poly$AREA)/400) %>% round(1),
           med_pix = (median(poly$AREA)/400) %>% round(1),
           max_pix = (max(poly$AREA)/400) %>% round(1),
           num_poly = NROW(poly))

# subset forested polygons only
poly_for <- poly[poly$POLYTYPE=='FOR',]

# from forested polygons only
int_df <- rbind(int_df,
                data.frame(name = 'interp_for',
                           min_pix = (min(poly_for$AREA)/400) %>% round(1),
                           mean_pix = (mean(poly_for$AREA)/400) %>% round(1),
                           med_pix = (median(poly_for$AREA)/400) %>% round(1),
                           max_pix = (max(poly_for$AREA)/400) %>% round(1),
                           num_poly = NROW(poly_for)))

# subset all non water polygons
poly_land <- poly[poly$POLYTYPE!='WAT',]

# from non water polygons only
int_df <- rbind(int_df,
                data.frame(name = 'interp_land',
                           min_pix = (min(poly_land$AREA)/400) %>% round(1),
                           mean_pix = (mean(poly_land$AREA)/400) %>% round(1),
                           med_pix = (median(poly_land$AREA)/400) %>% round(1),
                           max_pix = (max(poly_land$AREA)/400) %>% round(1),
                           num_poly = NROW(poly_land)))

# combine dfs from mean shift and interpreter
out_df <- rbind(ms_df, int_df)

# write df as csv
write.csv(out_df, 
          file = 'D:/ontario_inventory/test/plots/polygon_table.csv',
          row.names = F)

# output histograms from interpreter polygons
# area hist all poly types
ggplot(poly, aes(x = AREA/400)) +
  geom_density() +
  xlim(c(0,1500)) +
  ylim(c(0, 0.015)) +
  geom_vline(aes(xintercept = median(AREA/400, na.rm = T)), 
             linetype = "dashed", size = 0.6) +
  theme_bw() +
  xlab('Number of Pixels') +
  ylab('Density') +
  ggtitle('All Interpreter Polygons')

# save plot
ggsave('D:/ontario_inventory/test/plots/all_interp.png')

# area hist only forest polys
ggplot(poly_for, aes(x = AREA/400)) +
  geom_density() +
  xlim(c(0,1500)) +
  ylim(c(0, 0.015)) +
  geom_vline(aes(xintercept = median(AREA/400, na.rm = T)), 
             linetype = "dashed", size = 0.6) +
  theme_bw() +
  xlab('Number of Pixels') +
  ylab('Density') +
  ggtitle('Forested Interpreter Polygons')

# save plot
ggsave('D:/ontario_inventory/test/plots/forested_interp.png')

# area hist non water polys
ggplot(poly_land, aes(x = AREA/400)) +
  geom_density() +
  xlim(c(0,1500)) +
  ylim(c(0, 0.015)) +
  geom_vline(aes(xintercept = median(AREA/400, na.rm = T)), 
             linetype = "dashed", size = 0.6) +
  theme_bw() +
  xlab('Number of Pixels') +
  ylab('Density') +
  ggtitle('Non-Water Interpreter Polygons')

# save plot
ggsave('D:/ontario_inventory/test/plots/land_interp.png')

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

################################################
###COMBINE MISSING PIXELS INTO SINGLE POLYGON###
################################################

# first result

# load polygon dataset
p <- vect('D:/ontario_inventory/test/scale_100/scale_100_ms_10_10_100.shp')

# subset by polygons that only have one pixel (NA) and polygons that have more
p_na <- p[p$nbPixels==1,]
p_real <- p[p$nbPixels>1,]

# dissolve polygons that only have 1 pixels
p2 <- aggregate(p_na, by='nbPixels')

# add back into single file
p3 <- rbind(p_real, p2)

# write to file to check in QGIS
writeVector(p3, 'D:/ontario_inventory/test/scale_100_ms_10_10_100_AGG.shp')
  
# load roads, watercourses and waterbodies
roads <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Roads/RMF_roads.shp')
waterc <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Lakes and rivers/RMF_watercourses.shp')
waterb <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Lakes and rivers/RMF_waterbodies.shp')

# try using the "erase" to remove polygons or parts of polygons that are roads, waterc, or waterb
p4 <- p3 %>% erase(., roads) %>% erase(., waterc) %>% erase(., waterb)

# write to file to check in QGIS
writeVector(p4, 'D:/ontario_inventory/test/scale_100_ms_10_10_100_AGG_ERASE.shp')

# second result

# load polygon dataset
p <- vect('D:/ontario_inventory/test/scale_100/scale_100_ms_10_15_50.shp')

# subset by polygons that only have one pixel (NA) and polygons that have more
p_na <- p[p$nbPixels==1,]
p_real <- p[p$nbPixels>1,]

# dissolve polygons that only have 1 pixels
p2 <- aggregate(p_na, by='nbPixels')

# add back into single file
p3 <- rbind(p_real, p2)

# write to file to check in QGIS
writeVector(p3, 'D:/ontario_inventory/test/scale_100_ms_10_15_50_AGG.shp', overwrite=T)

# two datasets with roads, rivers, and water bodies masked

# load polygon dataset
p <- vect('D:/ontario_inventory/test/scale_100/scale_m_100_ms_10_15_50.shp')

# subset by polygons that only have one pixel (NA) and polygons that have more
p_na <- p[p$nbPixels==1,]
p_real <- p[p$nbPixels>1,]

# dissolve polygons that only have 1 pixels
p2 <- aggregate(p_na, by='nbPixels')

# add back into single file
p3 <- rbind(p_real, p2)

# write to file to check in QGIS
writeVector(p3, 'D:/ontario_inventory/test/scale_m_100_ms_10_15_50_AGG.shp', overwrite = T)

# load polygon dataset
p <- vect('D:/ontario_inventory/test/scale_100/scale_m_100_ms_10_10_100.shp')

# subset by polygons that only have one pixel (NA) and polygons that have more
p_na <- p[p$nbPixels==1,]
p_real <- p[p$nbPixels>1,]

# dissolve polygons that only have 1 pixels
p2 <- aggregate(p_na, by='nbPixels')

# add back into single file
p3 <- rbind(p_real, p2)

# write to file to check in QGIS
writeVector(p3, 'D:/ontario_inventory/test/scale_m_100_ms_10_10_100_AGG.shp', overwrite = T)

########################################
###TRY A SPECTRAL CLUSTERING APPROACH###
########################################

# clear environment
# rm(list=ls())

# load example polygons from MS
p <- vect('D:/ontario_inventory/test/scale_100_ms_10_15_50_AGG.shp')

# get centroids of polygons
cent <- centroids(p)

# convert to df
cent_df <- as.data.frame(cent)

# remove row with NA values
cent_df <- subset(cent_df, nbPixels!=1)

#