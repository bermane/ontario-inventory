# This code applies image segmentation to SPL data over the Romeo Malette Forest in Ontario
# in order to create forest stand polygons

# This updated version screens for non-forested pixels in the LiDAR dataset using FRI polytype
# based on feedback from project partners

# p95 is also being used instead of lorey's height to enable use where EFIs may 
# not be available and also since p95 is highly correlated to lorey's height.

# load packages
library(terra)
library(tidyverse)
library(exactextractr)
library(sf)
library(janitor)

#################################################
###LOAD MULTI BAND SPL RASTER FOR SEGMENTATION###
#################################################

# set names of SPL rasters to stack
p95 <- 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_p95.tif'
cc <- 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif'
cv <- 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_cv.tif'

# stack rasters
spl <- rast(c(p95, cc, cv))

# apply smoothing function on 5 cell square
spl[[1]] <- focal(spl[[1]], w=5, fun="mean")
spl[[2]] <- focal(spl[[2]], w=5, fun="mean")
spl[[3]] <- focal(spl[[3]], w=5, fun="mean")

##########################
###MASK ROADS AND WATER###
##########################

# # load roads, watercourses and waterbodies
# roads <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Roads/RMF_roads.shp') %>%
#   project(., spl)
# waterc <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Lakes and rivers/RMF_watercourses.shp') %>%
#   project(., spl)
# waterb <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Lakes and rivers/RMF_waterbodies.shp') %>%
#   project(., spl)
# 
# # subset roads to only mask the main types used by interpreter
# # RDTYPE = H, P, B
# roads <- roads[roads$RDTYPE %in% c('H', 'P', 'B'),]

# create spl template with all values equal to 1
spl_temp <- spl[[1]]
spl_temp[] <- 1

# # create roads polygon
# spl_r <- mask(spl_temp, roads, touches = T)
# npix <- sum(values(spl_r), na.rm = T)
# spl_r <- as.polygons(spl_r)
# names(spl_r) <- 'POLYTYPE'
# spl_r$POLYTYPE <- 'RDS'
# spl_r$nbPixels <- npix
# 
# # create waterbodies polygon
# spl_w <- mask(spl_temp, waterb, touches = T)
# npix <- sum(values(spl_w), na.rm = T)
# spl_w <- as.polygons(spl_w)
# names(spl_w) <- 'POLYTYPE'
# spl_w$POLYTYPE <- 'WAT'
# spl_w$nbPixels <- npix
# 
# # mask road and water body pixels to NA
# spl <- spl %>% 
#   mask(., roads, inverse = T, touches = T) %>% 
#   mask(., waterb, inverse = T, touches = T)
# # mask(., waterc, inverse = T) %>% 

################################
###MASK NON FORESTED POLYGONS###
################################

# mask non-forested polygons from segmentation
# add them back in as polygons with their polytype

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

# subset polygons that are WAT OR UCL
poly_sub <- poly[poly$POLYTYPE %in% c('WAT', 'UCL')]

# reproject to match lidar
poly_sub <- project(poly_sub, spl)

# loop through non-forested polygons, mask raster, and vectorize
for(i in seq_along(poly_sub)){
  pt <- poly_sub$POLYTYPE[i]
  if(i == 1){
    spl_pt <- spl_temp %>% crop(poly_sub[i], snap = 'out') %>%
      mask(poly_sub[i], touches = T)
    npix <- sum(values(spl_pt), na.rm = T)
    spl_pt <- as.polygons(spl_pt)
    names(spl_pt) <- 'POLYTYPE'
    spl_pt$POLYTYPE <- pt
    spl_pt$nbPixels <- npix
  } else{
    spl_hold <- spl_temp %>% crop(poly_sub[i], snap = 'out') %>%
      mask(poly_sub[i], touches = T)
    npix <- sum(values(spl_hold), na.rm = T)
    spl_hold <- as.polygons(spl_hold)
    names(spl_hold) <- 'POLYTYPE'
    spl_hold$POLYTYPE <- pt
    spl_hold$nbPixels <- npix
    spl_pt <- rbind(spl_pt, spl_hold)
  }
}

# reproject whole FRI to match lidar
poly <- project(poly, spl)

# mask lidar outside of FRI
spl <- mask(spl, poly, inverse = F, touches = T)

# mask WAT and UCL polygons
spl <- mask(spl, poly_sub, inverse = T, touches = T)

##################################
###MASK USING LANDCOVER AND p95###
##################################

# # load VLCE 2.0 landcover dataset from 2018
# lc <- rast('D:/ontario_inventory/VLCE/CA_forest_VLCE2_2018_CLIPPED.tif')
# 
# # project lc to crs of raster
# lc <- project(lc, spl, method = 'near')
# 
# # define treed classes
# #treed <- c('Mixed Wood', 'Broadleaf', 'Coniferous', 'Wetland-Treed')
# treed <- c(230, 220, 210, 81)
# 
# # set spl values to NA if not treed
# spl[!(lc %in% treed)] <- NA
# 
# # set spl values to NA if p95 < 5
# spl[spl[[1]] < 5] <- NA

########################################
###DEAL WITH MISSING DATA AND RESCALE###
########################################

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
writeRaster(spl, filename='D:/ontario_inventory/segmentation/spl_stack_mask_wat_ucl_polytype.tif', overwrite=T)

# save spl_pt
writeVector(spl_pt, 'D:/ontario_inventory/segmentation/mask_wat_ucl/wat_ucl_polygons.shp',
            overwrite = T)

# remove variables
rm(cc, cv, p95, spl, poly, poly_sub, spl_hold, spl_temp, spl_pt)

##############################
###RUN MEAN SHIFT ALGORITHM###
##############################

# load water and ucl polygons
spl_pt <- vect('D:/ontario_inventory/segmentation/mask_wat_ucl/wat_ucl_polygons.shp')

# SET PARAMETERS
rast_in <- 'D:/ontario_inventory/segmentation/spl_stack_for_polytype.tif'
out_p <- "D:/ontario_inventory/segmentation/mask_wat_ucl"
name_out <- "ms_10_10_100"
spat_r <- "10"
rang_r <- "10"
m_size <- "100"

# set working directory where temp files will be output
setwd('D:/temp')

# create function to run mean shift
meanshift_otb <- function(otb_path = "", raster_in = "", out_path = "", name ="", spatialr = "10",
                           ranger = "10", minsize = "100", tilesizex = "500", tilesizey = "500",
                           outmode = "vector", cleanup = "true", ram = "256"){
  # Set configuration
  conf <- paste("-in", raster_in, "-spatialr", spatialr, "-ranger", ranger,
                "-minsize", minsize, "-tilesizex", tilesizex, "-tilesizey", tilesizey,
                "-mode", outmode, "-mode.vector.out", paste(out_path, "/", name, ".shp", sep=""),
                "-cleanup", cleanup,"-ram", ram)
  
  # apply function in command line
  system(paste(otb_path, "/otbcli_LargeScaleMeanShift", " ", conf, sep=""))
  
  # save configuration for further use
  write.table(x = conf,file = paste(out_path,"/",name,"_conf.txt",sep=""),row.names = F, col.names = F)
}

# run mean shift
meanshift_otb(otb_path = "C:/OTB/bin",
              raster_in = rast_in,
              out_path = out_p,
              name = name_out,
              spatialr = spat_r,
              ranger = rang_r,
              minsize = m_size,
              ram = "1024")

# # usage, you can set any option listed above
# meanshift_otb(otb_path = "C:/OTB/bin",
#               raster_in = 'D:/ontario_inventory/segmentation/spl_stack.tif',
#               out_path = "D:/ontario_inventory/segmentation",
#               name = "ms_10_15_50",
#               spatialr = "10",
#               ranger = "15",
#               minsize = "50",
#               ram = "1024")
# 
# # usage, you can set any option listed above
# meanshift_otb(otb_path = "C:/OTB/bin",
#               raster_in = 'D:/ontario_inventory/segmentation/spl_stack.tif',
#               out_path = "D:/ontario_inventory/segmentation",
#               name = "ms_10_10_50",
#               spatialr = "10",
#               ranger = "10",
#               minsize = "50",
#               ram = "1024")

################################################
### ADD PRE-ALLOCATED POLYGONS AND LANDCOVER ###
################################################

# load name in and name out 
name_i <- str_c(out_p, '/', name_out, '.shp')
name_o <- str_c(out_p, '/', name_out, '_add_pt.shp')

# load polygon dataset
p <- vect('D:/ontario_inventory/segmentation/ms_10_10_100_for_polytype_mask_wat_ucl.shp')

# mask out missing pixels
# p1 <- p[p$nbPixels > 1]
p <- p[is.na(p$meanB0) == F]

# add non-FOR POLYTYPE polygons back in
p2 <- rbind(p, spl_pt)

# load VLCE 2.0 landcover dataset from 2018
lc <- rast('D:/ontario_inventory/VLCE/CA_forest_VLCE2_2018_CLIPPED.tif')

# project polygons to CRS of raster
p_lc <- project(p2, lc)

# convert to sf
p_lcsf <- st_as_sf(p_lc)

# extract landcover values
lc_vals <- exact_extract(lc, p_lcsf)

# set landcover class key
lc_key <- c(`20` = 'Water',
            `31` = 'Snow/Ice',
            `32` = 'Rock/Rubble',
            `33` = 'Exposed/Barren Land',
            `40` = 'Bryoids',
            `50` = 'Shrubland',
            `80` = 'Wetland',
            `81` = 'Wetland-Treed',
            `100` = 'Herbs',
            `210` = 'Coniferous',
            `220` = 'Broadleaf',
            `230` = 'Mixed Wood')

# # load mode function
# get_mode <- function(x) {
#   ux <- unique(x)
#   ux[which.max(tabulate(match(x, ux)))]
# }

# # load mode function for multiple modes
# get_mode2 <- function(x) {
#   ux <- unique(x)
#   tab <- tabulate(match(x, ux))
#   ux[tab == max(tab)]
# }
# 
# # set cov fraction
# cov_frac <- 0

# find dominant lc type in each polygon
# if there are multiple modes keep them
# apply over list
lc_mode <- sapply(lc_vals, function(x){
  x$value <- recode(x$value, !!!lc_key)
  x <- x %>% group_by(value) %>% summarize(sum = sum(coverage_fraction))
  m <- x$value[which(x$sum == max(x$sum))]
  # m <- get_mode2(x$value[x$coverage_fraction >= cov_frac])
  return(paste(m, collapse = " "))
})

# add to polygon dataset
p2$dom_lc <- lc_mode

# set landcover class key with single forested class
lc_key_for <- c(`20` = 'Water',
                `31` = 'Snow/Ice',
                `32` = 'Rock/Rubble',
                `33` = 'Exposed/Barren Land',
                `40` = 'Bryoids',
                `50` = 'Shrubland',
                `80` = 'Wetland',
                `81` = 'Forest',
                `100` = 'Herbs',
                `210` = 'Forest',
                `220` = 'Forest',
                `230` = 'Forest')

# find pixels with forest at least 50% of pixel
# apply over list
lc_dom_for <- sapply(lc_vals, function(x){
  x$value <- recode(x$value, !!!lc_key_for)
  x <- x %>% group_by(value) %>% summarize(sum = sum(coverage_fraction))
  m <- x$value[which(x$sum == max(x$sum))]
  if((length(m) == 1) & (m == 'Forest')[1]){
    if(x$sum[x$value == m]/sum(x$sum) >= 0.5){
      return('Yes')
    }else{return('No')}
  }else{return('No')}
})

# add to polygon dataset
p2$dom_for <- lc_dom_for

# load p95 lidar values
p95 <- rast('D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_p95.tif')

# project polygons to CRS of raster
p_p95 <- project(p2, p95)

# convert to sf
p_p95sf <- st_as_sf(p_p95)

#extract median values
p95_med <- exact_extract(p95, p_p95sf, fun = function(df){
  median(df$value[df$coverage_fraction == 1], na.rm = T)
}, summarize_df = T)

# add to polygon dataset
p2$p95 <- p95_med

# write to file to check in QGIS
writeVector(p2, name_o,
            overwrite = T)

##################################################################
###COMBINE MISSING PIXELS FROM SEGMENTATION INTO SINGLE POLYGON###
##################################################################

# # subset by polygons that only have one pixel (NA) and polygons that have more
# p_na <- p[p$nbPixels==1,]
# p_real <- p[p$nbPixels>1,]
# 
# # dissolve polygons that only have 1 pixels
# p2 <- aggregate(p_na, by='nbPixels')
# 
# # add back into single file
# p3 <- rbind(p_real, p2)
# 
# # write to file to check in QGIS
# writeVector(p3, 'D:/ontario_inventory/segmentation/ms_10_10_100_for_polytype_agg_na.shp')
# 
# # remove variables
# rm(p, p_na, p_real, p2, p3)

############################################
###EXTRACT POLYGON SUMMARY STATS AND PLOTS##
############################################

# extract polygon stats from qgis meanshift test runs
# load polygon and remove polygons with 1 pixel (eventually need to merge these maybe?)
# need a spatially contiguous file
p <- vect(name_o) %>% as.data.frame

# create dataframe with values wanted
ms_df <- data.frame(name = name,
                 min_ha = (min(p$nbPixels)*0.04) %>% round(2),
                 mean_ha = (mean(p$nbPixels)*0.04) %>% round(2),
                 med_ha = (median(p$nbPixels)*0.04) %>% round(2),
                 max_ha = (max(p$nbPixels)*0.04) %>% round(2),
                 num_poly = NROW(p))

# plot density
ggplot(data.frame(nbPixels = p$nbPixels), aes(x = nbPixels*0.04)) +
  geom_density() +
  xlim(c(0,100)) +
  ylim(c(0, 0.2)) +
  geom_vline(aes(xintercept = median(nbPixels*0.04, na.rm = T)), 
             linetype = "dashed", size = 0.6) +
  theme_bw() +
  xlab('Number of Hectares') +
  ylab('Density') +
  ggtitle(name) +
  theme(text = element_text(size = 20))

# save plot
ggsave(str_c('D:/ontario_inventory/segmentation/mask_wat_ucl/plots/', name_out, '.png'),
       width = 2100, height = 2100, units = 'px')

# subset reliable forested polygons
p2 <- p %>% filter(dom_for == 'Yes', is.na(POLYTYPE), p95 >= 5)

# rbind df with forested values
ms_df <- rbind(ms_df, data.frame(name = str_c(name, '_for'),
                           min_ha = (min(p2$nbPixels)*0.04) %>% round(2),
                           mean_ha = (mean(p2$nbPixels)*0.04) %>% round(2),
                           med_ha = (median(p2$nbPixels)*0.04) %>% round(2),
                           max_ha = (max(p2$nbPixels)*0.04) %>% round(2),
                           num_poly = NROW(p2)))

# plot density
ggplot(data.frame(nbPixels = p2$nbPixels), aes(x = nbPixels*0.04)) +
  geom_density() +
  xlim(c(0,100)) +
  ylim(c(0, 0.2)) +
  geom_vline(aes(xintercept = median(nbPixels*0.04, na.rm = T)), 
             linetype = "dashed", size = 0.6) +
  theme_bw() +
  xlab('Number of Hectares') +
  ylab('Density') +
  ggtitle(str_c('Forested ', name, ' Polygons')) +
  theme(text = element_text(size = 20))

# save plot
ggsave(str_c('D:/ontario_inventory/segmentation/mask_wat_ucl/plots/', name_out, '_for.png'),
       width = 2100, height = 2100, units = 'px')

# subset non masked WAT and UCL polygons
p3 <- p %>% filter(is.na(POLYTYPE))

# rbind df with forested values
ms_df <- rbind(ms_df, data.frame(name = str_c(name, '_mask_wat_ucl'),
                           min_ha = (min(p3$nbPixels)*0.04) %>% round(2),
                           mean_ha = (mean(p3$nbPixels)*0.04) %>% round(2),
                           med_ha = (median(p3$nbPixels)*0.04) %>% round(2),
                           max_ha = (max(p3$nbPixels)*0.04) %>% round(2),
                           num_poly = NROW(p3)))

# plot density
ggplot(data.frame(nbPixels = p3$nbPixels), aes(x = nbPixels*0.04)) +
  geom_density() +
  xlim(c(0,100)) +
  ylim(c(0, 0.2)) +
  geom_vline(aes(xintercept = median(nbPixels*0.04, na.rm = T)), 
             linetype = "dashed", size = 0.6) +
  theme_bw() +
  xlab('Number of Hectares') +
  ylab('Density') +
  ggtitle(str_c('Non WAT/UCL ', name, ' Polygons')) +
  theme(text = element_text(size = 20))

# save plot
ggsave(str_c('D:/ontario_inventory/segmentation/mask_wat_ucl/plots/', name_out, '_mask_wat_ucl.png'),
       width = 2100, height = 2100, units = 'px')

# load interpreter derived polygons to extract statistics
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp') %>%
  as.data.frame

# create df of interpreter polygons
# from full interpreter dataset
int_df <- data.frame(name = 'interp_full',
                     min_ha = (min(poly$AREA)/10000) %>% round(2),
                     mean_ha = (mean(poly$AREA)/10000) %>% round(2),
                     med_ha = (median(poly$AREA)/10000) %>% round(2),
                     max_ha = (max(poly$AREA)/10000) %>% round(2),
                     num_poly = NROW(poly))

# subset forested polygons only
poly_for <- poly[poly$POLYTYPE=='FOR',]

# from forested polygons only
int_df <- rbind(int_df,
                data.frame(name = 'interp_for',
                           min_ha = (min(poly_for$AREA)/10000) %>% round(2),
                           mean_ha = (mean(poly_for$AREA)/10000) %>% round(2),
                           med_ha = (median(poly_for$AREA)/10000) %>% round(2),
                           max_ha = (max(poly_for$AREA)/10000) %>% round(2),
                           num_poly = NROW(poly_for)))

# subset all non water/ucl polygons
poly_land <- poly[!(poly$POLYTYPE %in% c('WAT', 'UCL')),]

# from non water polygons only
int_df <- rbind(int_df,
                data.frame(name = 'interp_non_wat_ucl',
                           min_ha = (min(poly_land$AREA)/10000) %>% round(1),
                           mean_ha = (mean(poly_land$AREA)/10000) %>% round(1),
                           med_ha = (median(poly_land$AREA)/10000) %>% round(1),
                           max_ha = (max(poly_land$AREA)/10000) %>% round(1),
                           num_poly = NROW(poly_land)))

# combine dfs from mean shift and interpreter
out_df <- rbind(ms_df, int_df)

# write df as csv
write.csv(out_df,
          file = str_c('D:/ontario_inventory/segmentation/mask_wat_ucl/plots/', name_out, '_poly_tab.csv'),
          row.names = F)

# output histograms from interpreter polygons
# area hist all poly types
ggplot(poly, aes(x = AREA/10000)) +
  geom_density() +
  xlim(c(0,100)) +
  ylim(c(0, 0.2)) +
  geom_vline(aes(xintercept = median(AREA/10000, na.rm = T)),
             linetype = "dashed", size = 0.6) +
  theme_bw() +
  xlab('Number of Hectares') +
  ylab('Density') +
  ggtitle('All Interpreter Polygons') +
  theme(text = element_text(size = 20))

# save plot
ggsave('D:/ontario_inventory/segmentation/plots/all_interp.png',
       width = 2100, height = 2100, units = 'px')

# area hist only forest polys
ggplot(poly_for, aes(x = AREA/10000)) +
  geom_density() +
  xlim(c(0,100)) +
  ylim(c(0, 0.2)) +
  geom_vline(aes(xintercept = median(AREA/10000, na.rm = T)),
             linetype = "dashed", size = 0.6) +
  theme_bw() +
  xlab('Number of Hectares') +
  ylab('Density') +
  ggtitle('Forested Interpreter Polygons') +
  theme(text = element_text(size = 20))

# save plot
ggsave('D:/ontario_inventory/segmentation/plots/forested_interp.png',
       width = 2100, height = 2100, units = 'px')

# area hist non water/ucl polys
ggplot(poly_land, aes(x = AREA/10000)) +
  geom_density() +
  xlim(c(0,100)) +
  ylim(c(0, 0.2)) +
  geom_vline(aes(xintercept = median(AREA/10000, na.rm = T)),
             linetype = "dashed", size = 0.6) +
  theme_bw() +
  xlab('Number of Hectares') +
  ylab('Density') +
  ggtitle('Non-Water/Ucl Interpreter Polygons') +
  theme(text = element_text(size = 20))

# save plot
ggsave('D:/ontario_inventory/segmentation/plots/non_wat_ucl_interp.png',
       width = 2100, height = 2100, units = 'px')

#remove variables
rm(poly, poly_for, poly_land)

##################################################
### ANALYZE FORESTED POLYGONS NON-FOR POLYTYPE ###
##################################################

# load polygons
poly <- vect('D:/ontario_inventory/segmentation/ms_10_10_100_for_polytype_mask_wat_ucl_add_pt.shp')

# convert to df
p <- as.data.frame(poly)

# pre allocate lists
dom_lc <- list()
dom_for <- list()

# loop through unique POLYTYPES to create table of dom_lc and dom_for
for(i in unique(p$POLYTYPE)){
  
  # subset by polytype
  psub <- p %>% filter(POLYTYPE == i)
  
  # add dom_lc and dom_for tables to list
  dom_lc[[i]] <- psub %>% tabyl(dom_lc) %>% arrange(desc(n)) %>% mutate(percent = round(percent*100,2))
  dom_for[[i]] <- psub %>% tabyl(dom_for) %>% arrange(desc(n)) %>% mutate(percent = round(percent*100,2))
  
}
