# load packages
library(terra)
library(tidyverse)
library(exactextractr)
library(sf)
library(janitor)

# set parameters to loop through
n_o <- c("grm_15_01_09")
thresh_r <- c('15') 
spec_r <- c("0.1")
spat_r <- c("0.9")

# loop through segmentaiton
#for(i in seq_along(n_o)){
# for(i in 2:7){
  
  ##############################
  ###RUN MEAN SHIFT ALGORITHM###
  ##############################
  
  # load water and ucl polygons
  spl_pt <- vect('D:/ontario_inventory/segmentation/mask_wat_ucl/wat_ucl_polygons.shp')
  
  # SET PARAMETERS
  rast_in <- 'D:/ontario_inventory/segmentation/spl_stack_mask_wat_ucl_polytype.tif'
  out_p <- "D:/ontario_inventory/segmentation/grm"
  name_out <- n_o
  thresh <- thresh_r
  spec <- spec_r
  spat <- spat_r
  
  # set working directory where temp files will be output
  # later should add code to remove all temp files
  setwd('D:/temp')
  
  # create function to run generic region merging
  grm_otb <- function(otb_path = "", raster_in = "", out_path = "", name = "",
                      method = "bs", thresh = "", spec = "0.5", spat = "0.5"){
    
    # Set configuration
    conf <- paste("-in", raster_in, "-out", paste(out_path, "/", name, ".tif", sep=""),
                  "-criterion", method, "-threshold", thresh, "-cw", spec, "-sw", spat)
    
    # apply function in command line
    system(paste(otb_path, "/otbcli_GenericRegionMerging", " ", conf, sep=""))
    
    # save configuration for further use
    write.table(x = conf,file = paste(out_path,"/",name,"_conf.txt",sep=""),row.names = F, col.names = F)
  }
  
  # run mean shift
  grm_otb(otb_path = "C:/OTB/bin",
                raster_in = rast_in,
                out_path = out_p,
                name = name_out,
                thresh = thresh,
                spec = spec,
                spat = spat)
  
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
  
  ###########################
  ### MASK MISSING VALUES ###
  ###########################
  
  # load grm raster
  p <- rast(paste(out_p, "/", name_out, ".tif", sep=""))
  
  # load seg raster
  mask <- rast(rast_in) %>% .[[1]]
  
  # mask grm raster
  p <- mask(p, mask)
  
  # write grm raster
  writeRaster(p, paste(out_p, "/", name_out, ".tif", sep=""),
              overwrite = T)
  
  # convert to vector
  p_vec <- as.polygons(p)
  
  # write to disk
  writeVector(p_vec, paste(out_p, "/", name_out, ".shp", sep=""),
              overwrite = T)
  
#   ################################################
#   ### ADD PRE-ALLOCATED POLYGONS AND LANDCOVER ###
#   ################################################
#   
#   # load name in and name out 
#   name_i <- str_c(out_p, '/', name_out, '.shp')
#   name_o <- str_c(out_p, '/', name_out, '_add_pt.shp')
#   
#   # load polygon dataset
#   p <- vect(name_i)
#   
#   # mask out missing pixels
#   # p1 <- p[p$nbPixels > 1]
#   p <- p[is.na(p$meanB0) == F]
#   
#   # add non-FOR POLYTYPE polygons back in
#   p2 <- rbind(p, spl_pt)
#   
#   # load VLCE 2.0 landcover dataset from 2018
#   lc <- rast('D:/ontario_inventory/VLCE/CA_forest_VLCE2_2018_CLIPPED.tif')
#   
#   # project polygons to CRS of raster
#   p_lc <- project(p2, lc)
#   
#   # convert to sf
#   p_lcsf <- st_as_sf(p_lc)
#   
#   # extract landcover values
#   lc_vals <- exact_extract(lc, p_lcsf)
#   
#   # set landcover class key
#   lc_key <- c(`20` = 'Water',
#               `31` = 'Snow/Ice',
#               `32` = 'Rock/Rubble',
#               `33` = 'Exposed/Barren Land',
#               `40` = 'Bryoids',
#               `50` = 'Shrubland',
#               `80` = 'Wetland',
#               `81` = 'Wetland-Treed',
#               `100` = 'Herbs',
#               `210` = 'Coniferous',
#               `220` = 'Broadleaf',
#               `230` = 'Mixed Wood')
#   
#   # # load mode function
#   # get_mode <- function(x) {
#   #   ux <- unique(x)
#   #   ux[which.max(tabulate(match(x, ux)))]
#   # }
#   
#   # # load mode function for multiple modes
#   # get_mode2 <- function(x) {
#   #   ux <- unique(x)
#   #   tab <- tabulate(match(x, ux))
#   #   ux[tab == max(tab)]
#   # }
#   # 
#   # # set cov fraction
#   # cov_frac <- 0
#   
#   # find dominant lc type in each polygon
#   # if there are multiple modes keep them
#   # apply over list
#   lc_mode <- sapply(lc_vals, function(x){
#     x$value <- recode(x$value, !!!lc_key)
#     x <- x %>% group_by(value) %>% summarize(sum = sum(coverage_fraction))
#     m <- x$value[which(x$sum == max(x$sum))]
#     # m <- get_mode2(x$value[x$coverage_fraction >= cov_frac])
#     return(paste(m, collapse = " "))
#   })
#   
#   # add to polygon dataset
#   p2$dom_lc <- lc_mode
#   
#   # set landcover class key with single forested class
#   lc_key_for <- c(`20` = 'Water',
#                   `31` = 'Snow/Ice',
#                   `32` = 'Rock/Rubble',
#                   `33` = 'Exposed/Barren Land',
#                   `40` = 'Bryoids',
#                   `50` = 'Shrubland',
#                   `80` = 'Wetland',
#                   `81` = 'Forest',
#                   `100` = 'Herbs',
#                   `210` = 'Forest',
#                   `220` = 'Forest',
#                   `230` = 'Forest')
#   
#   # find pixels with forest at least 50% of pixel
#   # apply over list
#   lc_dom_for <- sapply(lc_vals, function(x){
#     x$value <- recode(x$value, !!!lc_key_for)
#     x <- x %>% group_by(value) %>% summarize(sum = sum(coverage_fraction))
#     m <- x$value[which(x$sum == max(x$sum))]
#     if((length(m) == 1) & (m == 'Forest')[1]){
#       if(x$sum[x$value == m]/sum(x$sum) >= 0.5){
#         return('Yes')
#       }else{return('No')}
#     }else{return('No')}
#   })
#   
#   # add to polygon dataset
#   p2$dom_for <- lc_dom_for
#   
#   # load p95 lidar values
#   p95 <- rast('D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_p95.tif')
#   
#   # project polygons to CRS of raster
#   p_p95 <- project(p2, p95)
#   
#   # convert to sf
#   p_p95sf <- st_as_sf(p_p95)
#   
#   #extract median values
#   p95_med <- exact_extract(p95, p_p95sf, fun = function(df){
#     median(df$value[df$coverage_fraction == 1], na.rm = T)
#   }, summarize_df = T)
#   
#   # add to polygon dataset
#   p2$p95 <- p95_med
#   
#   # write to file to check in QGIS
#   writeVector(p2, name_o,
#               overwrite = T)
#   
#   ##################################################################
#   ###COMBINE MISSING PIXELS FROM SEGMENTATION INTO SINGLE POLYGON###
#   ##################################################################
#   
#   # # subset by polygons that only have one pixel (NA) and polygons that have more
#   # p_na <- p[p$nbPixels==1,]
#   # p_real <- p[p$nbPixels>1,]
#   # 
#   # # dissolve polygons that only have 1 pixels
#   # p2 <- aggregate(p_na, by='nbPixels')
#   # 
#   # # add back into single file
#   # p3 <- rbind(p_real, p2)
#   # 
#   # # write to file to check in QGIS
#   # writeVector(p3, 'D:/ontario_inventory/segmentation/ms_10_10_100_for_polytype_agg_na.shp')
#   # 
#   # # remove variables
#   # rm(p, p_na, p_real, p2, p3)
#   
#   ############################################
#   ###EXTRACT POLYGON SUMMARY STATS AND PLOTS##
#   ############################################
#   
#   # extract polygon stats from qgis meanshift test runs
#   # load polygon and remove polygons with 1 pixel (eventually need to merge these maybe?)
#   # need a spatially contiguous file
#   p <- vect(name_o) %>% as.data.frame
#   
#   # create dataframe with values wanted
#   ms_df <- data.frame(name = name_out,
#                       min_ha = (min(p$nbPixels)*0.04) %>% round(2),
#                       mean_ha = (mean(p$nbPixels)*0.04) %>% round(2),
#                       med_ha = (median(p$nbPixels)*0.04) %>% round(2),
#                       max_ha = (max(p$nbPixels)*0.04) %>% round(2),
#                       num_poly = NROW(p))
#   
#   # plot density
#   ggplot(data.frame(nbPixels = p$nbPixels), aes(x = nbPixels*0.04)) +
#     geom_density() +
#     xlim(c(0,100)) +
#     ylim(c(0, 0.2)) +
#     geom_vline(aes(xintercept = median(nbPixels*0.04, na.rm = T)), 
#                linetype = "dashed", size = 0.6) +
#     theme_bw() +
#     xlab('Number of Hectares') +
#     ylab('Density') +
#     ggtitle(name_out) +
#     theme(text = element_text(size = 20))
#   
#   # save plot
#   ggsave(str_c('D:/ontario_inventory/segmentation/mask_wat_ucl/plots/', name_out, '.png'),
#          width = 2100, height = 2100, units = 'px')
#   
#   # subset reliable forested polygons
#   p2 <- p %>% filter(dom_for == 'Yes', is.na(POLYTYPE), p95 >= 5)
#   
#   # rbind df with forested values
#   ms_df <- rbind(ms_df, data.frame(name = str_c(name_out, '_for'),
#                                    min_ha = (min(p2$nbPixels)*0.04) %>% round(2),
#                                    mean_ha = (mean(p2$nbPixels)*0.04) %>% round(2),
#                                    med_ha = (median(p2$nbPixels)*0.04) %>% round(2),
#                                    max_ha = (max(p2$nbPixels)*0.04) %>% round(2),
#                                    num_poly = NROW(p2)))
#   
#   # plot density
#   ggplot(data.frame(nbPixels = p2$nbPixels), aes(x = nbPixels*0.04)) +
#     geom_density() +
#     xlim(c(0,100)) +
#     ylim(c(0, 0.2)) +
#     geom_vline(aes(xintercept = median(nbPixels*0.04, na.rm = T)), 
#                linetype = "dashed", size = 0.6) +
#     theme_bw() +
#     xlab('Number of Hectares') +
#     ylab('Density') +
#     ggtitle(str_c('Forested ', name_out, ' Polygons')) +
#     theme(text = element_text(size = 20))
#   
#   # save plot
#   ggsave(str_c('D:/ontario_inventory/segmentation/mask_wat_ucl/plots/', name_out, '_for.png'),
#          width = 2100, height = 2100, units = 'px')
#   
#   # subset non masked WAT and UCL polygons
#   p3 <- p %>% filter(is.na(POLYTYPE))
#   
#   # rbind df with forested values
#   ms_df <- rbind(ms_df, data.frame(name = str_c(name_out, '_mask_wat_ucl'),
#                                    min_ha = (min(p3$nbPixels)*0.04) %>% round(2),
#                                    mean_ha = (mean(p3$nbPixels)*0.04) %>% round(2),
#                                    med_ha = (median(p3$nbPixels)*0.04) %>% round(2),
#                                    max_ha = (max(p3$nbPixels)*0.04) %>% round(2),
#                                    num_poly = NROW(p3)))
#   
#   # plot density
#   ggplot(data.frame(nbPixels = p3$nbPixels), aes(x = nbPixels*0.04)) +
#     geom_density() +
#     xlim(c(0,100)) +
#     ylim(c(0, 0.2)) +
#     geom_vline(aes(xintercept = median(nbPixels*0.04, na.rm = T)), 
#                linetype = "dashed", size = 0.6) +
#     theme_bw() +
#     xlab('Number of Hectares') +
#     ylab('Density') +
#     ggtitle(str_c('Non WAT/UCL ', name_out, ' Polygons')) +
#     theme(text = element_text(size = 20))
#   
#   # save plot
#   ggsave(str_c('D:/ontario_inventory/segmentation/mask_wat_ucl/plots/', name_out, '_mask_wat_ucl.png'),
#          width = 2100, height = 2100, units = 'px')
#   
#   # load interpreter derived polygons to extract statistics
#   poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp') %>%
#     as.data.frame
#   
#   # create df of interpreter polygons
#   # from full interpreter dataset
#   int_df <- data.frame(name = 'interp_full',
#                        min_ha = (min(poly$AREA)/10000) %>% round(2),
#                        mean_ha = (mean(poly$AREA)/10000) %>% round(2),
#                        med_ha = (median(poly$AREA)/10000) %>% round(2),
#                        max_ha = (max(poly$AREA)/10000) %>% round(2),
#                        num_poly = NROW(poly))
#   
#   # subset forested polygons only
#   poly_for <- poly[poly$POLYTYPE=='FOR',]
#   
#   # from forested polygons only
#   int_df <- rbind(int_df,
#                   data.frame(name = 'interp_for',
#                              min_ha = (min(poly_for$AREA)/10000) %>% round(2),
#                              mean_ha = (mean(poly_for$AREA)/10000) %>% round(2),
#                              med_ha = (median(poly_for$AREA)/10000) %>% round(2),
#                              max_ha = (max(poly_for$AREA)/10000) %>% round(2),
#                              num_poly = NROW(poly_for)))
#   
#   # subset all non water/ucl polygons
#   poly_land <- poly[!(poly$POLYTYPE %in% c('WAT', 'UCL')),]
#   
#   # from non water polygons only
#   int_df <- rbind(int_df,
#                   data.frame(name = 'interp_non_wat_ucl',
#                              min_ha = (min(poly_land$AREA)/10000) %>% round(1),
#                              mean_ha = (mean(poly_land$AREA)/10000) %>% round(1),
#                              med_ha = (median(poly_land$AREA)/10000) %>% round(1),
#                              max_ha = (max(poly_land$AREA)/10000) %>% round(1),
#                              num_poly = NROW(poly_land)))
#   
#   # combine dfs from mean shift and interpreter
#   out_df <- rbind(ms_df, int_df)
#   
#   # write df as csv
#   write.csv(out_df,
#             file = str_c('D:/ontario_inventory/segmentation/mask_wat_ucl/plots/', name_out, '_poly_tab.csv'),
#             row.names = F)
#   
# #}
# 
