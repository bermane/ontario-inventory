# This code applies image segmentation to SPL data over the Nipissing Area in Ontario
# in order to create forest stand polygons

# This updated version masks roads and water (and other unclassified features)
# based on the FRI polytype in the dataset

# p95 is also being used instead of lorey's height to enable use where EFIs may 
# not be available and also since p95 is highly correlated to lorey's height.

# load packages
library(terra)
library(tidyverse)
library(exactextractr)
library(sf)
library(janitor)
library(berryFunctions)
library(lwgeom)
library(magrittr)

#################################################
###LOAD MULTI BAND SPL RASTER FOR SEGMENTATION###
#################################################

# set names of SPL rasters to stack
p95 <- 'D:/ontario_inventory/nipissing/als/NIPISSING_Metrics/NIPISSING_p95.tif'
cc <- 'D:/ontario_inventory/nipissing/als/NIPISSING_Metrics/NIPISSING_cov.tif'
cv <- 'D:/ontario_inventory/nipissing/als/NIPISSING_Metrics/NIPISSING_cv.tif'

# stack rasters
spl <- rast(c(p95, cc, cv))

# apply smoothing function on 5 cell square
spl[[1]] <- focal(spl[[1]], w=5, fun="mean")
spl[[2]] <- focal(spl[[2]], w=5, fun="mean")
spl[[3]] <- focal(spl[[3]], w=5, fun="mean")

# create spl template with all values equal to 1
spl_temp <- spl[[1]]
spl_temp[] <- 1

################
###MASK ROADS###
################

# load roads layer
roads <- vect('D:/ontario_inventory/nipissing/roads/NFRoad_Export.shp')

# reproject to match lidar
roads <- project(roads, spl)

# create roads polygon
spl_r <- mask(spl_temp, roads, touches = T)
npix <- sum(values(spl_r), na.rm = T)
spl_r <- as.polygons(spl_r)
names(spl_r) <- 'POLYTYPE'
spl_r$POLYTYPE <- 'RDS'
spl_r$nbPixels <- npix

# mask road and water body pixels to NA
spl <- spl %>%
  mask(., roads, inverse = T, touches = T)

#########################
###MASK WATER POLYGONS###
#########################

# for nipissing only masking water polygons since big UCL polygons
# will try masking roads separately from roads layer

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/nipissing/nipissing_fri.shp')

# subset polygons that are WAT OR UCL
poly_sub <- poly[poly$POLYTYPE %in% c('WAT')]

# reproject to match lidar
poly_sub <- project(poly_sub, spl)

# loop through non-forested polygons, mask raster, and vectorize
#for(i in seq_along(poly_sub)){
for(i in 1:length(poly_sub)){
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
    if(is.error(spl_temp %>% crop(poly_sub[i], snap = 'out') %>%
                mask(poly_sub[i], touches = T)) == F){
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
}

# reproject whole FRI to match lidar
poly <- project(poly, spl)

# mask lidar outside of FRI
spl <- mask(spl, poly, inverse = F, touches = T)

# mask WAT polygons
spl <- mask(spl, poly_sub, inverse = T, touches = T)

#############################################
###COMBINE ROAD AND WATER POLYGON DATASETS###
#############################################

spl_pt <- rbind(spl_pt, spl_r)

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
writeRaster(spl, filename = 'D:/ontario_inventory/nipissing/segmentation/spl_stack.tif', overwrite=T)

# save spl_pt
writeVector(spl_pt, 'D:/ontario_inventory/nipissing/segmentation/water_roads_polygons.shp',
            overwrite = T)

# remove variables
rm(cc, cv, p95, spl, poly, poly_sub, spl_hold, spl_temp, spl_pt, spl_r,
   i, npix, pt, roads)

#######################
###RUN GRM ALGORITHM###
#######################

# load water and roads polygons
spl_pt <- vect('D:/ontario_inventory/nipissing/segmentation/water_roads_polygons.shp')

# SET PARAMETERS
rast_in <- 'D:/ontario_inventory/nipissing/segmentation/spl_stack.tif'
out_p <- 'D:/ontario_inventory/nipissing/segmentation'
name_out <- 'grm_10_01_05_segments'
thresh <- "10"
spec <- "0.1"
spat <- "0.5"

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

# run grm
grm_otb(otb_path = "C:/OTB/bin",
        raster_in = rast_in,
        out_path = out_p,
        name = name_out,
        thresh = thresh,
        spec = spec,
        spat = spat)

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

# convert to vector based on cell value
vec <- as.polygons(p)

# create table of number of pixels in each polygon
num <- as.vector(values(p))
num_pix <- tabyl(num)

# drop na row
num_pix <- na.omit(num_pix)

# get pixel ids from vector
vec_dat <- tibble(id = values(vec)[,1])
colnames(vec_dat) <- 'id'

# loop through values and add to vector data
vec_dat$nbPixels <- NA
for(i in 1:NROW(vec_dat)){
  vec_dat$nbPixels[i] <- num_pix$n[num_pix$num == vec_dat$id[i]]
}

# remove current column of data and add id
# add nbPixels to vector
vec <- vec[,-1]
vec$id <- vec_dat$id
vec$nbPixels <- vec_dat$nbPixels

################################################
### ADD PRE-ALLOCATED POLYGONS AND LANDCOVER ###
################################################

# load polygon dataset
p <- vec

# load water and roads polygons
spl_pt <- vect('D:/ontario_inventory/nipissing/segmentation/water_roads_polygons.shp')

# reproject segmented polygons to ensure same crs
p <- project(p, spl_pt)

# add non-FOR POLYTYPE polygons back in
p2 <- rbind(p, spl_pt)

# load VLCE 2.0 landcover dataset from 2018
lc <- rast('D:/ontario_inventory/VLCE/CA_forest_VLCE2_2018.tif')

# project polygons to CRS of raster
p_lc <- project(p2, lc)

# crop raster
lc <- crop(lc, p_lc)

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
p95 <- rast('D:/ontario_inventory/nipissing/als/NIPISSING_Metrics/NIPISSING_p95.tif')

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

# load canopy cover lidar values
cc <- rast('D:/ontario_inventory/nipissing/als/NIPISSING_Metrics/NIPISSING_cov.tif')

# project polygons to CRS of raster
p_cc <- project(p2, cc)

# convert to sf
p_ccsf <- st_as_sf(p_cc)

#extract median values
cc_med <- exact_extract(cc, p_ccsf, fun = function(df){
  median(df$value[df$coverage_fraction == 1], na.rm = T)
}, summarize_df = T)

# add to polygon dataset
p2$cc <- cc_med

# convert to sf
p2_sf <- st_as_sf(p2)

# calculate perimeter
p2$perim <- st_perimeter(p2_sf)

# calculate area
p2$area <- st_area(p2_sf)

# calculate perimeter to area ratio
p2$p_to_a <- p2$perim/p2$area
p2$p_to_a <- round(p2$p_to_a, 3)

# write to file to check in QGIS
writeVector(p2, paste(str_replace(out_p, 'segments', 'shp'), "/", name_out, ".shp", sep=""),
            overwrite = T)

###########################################
### EXTRACT FINAL POLYGON SUMMARY STATS ###
###########################################

# create list of polygon files, names and parameters
file <- 'D:/ontario_inventory/nipissing/segmentation/grm_10_01_05_segments.shp'
out_loc <- 'D:/ontario_inventory/nipissing/segmentation/stats/'
grm_input <- 'D:/ontario_inventory/nipissing/segmentation/spl_stack.tif'
alg <- 'GRM'
param <- '10_0.1_0.5'
name <- '10_01_05'

# create standard error function
se <- function(x) sd(x)/sqrt(length(x))

  # load file
  p <- vect(file)
  
  # convert to sf
  p_sf <- st_as_sf(p)
  
  # calculate perimeter
  p$perim <- st_perimeter(p_sf) %>% as.numeric
  
  # calculate area
  p$area <- st_area(p_sf) %>% as.numeric
  
  # calculate perimeter to area ratio
  p$p_to_a <- p$perim/p$area
  p$p_to_a <- round(p$p_to_a, 3)
  
  # subset non masked WAT and RD polygons
  p2_sf <- p[is.na(p$POLYTYPE)] %>% st_as_sf
  p2 <- p[is.na(p$POLYTYPE)] %>% as.data.frame
  
  # subset reliable forested polygons
  # p3 <- p %>% filter(dom_for == 'Yes', is.na(POLYTYPE), p95 >= 5)
  
  # calculate msi
  p2$msi <- p2$perim/sqrt(pi * p2$area)
  
  # load original raster input file
  ras <- rast(grm_input)
  
  # rename bands
  names(ras) <- c('p95', 'cc', 'cv')
  
  # extract pixel values
  pvals <- exact_extract(ras, p2_sf)
  
  # calculate SSE
  sse <- sapply(pvals, FUN = function(x){
    p95_mean <- mean(x$p95, na.rm = T)
    cc_mean <- mean(x$cc, na.rm = T)
    cv_mean <- mean(x$cv, na.rm = T)
    
    return(c(sum((x$p95 - p95_mean)^2, na.rm = T),
             sum((x$cc - cc_mean)^2, na.rm = T),
             sum((x$cv - cv_mean)^2, na.rm = T)))
  })
  
  # transpose
  sse <- t(sse)
  
  # calculate final sums
  sse <- colSums(sse)
  
  # unlist values
  pvals2 <- do.call(rbind, pvals)
  
  # calculate global mean values
  p95_mean <- mean(pvals2$p95, na.rm = T)
  cc_mean <- mean(pvals2$cc, na.rm = T)
  cv_mean <- mean(pvals2$cv, na.rm = T)
  
  rm(pvals2)
  
  # calculate SST
  sst <- sapply(pvals, FUN = function(x){
    return(c(sum((x$p95 - p95_mean)^2, na.rm = T),
             sum((x$cc - cc_mean)^2, na.rm = T),
             sum((x$cv - cv_mean)^2, na.rm = T)))
  })
  
  # transpose
  sst <- t(sst)
  
  # calculate final sums
  sst <- colSums(sst)
  
  # calculate r2 values
  r2_p95 <- 1 - (sse[1]/sst[1]) %>% round(3)
  r2_cc <- 1 - (sse[2]/sst[2]) %>% round(3)
  r2_cv <- 1 - (sse[3]/sst[3]) %>% round(3)
  r2_all <- (sum(r2_p95, r2_cc, r2_cv) / 3) %>% round(3)
  
  # create dataframe with values wanted
  df <- data.frame(alg = alg,
                      param = param,
                      min_pix = (min(p2$nbPixels)),
                      max_pix = (max(p2$nbPixels)),
                      mean_pix = (mean(p2$nbPixels)),
                      med_pix = (median(p2$nbPixels)),
                      num_poly = NROW(p2),
                      mean_area = mean(p2$area),
                      se_area = se(p2$area),
                      sd_area = sd(p2$area),
                      mean_perim = mean(p2$perim),
                      se_perim = se(p2$perim),
                      sd_perim = sd(p2$perim),
                      mean_p_a = mean(p2$p_to_a),
                      se_p_a = se(p2$p_to_a),
                      sd_p_a = sd(p2$p_to_a),
                      mean_msi = mean(p2$msi),
                      se_msi = se(p2$msi),
                      sd_msi = sd(p2$msi),
                      r2_p95 = r2_p95,
                      r2_cc = r2_cc,
                      r2_cv = r2_cv,
                      r2_all = r2_all)
  
  # round numeric columns
  df %<>% 
    mutate_at(c('min_pix',
                'max_pix', 
                'mean_pix', 
                'med_pix', 
                'mean_area', 
                'se_area', 
                'sd_area',
                'mean_perim',
                'se_perim',
                'sd_perim'), 
              function(x) round(x, 2)) %>%
    mutate_at(c('mean_p_a', 
                'se_p_a', 
                'sd_p_a', 
                'mean_msi', 
                'se_msi', 
                'sd_msi'), 
              function(x) round(x, 4))
  
  #############################
  ### CALCULATE FINAL PLOTS ###
  #############################
  
  # plot density
  ggplot(data.frame(nbPixels = p2$nbPixels), aes(x = nbPixels)) +
    geom_density(fill = 'grey') +
    xlim(c(0,1000)) +
    ylim(c(0, 0.015)) +
    geom_vline(aes(xintercept = median(nbPixels)), 
               linetype = "dashed", size = 0.6) +
    theme_bw() +
    xlab('Number of Pixels') +
    ylab('Density') +
    ggtitle(str_c(alg)) +
    theme(text = element_text(size = 20))
  
  ggsave(str_c(out_loc, alg, '_final_pix_dens.png'),
         width = 2100, height = 2100, units = 'px')
  
  # # plot perimeter to area ratio
  # ggplot(data.frame(p_to_a = as.numeric(p2$p_to_a)), aes(x = p_to_a)) +
  #   geom_density() +
  #   xlim(c(0,0.2)) +
  #   ylim(c(0, 100)) +
  #   geom_vline(aes(xintercept = median(p_to_a)),
  #              linetype = "dashed", size = 0.6) +
  #   theme_bw() +
  #   xlab('Perimeter to Area Ratio') +
  #   ylab('Density') +
  #   ggtitle(str_c(alg[i], ' ', param[i])) +
  #   theme(text = element_text(size = 20))
  # 
  # ggsave(str_c('D:/ontario_inventory/segmentation/results_paper/plots/', alg[i], '_', name[i], '_ptoa_dens.png'),
  #        width = 2100, height = 2100, units = 'px')
  
  # plot shape index
  ggplot(data.frame(msi = as.numeric(p2$msi)), aes(x = msi)) +
    geom_density(fill = 'grey') +
    xlim(c(0,10)) +
    ylim(c(0, 1.5)) +
    geom_vline(aes(xintercept = median(msi)),
               linetype = "dashed", size = 0.6) +
    theme_bw() +
    xlab('Shape Index') +
    ylab('Density') +
    ggtitle(str_c(alg)) +
    theme(text = element_text(size = 20))
  
  ggsave(str_c(out_loc, alg, '_final_shape_dens.png'),
         width = 2100, height = 2100, units = 'px')
  
  #####################
  ### ADD FRI STATS ###
  #####################
  
  # load interpreter derived polygons to extract statistics
  p <- vect('D:/ontario_inventory/nipissing/nipissing_fri.shp')
  
  # convert to sf
  p_sf <- st_as_sf(p)
  
  # calculate perimeter
  p$perim <- st_perimeter(p_sf) %>% as.numeric
  
  # calculate area
  p$area <- st_area(p_sf) %>% as.numeric
  
  # calculate perimeter to area ratio
  p$p_to_a <- p$perim/p$area
  p$p_to_a <- round(p$p_to_a, 3)
  
  # subset all non water/ucl polygons
  p2_sf <- p[!(p$POLYTYPE %in% c('WAT', 'UCL'))] %>% st_as_sf
  p2 <- p[!(p$POLYTYPE %in% c('WAT', 'UCL'))] %>% as.data.frame
  
  # calculate msi
  p2$msi <- p2$perim/sqrt(pi * p2$area)
  
  # load original raster input file
  ras <- rast(grm_input)
  
  # rename bands
  names(ras) <- c('p95', 'cc', 'cv')
  
  # extract pixel values
  pvals <- exact_extract(ras, p2_sf)
  
  # calculate SSE
  sse <- sapply(pvals, FUN = function(x){
    
    # subset values based on coverage fraction
    x %<>% filter(coverage_fraction >= 0.5)
    
    p95_mean <- mean(x$p95, na.rm = T)
    cc_mean <- mean(x$cc, na.rm = T)
    cv_mean <- mean(x$cv, na.rm = T)
    
    return(c(sum((x$p95 - p95_mean)^2, na.rm = T),
             sum((x$cc - cc_mean)^2, na.rm = T),
             sum((x$cv - cv_mean)^2, na.rm = T)))
  })
  
  # transpose
  sse <- t(sse)
  
  # calculate final sums
  sse <- colSums(sse)
  
  # unlist values
  pvals2 <- do.call(rbind, pvals)
  
  # subset values based on coverage fraction
  pvals2 %<>% filter(coverage_fraction >= 0.5)
  
  # calculate global mean values
  p95_mean <- mean(pvals2$p95, na.rm = T)
  cc_mean <- mean(pvals2$cc, na.rm = T)
  cv_mean <- mean(pvals2$cv, na.rm = T)
  
  rm(pvals2)
  
  # calculate SST
  sst <- sapply(pvals, FUN = function(x){
    # subset values based on coverage fraction
    x %<>% filter(coverage_fraction >= 0.5)
    
    return(c(sum((x$p95 - p95_mean)^2, na.rm = T),
             sum((x$cc - cc_mean)^2, na.rm = T),
             sum((x$cv - cv_mean)^2, na.rm = T)))
  })
  
  # transpose
  sst <- t(sst)
  
  # calculate final sums
  sst <- colSums(sst)
  
  # calculate r2 values
  r2_p95 <- 1 - (sse[1]/sst[1]) %>% round(3)
  r2_cc <- 1 - (sse[2]/sst[2]) %>% round(3)
  r2_cv <- 1 - (sse[3]/sst[3]) %>% round(3)
  r2_all <- (sum(r2_p95, r2_cc, r2_cv) / 3) %>% round(3)
  
  # create dataframe with values wanted
  ms_df <- data.frame(alg = 'FRI',
                      param = 'FRI',
                      min_pix = (min(p2$area / 400)),
                      max_pix = (max(p2$area / 400)),
                      mean_pix = (mean(p2$area / 400)),
                      med_pix = (median(p2$area / 400)),
                      num_poly = NROW(p2),
                      mean_area = mean(p2$area),
                      se_area = se(p2$area),
                      sd_area = sd(p2$area),
                      mean_perim = mean(p2$perim),
                      se_perim = se(p2$perim),
                      sd_perim = sd(p2$perim),
                      mean_p_a = mean(p2$p_to_a),
                      se_p_a = se(p2$p_to_a),
                      sd_p_a = sd(p2$p_to_a),
                      mean_msi = mean(p2$msi),
                      se_msi = se(p2$msi),
                      sd_msi = sd(p2$msi),
                      r2_p95 = r2_p95,
                      r2_cc = r2_cc,
                      r2_cv = r2_cv,
                      r2_all = r2_all)
  
  # round numeric columns
  ms_df %<>% 
    mutate_at(c('min_pix',
                'max_pix', 
                'mean_pix', 
                'med_pix'),
              function(x) round(x)) %>%
    mutate_at(c('mean_area', 
                'se_area', 
                'sd_area',
                'mean_perim',
                'se_perim',
                'sd_perim'), 
              function(x) round(x, 2)) %>%
    mutate_at(c('mean_p_a', 
                'se_p_a', 
                'sd_p_a', 
                'mean_msi', 
                'se_msi', 
                'sd_msi'), 
              function(x) round(x, 4))
  
  # bind df
  df <- rbind(df, ms_df)

  # write df as csv
  write.csv(df,
            file = str_c(out_loc, 'summary_stats.csv'),
            row.names = F)
  
  #############################
  ### CALCULATE FINAL PLOTS ###
  #############################
  
  # plot density
  ggplot(data.frame(nbPixels = p2$area / 400), aes(x = nbPixels)) +
    geom_density(fill = 'grey') +
    xlim(c(0,1000)) +
    ylim(c(0, 0.015)) +
    geom_vline(aes(xintercept = median(nbPixels)), 
               linetype = "dashed", size = 0.6) +
    theme_bw() +
    xlab('Number of Pixels') +
    ylab('Density') +
    ggtitle('FRI') +
    theme(text = element_text(size = 20))
  
  ggsave(str_c(out_loc, 'fri_final_pix_dens.png'),
         width = 2100, height = 2100, units = 'px')
  
  # # plot perimeter to area ratio
  # ggplot(data.frame(p_to_a = as.numeric(p2$p_to_a)), aes(x = p_to_a)) +
  #   geom_density() +
  #   xlim(c(0,0.2)) +
  #   ylim(c(0, 100)) +
  #   geom_vline(aes(xintercept = median(p_to_a)),
  #              linetype = "dashed", size = 0.6) +
  #   theme_bw() +
  #   xlab('Perimeter to Area Ratio') +
  #   ylab('Density') +
  #   ggtitle('FRI') +
  #   theme(text = element_text(size = 20))
  # 
  # ggsave(str_c('D:/ontario_inventory/segmentation/results_paper/plots/fri_ptoa_dens.png'),
  #        width = 2100, height = 2100, units = 'px')
  
  # plot shape index
  ggplot(data.frame(msi = as.numeric(p2$msi)), aes(x = msi)) +
    geom_density(fill = 'grey') +
    xlim(c(0,10)) +
    ylim(c(0, 1.5)) +
    geom_vline(aes(xintercept = median(msi)),
               linetype = "dashed", size = 0.6) +
    theme_bw() +
    xlab('Shape Index') +
    ylab('Density') +
    ggtitle('FRI') +
    theme(text = element_text(size = 20))
  
  ggsave(str_c(out_loc, 'fri_final_shape_dens.png'),
         width = 2100, height = 2100, units = 'px')
  