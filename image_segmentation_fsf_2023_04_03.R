# This code applies image segmentation to SPL data in Ontario
# in order to create a vector-based inventory of the landscape

##################################
### INSTALL PACKAGES IF NEEDED ###
##################################

# install.packages(c('terra',
#                    'tidyverse',
#                    'exactextractr',
#                    'sf',
#                    'janitor',
#                    'berryFunctions',
#                    'lwgeom',
#                    'magrittr'))

# make sure to have OTB installed from here:
# https://www.orfeo-toolbox.org/

####################################
### SET CODE AND FILE PARAMETERS ###
####################################

# set file names for SPL input variables
# these should be gridded raster data with the same CRS and extent
# p95 = 95th percentile of Height returns
# cc = canopy cover (% of FIRST returns above 2 m)
# cv = coefficient of variation (sd of height/mean height)
p95_f <- 'D:/ontario_inventory/FSF/ALS/zq95.img'
cc_f <- 'D:/ontario_inventory/FSF/ALS/cov_2m.img'
cv_f <- 'D:/ontario_inventory/FSF/ALS/cv.img'

# set file location of roads shape file (spatial lines)
# roads will be polygonized, masked from segmentation
# and re-added to final dataset as polygons
roads_f <-
  'D:/ontario_inventory/FSF/roads/FS_roads.shp'

# set file location of FRI polygons shape file
# FRI POLYTYPE should have a "WAT" classification to mask water polygons
fri <- 'D:/ontario_inventory/FSF/FRI/FSF_opi_polygon_CSRS_NAD83_17.shp'

# set output folder for files generated
# make sure no "/" at end of folder location!
out_dir <- 'C:/Users/bermane/Desktop/FSF'

# set folder location of OTB (where you installed OTB earlier)
otb_dir <- "C:/OTB/bin"

# set GRM segmentation parameters
# the default are listed below
# refer to paper or OTB GRM webpage for description of params
thresh <- "10"
spec <- "0.1"
spat <- "0.5"

# set file location of 2018 VLCE 2.0 landcover data
# using 2018 because it is the year of Romeo ALS acquisition
# can change based on ALS acquisition year
# download here:
# https://opendata.nfis.org/mapserver/nfis-change_eng.html
lc_f <- 'D:/ontario_inventory/VLCE/CA_forest_VLCE2_2018.tif'

###########################################################################
###########################################################################
### DO NOT EDIT CODE BELOW HERE !!! ### DO NOT EDIT CODE BELOW HERE !!! ###
###########################################################################
###########################################################################

#####################
### LOAD PACKAGES ###
#####################

# load packages
library(terra)
library(tidyverse)
library(exactextractr)
library(sf)
library(janitor)
library(berryFunctions)
library(lwgeom)
library(magrittr)

###################################################
### LOAD MULTI BAND SPL RASTER FOR SEGMENTATION ###
###################################################

# stack rasters
spl <- rast(c(p95_f, cc_f, cv_f))

# apply smoothing function on 5 cell square
spl[[1]] <- focal(spl[[1]], w = 5, fun = "mean")
spl[[2]] <- focal(spl[[2]], w = 5, fun = "mean")
spl[[3]] <- focal(spl[[3]], w = 5, fun = "mean")

# create spl template with all values equal to 1
spl_temp <- spl[[1]]
spl_temp[] <- 1

##################
### MASK ROADS ###
##################

# load roads layer
roads <- vect(roads_f)

# reproject to match lidar
roads <- project(roads, spl)

# create roads polygon
spl_r <- mask(spl_temp, roads, touches = T)
npix <- sum(values(spl_r), na.rm = T)
spl_r <- as.polygons(spl_r)
names(spl_r) <- 'POLYTYPE'
spl_r$POLYTYPE <- 'RDS'
spl_r$nbPixels <- npix

# mask road pixels to NA
spl <- spl %>%
  mask(., roads, inverse = T, touches = T)

###########################
### MASK WATER POLYGONS ###
###########################

# water polygons from the FRI are masked and readded after segmenation

# load photo interpreted polygons
poly <- vect(fri)

# subset polygons that are WAT
poly_sub <- poly[poly$POLYTYPE %in% c('WAT')]

# reproject to match lidar
poly_sub <- project(poly_sub, spl)

# loop through water polygons, mask raster, and vectorize
for (i in 1:length(poly_sub)) {
  pt <- poly_sub$POLYTYPE[i]
  if (i == 1) {
    spl_pt <- spl_temp %>% crop(poly_sub[i], snap = 'out') %>%
      mask(poly_sub[i], touches = T)
    npix <- sum(values(spl_pt), na.rm = T)
    spl_pt <- as.polygons(spl_pt)
    names(spl_pt) <- 'POLYTYPE'
    spl_pt$POLYTYPE <- pt
    spl_pt$nbPixels <- npix
  } else{
    if (is.error(spl_temp %>% crop(poly_sub[i], snap = 'out') %>%
                 mask(poly_sub[i], touches = T)) == F) {
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

###############################################
### COMBINE ROAD AND WATER POLYGON DATASETS ###
###############################################

spl_pt <- rbind(spl_pt, spl_r)

##########################################
### DEAL WITH MISSING DATA AND RESCALE ###
##########################################

# if any band is missing values set all to NA
spl[is.na(spl[[1]])] <- NA
spl[is.na(spl[[2]])] <- NA
spl[is.na(spl[[3]])] <- NA

# create function to rescale values from 0 to 100 using 1 and 99 percentile
scale_100 <- function(x) {
  # calculate 1st and 99th percentile of input raster
  perc <-
    values(x, mat = F) %>% quantile(., probs = c(0.01, 0.99), na.rm = T)
  
  # rescale raster using 1st and 99th %
  x <- (x - perc[1]) / (perc[2] - perc[1]) * 100
  
  #reset values below 0 and above 100
  x[x < 0] <- 0
  x[x > 100] <- 100
  
  return(x)
}

# rescale rasters from 0 to 100
spl[[1]] <- scale_100(spl[[1]])
spl[[2]] <- scale_100(spl[[2]])
spl[[3]] <- scale_100(spl[[3]])

# check if temp dir exists and create
if (dir.exists(file.path(out_dir, 'temp')) == F) {
  dir.create(file.path(out_dir, 'temp'))
}

# write raster to tif
writeRaster(spl,
            filename = str_c(out_dir, '/temp/spl_stack.tif'),
            overwrite = T)

# write spl_pt
writeVector(spl_pt,
            str_c(out_dir, '/temp/water_roads_polygons.shp'),
            overwrite = T)

#######################
###RUN GRM ALGORITHM###
#######################

# SET PARAMETERS
rast_in <- str_c(out_dir, '/temp/spl_stack.tif')
out_p <- str_c(out_dir, '/temp')
name_out <- str_c(
  'grm_',
  thresh,
  '_',
  gsub(".", "", spec, fixed = TRUE),
  '_',
  gsub(".", "", spat, fixed = TRUE)
)

# create function to run generic region merging
grm_otb <-
  function(otb_path = "",
           raster_in = "",
           out_path = "",
           name = "",
           method = "bs",
           thresh = "",
           spec = "0.5",
           spat = "0.5") {
    # Set configuration
    conf <-
      paste(
        "-in",
        raster_in,
        "-out",
        paste(out_path, "/", name, ".tif", sep = ""),
        "-criterion",
        method,
        "-threshold",
        thresh,
        "-cw",
        spec,
        "-sw",
        spat
      )
    
    # apply function in command line
    system(paste(otb_path, "/otbcli_GenericRegionMerging", " ", conf, sep =
                   ""))
    
    # save configuration for further use
    write.table(
      x = conf,
      file = paste(out_path, "/", name, "_conf.txt", sep = ""),
      row.names = F,
      col.names = F
    )
  }

# run grm
grm_otb(
  otb_path = otb_dir,
  raster_in = rast_in,
  out_path = out_p,
  name = name_out,
  thresh = thresh,
  spec = spec,
  spat = spat
)

###########################
### MASK MISSING VALUES ###
###########################

# load grm raster
p <- rast(str_c(out_p, "/", name_out, ".tif"))

# load seg raster
mask <- rast(rast_in) %>% .[[1]]

# mask grm raster
p <- mask(p, mask)

# write grm raster
writeRaster(p, paste(out_p, "/", name_out, ".tif", sep = ""),
            overwrite = T)

# convert to vector based on cell value
vec <- as.polygons(p)

# create table of number of pixels in each polygon
num <- as.vector(values(p))
num_pix <- tabyl(num)

# drop na row
num_pix <- na.omit(num_pix)

# get pixel ids from vector
vec_dat <- tibble(id = values(vec)[, 1])
colnames(vec_dat) <- 'id'

# loop through values and add to vector data
vec_dat$nbPixels <- NA
for (i in 1:NROW(vec_dat)) {
  vec_dat$nbPixels[i] <- num_pix$n[num_pix$num == vec_dat$id[i]]
}

# remove current column of data and add id
# add nbPixels to vector
vec <- vec[, -1]
vec$id <- vec_dat$id
vec$nbPixels <- vec_dat$nbPixels

##################################
### ADD PRE-ALLOCATED POLYGONS ###
##################################

# load polygon dataset
p <- vec

# reproject segmented polygons to ensure same crs
p <- project(p, spl_pt)

# add non-FOR POLYTYPE polygons back in
p2 <- rbind(p, spl_pt)

#####################
### ADD LANDCOVER ###
#####################

# load VLCE 2.0 landcover dataset
lc <- rast(lc_f)

# project polygons to CRS of raster
p_lc <- project(p2, lc)

# crop raster
lc <- crop(lc, p_lc)

# convert to sf
p_lcsf <- st_as_sf(p_lc)

# extract landcover values
lc_vals <- exact_extract(lc, p_lcsf)

# set landcover class key
lc_key <- c(`0` = 'NA',
            `20` = 'Water',
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
lc_key_for <- c(`0` = 'NA',
                `20` = 'Water',
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

# ######################
# ### ADD P95 and CC ###
# ######################
#
# # load p95 lidar values
# p95 <- rast(p95_f)
#
# # project polygons to CRS of raster
# p_p95 <- project(p2, p95)
#
# # convert to sf
# p_p95sf <- st_as_sf(p_p95)
#
# #extract median values
# p95_med <- exact_extract(p95, p_p95sf, fun = function(df){
#   median(df$value[df$coverage_fraction == 1], na.rm = T)
# }, summarize_df = T)
#
# # add to polygon dataset
# p2$p95 <- p95_med
#
# # load canopy cover lidar values
# cc <- rast(cc_f)
#
# # project polygons to CRS of raster
# p_cc <- project(p2, cc)
#
# # convert to sf
# p_ccsf <- st_as_sf(p_cc)
#
# #extract median values
# cc_med <- exact_extract(cc, p_ccsf, fun = function(df){
#   median(df$value[df$coverage_fraction == 1], na.rm = T)
# }, summarize_df = T)
#
# # add to polygon dataset
# p2$cc <- cc_med

##############################
### ADD AREA AND PERIMETER ###
##############################

# convert to sf
p2_sf <- st_as_sf(p2)

# calculate perimeter
p2$perim <- st_perimeter(p2_sf) %>% as.numeric

# calculate area
p2$area <- st_area(p2_sf) %>% as.numeric

# write to file
writeVector(p2, str_c(out_dir, "/", name_out, ".shp"),
            overwrite = T)

###########################################
### EXTRACT FINAL POLYGON SUMMARY STATS ###
###########################################

# create list of polygon files, names and parameters
file <- str_c(out_dir, "/", name_out, ".shp")
out_loc <- out_dir
grm_input <- str_c(out_dir, '/temp/spl_stack.tif')
name <- name_out

# create standard error function
se <- function(x)
  sd(x) / sqrt(length(x))

# load file
p <- vect(file)

# convert to sf
p_sf <- st_as_sf(p)

# subset non masked WAT and RD polygons
p2_sf <- p[is.na(p$POLYTYPE)] %>% st_as_sf
p2 <- p[is.na(p$POLYTYPE)] %>% as.data.frame

# calculate perimeter to area ratio
p2$p_to_a <- p2$perim / p2$area
p2$p_to_a <- round(p2$p_to_a, 3)

# calculate msi
p2$msi <- p2$perim / sqrt(pi * p2$area)

# load original raster input file
ras <- rast(grm_input)

# rename bands
names(ras) <- c('p95', 'cc', 'cv')

# extract pixel values
pvals <- exact_extract(ras, p2_sf)

# calculate SSE
sse <- sapply(
  pvals,
  FUN = function(x) {
    p95_mean <- mean(x$p95, na.rm = T)
    cc_mean <- mean(x$cc, na.rm = T)
    cv_mean <- mean(x$cv, na.rm = T)
    
    return(c(sum((x$p95 - p95_mean) ^ 2, na.rm = T),
             sum((x$cc - cc_mean) ^ 2, na.rm = T),
             sum((x$cv - cv_mean) ^ 2, na.rm = T)))
  }
)

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
sst <- sapply(
  pvals,
  FUN = function(x) {
    return(c(sum((x$p95 - p95_mean) ^ 2, na.rm = T),
             sum((x$cc - cc_mean) ^ 2, na.rm = T),
             sum((x$cv - cv_mean) ^ 2, na.rm = T)))
  }
)

# transpose
sst <- t(sst)

# calculate final sums
sst <- colSums(sst)

# calculate r2 values
r2_p95 <- 1 - (sse[1] / sst[1]) %>% round(3)
r2_cc <- 1 - (sse[2] / sst[2]) %>% round(3)
r2_cv <- 1 - (sse[3] / sst[3]) %>% round(3)
r2_all <- (sum(r2_p95, r2_cc, r2_cv) / 3) %>% round(3)

# create dataframe with values wanted
df <- data.frame(
  alg = name,
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
  r2_all = r2_all
)

# round numeric columns
df %<>%
  mutate_at(c(
    'min_pix',
    'max_pix',
    'mean_pix',
    'med_pix',
    'mean_area',
    'se_area',
    'sd_area',
    'mean_perim',
    'se_perim',
    'sd_perim'
  ),
  function(x)
    round(x, 2)) %>%
  mutate_at(c('mean_p_a',
              'se_p_a',
              'sd_p_a',
              'mean_msi',
              'se_msi',
              'sd_msi'),
            function(x)
              round(x, 4))

#############################
### CALCULATE FINAL PLOTS ###
#############################

# plot density
ggplot(data.frame(nbPixels = p2$nbPixels), aes(x = nbPixels)) +
  geom_density(fill = 'grey') +
  xlim(c(0, 1000)) +
  ylim(c(0, 0.015)) +
  geom_vline(aes(xintercept = median(nbPixels)),
             linetype = "dashed",
             linewidth = 0.6) +
  theme_bw() +
  xlab('Number of Pixels') +
  ylab('Density') +
  ggtitle(str_c(name)) +
  theme(text = element_text(size = 20))

ggsave(
  str_c(out_loc, '/', name, '_final_pix_dens.png'),
  width = 2100,
  height = 2100,
  units = 'px'
)

# plot shape index
ggplot(data.frame(msi = as.numeric(p2$msi)), aes(x = msi)) +
  geom_density(fill = 'grey') +
  xlim(c(0, 10)) +
  ylim(c(0, 1.5)) +
  geom_vline(aes(xintercept = median(msi)),
             linetype = "dashed",
             linewidth = 0.6) +
  theme_bw() +
  xlab('Shape Index') +
  ylab('Density') +
  ggtitle(str_c(name)) +
  theme(text = element_text(size = 20))

ggsave(
  str_c(out_loc, '/', name, '_final_shape_dens.png'),
  width = 2100,
  height = 2100,
  units = 'px'
)

#####################
### ADD FRI STATS ###
#####################

# load interpreter derived polygons to extract statistics
p <- vect(fri)

# convert to sf
p_sf <- st_as_sf(p)

# calculate perimeter
p$perim <- st_perimeter(p_sf) %>% as.numeric

# calculate area
p$area <- st_area(p_sf) %>% as.numeric

# calculate perimeter to area ratio
p$p_to_a <- p$perim / p$area
p$p_to_a <- round(p$p_to_a, 3)

# subset all non water/ucl polygons
p2_sf <- p[!(p$POLYTYPE %in% c('WAT', 'UCL'))] %>% st_as_sf
p2 <- p[!(p$POLYTYPE %in% c('WAT', 'UCL'))] %>% as.data.frame

# calculate msi
p2$msi <- p2$perim / sqrt(pi * p2$area)

# load original raster input file
ras <- rast(grm_input)

# rename bands
names(ras) <- c('p95', 'cc', 'cv')

# extract pixel values
pvals <- exact_extract(ras, p2_sf)

# calculate SSE
sse <- sapply(
  pvals,
  FUN = function(x) {
    # subset values based on coverage fraction
    x %<>% filter(coverage_fraction >= 0.5)
    
    p95_mean <- mean(x$p95, na.rm = T)
    cc_mean <- mean(x$cc, na.rm = T)
    cv_mean <- mean(x$cv, na.rm = T)
    
    return(c(sum((x$p95 - p95_mean) ^ 2, na.rm = T),
             sum((x$cc - cc_mean) ^ 2, na.rm = T),
             sum((x$cv - cv_mean) ^ 2, na.rm = T)))
  }
)

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
sst <- sapply(
  pvals,
  FUN = function(x) {
    # subset values based on coverage fraction
    x %<>% filter(coverage_fraction >= 0.5)
    
    return(c(sum((x$p95 - p95_mean) ^ 2, na.rm = T),
             sum((x$cc - cc_mean) ^ 2, na.rm = T),
             sum((x$cv - cv_mean) ^ 2, na.rm = T)))
  }
)

# transpose
sst <- t(sst)

# calculate final sums
sst <- colSums(sst)

# calculate r2 values
r2_p95 <- 1 - (sse[1] / sst[1]) %>% round(3)
r2_cc <- 1 - (sse[2] / sst[2]) %>% round(3)
r2_cv <- 1 - (sse[3] / sst[3]) %>% round(3)
r2_all <- (sum(r2_p95, r2_cc, r2_cv) / 3) %>% round(3)

# create dataframe with values wanted
ms_df <- data.frame(
  alg = 'FRI',
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
  r2_all = r2_all
)

# round numeric columns
ms_df %<>%
  mutate_at(c('min_pix',
              'max_pix',
              'mean_pix',
              'med_pix'),
            function(x)
              round(x)) %>%
  mutate_at(c(
    'mean_area',
    'se_area',
    'sd_area',
    'mean_perim',
    'se_perim',
    'sd_perim'
  ),
  function(x)
    round(x, 2)) %>%
  mutate_at(c('mean_p_a',
              'se_p_a',
              'sd_p_a',
              'mean_msi',
              'se_msi',
              'sd_msi'),
            function(x)
              round(x, 4))

# bind df
df <- rbind(df, ms_df)

# write df as csv
write.csv(df,
          file = str_c(out_loc, '/summary_stats.csv'),
          row.names = F)

#############################
### CALCULATE FINAL PLOTS ###
#############################

# plot density
ggplot(data.frame(nbPixels = p2$area / 400), aes(x = nbPixels)) +
  geom_density(fill = 'grey') +
  xlim(c(0, 1000)) +
  ylim(c(0, 0.015)) +
  geom_vline(aes(xintercept = median(nbPixels)),
             linetype = "dashed",
             linewidth = 0.6) +
  theme_bw() +
  xlab('Number of Pixels') +
  ylab('Density') +
  ggtitle('FRI') +
  theme(text = element_text(size = 20))

ggsave(
  str_c(out_loc, '/fri_final_pix_dens.png'),
  width = 2100,
  height = 2100,
  units = 'px'
)


# plot shape index
ggplot(data.frame(msi = as.numeric(p2$msi)), aes(x = msi)) +
  geom_density(fill = 'grey') +
  xlim(c(0, 10)) +
  ylim(c(0, 1.5)) +
  geom_vline(aes(xintercept = median(msi)),
             linetype = "dashed",
             linewidth = 0.6) +
  theme_bw() +
  xlab('Shape Index') +
  ylab('Density') +
  ggtitle('FRI') +
  theme(text = element_text(size = 20))

ggsave(
  str_c(out_loc, '/fri_final_shape_dens.png'),
  width = 2100,
  height = 2100,
  units = 'px'
)
