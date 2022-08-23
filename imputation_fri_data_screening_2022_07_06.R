# This code loads photo interpreted polygons from the Romeo Malette Forest in Ontario
# and applies a number of data screens to curate an accurate and reliable dataset
# this version removes all polygons that are not listed as productive forest in the FRI
# because we want to calibrate the model to accurately predict forest structure and 
# impute forest variables. Including other land cover types decreases the performance
# of the prediction model

# load packages
library(terra)
library(tidyverse)
library(exactextractr)
library(sf)
library(magrittr)

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

# convert to df
dat <- as.data.frame(poly)

###########################
###DATA SCREENING PART 0###
###########################

# remove non forested polygons
dat <- filter(dat, POLYTYPE == 'FOR')

# create smaller polygon set without water and unclassified polytypes
poly_dat <- poly[poly$POLYTYPE == 'FOR']

# export forested polygons only for maps
# writeVector(poly_dat, filename = 'D:/ontario_inventory/imputation/distributions_for_only/vector/fri_polygons_for_only.shp', overwrite = T)

###########################
###DATA SCREENING PART 1###
###########################

# A. Polygon area < 50 ha (500000 m^2)

# check number of polygons with area < 50 ha
#NROW(dat[dat$AREA < 500000,])

# subset data
#dat_1a <- dat[dat$AREA < 500000,]

# look AT HISTOGRAM OF VERY LARGE POLYGONS AND GET RID OF THEM and ONLY FOR POLYGONS!!!

# let's skip this skip for now.

# B. Within-polygon height and canopy cover coefficient of variation < 0.5

# load canopy cover and p95 from LiDAR
lidar <- c('cc' ='D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif',
               'p95' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_p95.tif')

# load rasters
lidar_ras <- rast(lidar)

# project poly to crs of rasters
poly_ras <- project(poly_dat, lidar_ras)

# convert to sf
poly_ras <- st_as_sf(poly_ras)

# extract coeff of variation -- no edge pixels
#extract median values
vec <- exact_extract(lidar_ras, poly_ras, function(values, coverage_fraction){
  values <- values[coverage_fraction == 1,]
  apply(values, 2, function(x) {sd(x, na.rm = T)/mean(x, na.rm = T)})
})

# transpose matrix and make data.frame
vec <- t(vec) %>% as.data.frame

# change column names
colnames(vec) <- c('cc_cv', 'p95_cv')

# add new column into dat
dat <- cbind(dat, vec)

# check number of polygons with cv < 0.5
NROW(dat[dat$cc_cv < 0.5,])
NROW(dat[dat$p95_cv < 0.5,])

# subset dat based on cv
dat_1b <- dat[dat$cc_cv < 0.5 & dat$p95_cv < 0.5,]

# C. number of land cover types within a polygon < 2

# load VLCE 2.0 landcover dataset from 2018
lc <- rast('D:/ontario_inventory/VLCE/CA_forest_VLCE2_2018_CLIPPED.tif')

# project poly to crs of raster
poly_lc <- project(poly_dat, lc)

# convert to sf
poly_lcsf <- st_as_sf(poly_lc)

# extract landcover values
lc_poly <- exact_extract(lc, poly_lcsf)

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

# find number of unique lc types in each polygon
# apply over list
lc_uni <- sapply(lc_poly, function(x){
  x$value <- recode(x$value, !!!lc_key)
  x %<>% filter(coverage_fraction >= 0.5)
  return(length(unique(x$value)))
})

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
lc_dom_for <- sapply(lc_poly, function(x){
  x$value <- recode(x$value, !!!lc_key_for)
  x <- x %>% group_by(value) %>% summarize(sum = sum(coverage_fraction))
  m <- x$value[which(x$sum == max(x$sum))]
  if((length(m) == 1) & (m == 'Forest')[1]){
    if(x$sum[x$value == m]/sum(x$sum) >= 0.5){
      return('Yes')
    }else{return('No')}
  }else{return('No')}
})

# add new columns into dat
dat <- dat %>% add_column(num_uniq = lc_uni,
                          dom_for = lc_dom_for)

# rm variables
rm(lc, lc_poly, poly_lc, lc_key, lc_key_for,
   lc_dom_for, lc_uni, poly_lcsf)

# check number of polygons with LC types < 2
# perhaps Chen meant LC types <= 2??
NROW(dat[dat$num_uniq == 1,])
NROW(dat[dat$num_uniq %in% c(1, 2),])

# # check number of polygons with forested classes > 95%
# NROW(dat[dat$perc_for >= .95,])

# try dominant forest approach
# since we will be imputing into polygons that are dominant forested
NROW(dat[dat$dom_for == 'Yes',])

# subset dat based on lc info
dat_1c <- dat[dat$num_uniq %in% c(1,2) & dat$dom_for == 'Yes',]

###########################
###DATA SCREENING PART 2###
###########################

# first load all data, take median values, and enter into dat
# Load Lorey's Height, cc, and ba from LiDAR
lidar <- c('cc' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif',
           'ba' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_ba_ha.tif',
           'p95' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_p95.tif')

# load LiDAR rasters as raster stack
lidar_ras <- rast(lidar)

# project poly to crs of raster
poly_ras <- project(poly_dat, lidar_ras)

# convert to sf
poly_ras <- st_as_sf(poly_ras)

#extract median values
vec <- exact_extract(lidar_ras, poly_ras, function(values, coverage_fraction){
  values <- values[coverage_fraction == 1,]
  apply(values, 2, function(x) median(x, na.rm = T))
})

# transpose matrix and make data.frame
vec <- t(vec) %>% as.data.frame

# change column names
colnames(vec) <- names(lidar)

# add new column into dat
dat <- cbind(dat, vec)

# How many rows have p95 > 5, which is definition of forest
NROW(dat[dat$p95 >= 5,])

# require p95 > 5
dat_2_p95 <- dat[dat$p95 >= 5,]

# How many rows have cc >= 10, which is definition of forest
NROW(dat[dat$cc >= 10,])

# require cc > 10
dat_2_cc <- dat[dat$cc >= 10,]

# A. Linear regression of Height vs. P95
# we want to exclude using EFI variables so more applicable anywhere

# set the percentile of residuals we want to get rid of
perc_resid <- c(.10, .90)

# remove rows with missing values
dat_2a <- dat[is.na(dat$HT) == F & is.na(dat$p95) == F,]

# plot variables against each other
plot(dat_2a$p95, dat_2a$HT, xlab = 'p95', ylab = 'HT', 
     main = 'HT vs. p95 in FRI Polygons')

# run simple linear model
lm_ht <- lm(HT ~ p95, data = dat_2a)
summary(lm_ht)

# bind residuals to data
dat_2a <- dat_2a %>% add_column(lm_ht_resid = resid(lm_ht))

# find 10th and 90th percentile
lm_ht_perc <- quantile(dat_2a$lm_ht_resid, probs = perc_resid)

# remove lower 10 and upper 10 from dat
dat_2a <- dat_2a[dat_2a$lm_ht_resid > lm_ht_perc[1] & dat_2a$lm_ht_resid < lm_ht_perc[2],]

# re plot relationship
plot(dat_2a$p95, dat_2a$HT, xlab = 'p95', ylab = 'HT', 
     main = 'HT vs. p95 in FRI Polygons 10/90 Residual Screened')

# run simple linear model on remaining data
lm_ht2 <- lm(HT ~ p95, data = dat_2a)
summary(lm_ht2)

# # B. Linear regression of Canopy Closure vs. LiDAR Canopy Cover
# 
# # remove rows with missing values
# dat_2b <- dat[is.na(dat$CC) == F & is.na(dat$cc) == F,]
# 
# # plot variables against each other
# plot(dat_2b$cc, dat_2b$CC)
# 
# # run simple linear model
# lm_cc <- lm(CC ~ cc, data = dat_2b)
# summary(lm_cc)
# 
# # bind residuals to data
# dat_2b <- dat_2b %>% add_column(lm_cc_resid = resid(lm_cc))
# 
# # find 10th and 90th percentile
# lm_cc_perc <- quantile(dat_2b$lm_cc_resid, probs = perc_resid)
# 
# # remove lower 10 and upper 10 from dat
# dat_2b <- dat_2b[dat_2b$lm_cc_resid > lm_cc_perc[1] & dat_2b$lm_cc_resid < lm_cc_perc[2],]
# 
# # re plot relationship
# plot(dat_2b$cc, dat_2b$CC)
# 
# # run simple linear model on remaining data
# lm_cc2 <- lm(CC ~ cc, data = dat_2b)
# summary(lm_cc2)
# 
# # C. Linear regression of Basal Area vs. LiDAR Basal Area
# 
# # remove rows with missing values
# dat_2c <- dat[is.na(dat$BA) == F & is.na(dat$ba) == F,]
# 
# # plot variables against each other
# plot(dat_2c$ba, dat_2c$BA)
# 
# # run simple linear model
# lm_ba <- lm(BA ~ ba, data = dat_2c)
# summary(lm_ba)
# 
# # bind residuals to data
# dat_2c <- dat_2c %>% add_column(lm_ba_resid = resid(lm_ba))
# 
# # find 10th and 10th percentile
# lm_ba_perc <- quantile(dat_2c$lm_ba_resid, probs = perc_resid)
# 
# # remove lower 10 and upper 10 from dat
# dat_2c <- dat_2c[dat_2c$lm_ba_resid > lm_ba_perc[1] & dat_2c$lm_ba_resid < lm_ba_perc[2],]
# 
# # re plot relationship
# plot(dat_2c$ba, dat_2c$BA)
# 
# # run simple linear model on remaining data
# lm_ba2 <- lm(BA ~ ba, data = dat_2c)
# summary(lm_ba2)

####################################
### COMBINE INTERSECTION OF DATA ###
####################################

# combine results from part 1
dat_1 <- intersect(dat_1b %>% subset(select = POLYID), dat_1c %>% subset(select = POLYID))

# combine results from part 2 (POLYID column only)
# dat_2 <- intersect(dat_2a %>% subset(select = POLYID), dat_2b %>% subset(select = POLYID)) %>%
#   intersect(., dat_2c %>% subset(select = POLYID)) %>%
#   intersect(., dat_2_pc %>% subset(select = POLYID))

# only use results from part 2a since the other models are not great
dat_2 <- intersect(dat_2_p95 %>% subset(select = POLYID), dat_2_cc %>% subset(select = POLYID))
dat_2 <- intersect(dat_2, dat_2a %>% subset(select = POLYID))

# combine all
dat_screen <- intersect(dat_1, dat_2)

# reload attributes
dat_screen <- dat[dat$POLYID %in% dat_screen$POLYID,]

# write to disk
write.csv(dat_screen, file = 'D:/ontario_inventory/imputation/fri_data_screen_1bc_2a_2p95_2cc_10perc_2022_07_06.csv', row.names = F)

# save working image to speed up future changes
# save.image('D:/ontario_inventory/imputation/imputation_fri_data_screening_2022_07_06.RData')

# # load image
# load('D:/ontario_inventory/imputation/imputation_fri_data_screening_2022_07_06.RData')
# 
# # check lm of ht and p95 with screened data
# 
# # plot variables against each other
# plot(dat_screen$p95, dat_screen$HT, xlab = 'p95', ylab = 'HT', 
#      main = 'HT vs. p95 in FRI Polygons')
# 
# # run simple linear model
# lm_ht <- lm(HT ~ p95, data = dat_screen)
# summary(lm_ht)
# 
# w# Make predictions
# predictions <- lm_ht %>% predict(dat_screen)
# 
# # calculate rmse different variables
# rmse <- function(obs, est) sqrt(mean((obs - est)^2))
# 
# # Model performance
# # (a) Prediction error, RMSE
# rmse(dat_screen$HT, predictions)
# 
# # (b) R-square
# R2(predictions, test.data$sales)

