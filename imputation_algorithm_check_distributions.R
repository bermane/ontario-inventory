# This code runs an imputation algorithm on LiDAR attributes connecting them to interpreter
# derived response attributes resulting in a model that can be applied to forest stands
# segmented from LiDAR rasters

# The purpose is to compare the distributions of various attributes (HT, CC, BA, AGE, POLYTYPE, WG)
# from FRI polygons, LiDAR attributes in FRI polygons, imputed vars into FRI polys, imputed vars into LiDAR polys

# load packages
library(terra)
library(tidyverse)
library(yaImpute)
library(randomForest)

###################################################
###LOAD POLYGON DATASET AND ADD LIDAR ATTRIBUTES###
###################################################

# # load photo interpreted polygons
# poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')
#
# # convert to df
# dat <- as.data.frame(poly)

# load LiDAR datasets we need to extract over polygons
# create named vector with variable names and data links

# lidar <- c('cc' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif',
#            'avg' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_avg.tif',
#            'max' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_max.tif',
#            'qav' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_qav.tif',
#            'ske' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_ske.tif',
#            'kur' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_kur.tif',
#            'cv' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_cv.tif',
#            'lor' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_lor.tif',
#            'ba' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_ba_ha.tif',
#            'qmdbh' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_qmdbh.tif',
#            'dens' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_dens.tif',
#            'agb' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_AGB_ha.tif')
#
# # loop through LiDAR datasets and add to main dataframe
# for(i in 1:length(lidar)){
#
#   # load raster variable
#   ras <- rast(lidar[i])
#
#   # project poly to crs of raster
#   poly_ras <- project(poly, ras)
#
#   # extract median values within each polygon
#   ras_med <- terra::extract(ras, poly_ras, fun = function(x){median(x, na.rm = T)})
#
#   # add new column into dat
#   dat <- dat %>% add_column(ras = ras_med[,2])
#
#   # change column name
#   colnames(dat)[NCOL(dat)] <- names(lidar[i])
#
#   # clean up
#   rm(ras, poly_ras, ras_med)
# }
#
# # clean up
# rm(i, lidar)
# 
# # save extracted dataframe for fast rebooting
# save.image('D:/ontario_inventory/imputation/imputation_df.RData')

# load extracted data frame
load('D:/ontario_inventory/imputation/imputation_df.RData')

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

# save dat into different object
dat_full <- dat

# change missing WG values to "UCL" so they aren't blank
dat_full$WG[dat_full$WG == ""] <- "UCL"

# create further separation of data to impute
# remove WAT POLYTYPE
dat_imp <- dat_full[dat_full$POLYTYPE != "WAT",]

# change all non-numeric variables to factor
dat_imp[sapply(dat_imp, is.character)] <- lapply(dat_imp[sapply(dat_imp, is.character)], 
                                                   as.factor)

# create ID col based on row names
# dat_full$ID <- rownames(dat_full)

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen_2_30_70.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# SPL vars only

# change missing WG values to "UCL" so they aren't blank
dat$WG[dat$WG == ""] <- "UCL"

# remove WAT POLYTYPE
# already did though in screening
# dat <- dat[dat$POLYTYPE != "WAT",]

# change all non-numeric variables to factor
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], 
                                       as.factor)

#load names of lidar columns added to data above
colnames(dat[,112:123])

# build dataframe of reference variables
ref_vars <- c('cc', 'p80', 'avg', 'qav', 'cv', 'kur', 'max', 'ske')
x <- dat[, ref_vars]

# build dataframe of target variables
tar_vars <- c('HT', 'CC', 'BA', 'POLYTYPE', 'WG')
y <- dat[, tar_vars]

# build dataframe of ancillary data
anci <- dat[, !(colnames(dat) %in% ref_vars)]

# run imputation algorithm
rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y))

# plot output
plot(rf, vars = yvars(rf))

####################################################################
###LOAD LIDAR DERIVED POLYGONS AND POPULATE WITH LIDAR ATTRIBUTES###
####################################################################

# load lidar derived polygons
poly_lidar <- vect('D:/ontario_inventory/segmentation/ms_10_10_100_agg_na.shp')
# 
# # currently testing run on all polygons
# # take a subset to test imputation
# # poly_lidar <- poly_lidar[1:1000,]
# 
# # convert to df
# dat_lidar <- as.data.frame(poly_lidar)
# 
# # only keep label and nbPixels columns
# dat_lidar <- dat_lidar %>% select(label, nbPixels)
# 
# # we need to set rownames of lidar polygons so they don't overlap reference data
# # find final rowname from reference and add 1
# start_r <- as.numeric(rownames(dat)[NROW(dat)]) + 1
# end_r <- start_r + NROW(dat_lidar) - 1
# 
# # change row names so don't overlap with reference dataset
# rownames(dat_lidar) <- start_r:end_r
# 
# # load LiDAR datasets we need to extract over polygons
# # create named vector with variable names and data links
# 
# lidar <- c('cc' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif',
#            'avg' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_avg.tif',
#            'max' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_max.tif',
#            'p80' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_p80.tif',
#            'qav' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_qav.tif',
#            'ske' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_ske.tif',
#            'kur' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_kur.tif',
#            'cv' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_cv.tif',
#            'lor' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_lor.tif',
#            'ba' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_ba_ha.tif',
#            'qmdbh' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_qmdbh.tif',
#            'dens' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_dens.tif',
#            'agb' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_AGB_ha.tif',
#            'top_height' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_top_height.tif',
#            'v' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_V_ha.tif',
#            'v_merch' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_Vmerch_ha.tif')
# 
# # loop through LiDAR datasets and add to main dataframe
# for(i in 1:length(lidar)){
# 
#   # load raster variable
#   ras <- rast(lidar[i])
# 
#   # project poly to crs of raster
#   poly_ras <- project(poly_lidar, ras)
# 
#   # extract median values within each polygon
#   ras_med <- terra::extract(ras, poly_ras, fun = function(x){median(x, na.rm = T)})
# 
#   # add new column into dat
#   dat_lidar <- dat_lidar %>% add_column(ras = ras_med[,2])
# 
#   # change column name
#   colnames(dat_lidar)[NCOL(dat_lidar)] <- names(lidar[i])
# 
#   # clean up
#   rm(ras, poly_ras, ras_med)
# }
# 
# # clean up
# rm(i, lidar)
# 
# # save loaded LiDAR attributes
# save(dat_lidar, file = 'D:/ontario_inventory/imputation/seg_df_ms_10_10_100_agg_na.RData')

# load LiDAR attributes
load('D:/ontario_inventory/imputation/seg_df_ms_10_10_100_agg_na.RData')

# run imputation over new polygons
imp_lidar <- predict(object = rf,
                     newdata = dat_lidar,
                     ancillaryData = anci,
                     observed = F)

# find any rows that were removed
rm_rows <- setdiff(rownames(dat_lidar), rownames(imp_lidar))

# re add any rows that were removed
if(length(rm_rows) != 0){
  add_rows <- imp_lidar[1:length(rm_rows),]
  add_rows[] <- NA
  rownames(add_rows) <- rm_rows
  imp_lidar <- rbind(imp_lidar, add_rows)
}

# re order the rows of imp_lidar to match dat_lidar
imp_lidar <- imp_lidar[order(as.numeric(row.names(imp_lidar))),]

# remove columns from imputed dataset if they exist from lidar
imp_lidar <- imp_lidar[, !(colnames(imp_lidar) %in% colnames(dat_lidar))]

# add additional columns from lidar dataset that were not imputed
imp_lidar <- imp_lidar %>% add_column(dat_lidar)

# change back to original row names
rownames(imp_lidar) <- 1:NROW(imp_lidar)

# repopulate polygons with data
values(poly_lidar) <- as.data.frame(imp_lidar)

# do we need to re add big polygon of NAs?

# write polygons with imputed data
writeVector(poly_lidar, filename = 'D:/ontario_inventory/imputation/distributions/vector/ms_10_10_100_agg_na_IMP_SPL.shp')

# write dataframe with imputed data
write.csv(imp_lidar, file = 'D:/ontario_inventory/imputation/distributions/vector/ms_10_10_100_agg_na_IMP_SPL.csv')

#################################
###COMPARISON WITH FRI DATASET###
#################################

# we have entire raw FRI dataset (with lidar attributes)
# also need to apply the imputation to the entire dataset of FRI polygons

# build dataframe of ancillary data we actually want for the comparison
anci <- dat[, c('HT', 'CC', 'BA', 'AGE', 'POLYID', 'POLYTYPE', 'WG')]

# run imputation over new polygons
# we want to keep observed cols to see relationship
imp_full <- predict(object = rf,
                    newdata = dat_imp,
                    ancillaryData = anci,
                    observed = F)

# add observed cols
imp_full <- imp_full %>% add_column(HT_o = dat_full$HT[rownames(dat_full) %in% rownames(imp_full)],
                                    CC_o = dat_full$CC[rownames(dat_full) %in% rownames(imp_full)],
                                    BA_o = dat_full$BA[rownames(dat_full) %in% rownames(imp_full)],
                                    AGE_o = dat_full$AGE[rownames(dat_full) %in% rownames(imp_full)],
                                    POLYID_o = dat_full$POLYID[rownames(dat_full) %in% rownames(imp_full)],
                                    POLYTYPE_o = dat_full$POLYTYPE[rownames(dat_full) %in% rownames(imp_full)],
                                    WG_o = dat_full$WG[rownames(dat_full) %in% rownames(imp_full)])

# to compare the full distribution of the dataset we also want the imputed values from the
# data actually used to build the model
imp_screen <- impute(object = rf,
                  ancillaryData = anci,
                  observed = T)

# rename observed cols
imp_screen <- rename_with(imp_screen, ~ gsub('.o', '_o', .x, fixed = T))

# rbind imputed data
imp_full <- rbind(imp_full, imp_screen)

# sort by row name
imp_full <- imp_full[order(as.numeric(row.names(imp_full))),]

# find any rows that were removed during imputation
rm_rows <- dat_full[rownames(dat_full) %in% setdiff(rownames(dat_full), rownames(imp_full)),]

# format to match imp_full df
rm_rows <- data.frame(HT = NA,
                      CC = NA,
                      BA = NA,
                      AGE = NA,
                      POLYID = NA,
                      POLYTYPE = NA,
                      WG = NA,
                      HT_o = rm_rows$HT,
                      CC_o = rm_rows$CC,
                      BA_o = rm_rows$BA,
                      AGE_o = rm_rows$AGE,
                      POLYID_o = rm_rows$POLYID,
                      POLYTYPE_o = rm_rows$POLYTYPE,
                      WG_o = rm_rows$WG)

rownames(rm_rows) <- rm_rows$POLYID_o

# rbind to imp_full
imp_full <- rbind(imp_full, rm_rows)

# sort by row name
imp_full <- imp_full[order(as.numeric(row.names(imp_full))),]

# add column of original lidar attributes we are looking at
imp_full <- imp_full %>% add_column(cc = dat_full$cc[dat_full$POLYID %in% imp_full$POLYID_o],
                                    ba = dat_full$ba[dat_full$POLYID %in% imp_full$POLYID_o],
                                    lor = dat_full$lor[dat_full$POLYID %in% imp_full$POLYID_o])

# create new polygon dataset for imp data
poly_imp <- poly

# fill with data
values(poly_imp) <- as.data.frame(imp_full)

# write polygons with imputed data
writeVector(poly_imp, filename = 'D:/ontario_inventory/imputation/distributions/vector/fri_polygons_imputed.shp')

# write dataframe with imputed data
write.csv(imp_full, file = 'D:/ontario_inventory/imputation/distributions/vector/fri_polygons_imputed.csv')

#########################################
###COMPARE DISTRIBUTIONS OF ATTRIBUTES###
#########################################

#need to choose a color palette colorblind friendly
#https://personal.sron.nl/~pault/#sec:qualitative
#'#4477AA', '#000000', '#228833', '#CCBB44', '#EE6677', '#AA3377'

colscale_fri <- scale_color_manual(name = "Dataset", values = c('FRI Polygons' = '#4477AA', 'LiDAR Imp Polygons' = '#EE6677',
                                                                'FRI Imp Polygons' = '#228833'))

# set line size
sz <- 1.2

# Height
ggplot() +
  geom_density(data = imp_full, mapping = aes(x = lor, color = 'FRI (Lor Height)'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = lor, color = 'LiDAR (Lor Height)'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = HT_o, color = 'FRI (HT)'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = HT, color = 'LiDAR Imp (HT)'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = HT, color = 'FRI Imp (HT)'), size = sz) +
  theme_bw() + 
  ggtitle("Height Density Plot") +
  xlab("Height (m)") + ylab("Density") + 
  scale_color_manual(name = "Dataset", values = c('FRI (Lor Height)' = '#CCBB44', 'LiDAR (Lor Height)' = '#AA3377',
                                                  'FRI (HT)' = '#4477AA', 'LiDAR Imp (HT)' = '#EE6677',
                                                  'FRI Imp (HT)' = '#228833'))
ggsave(filename = 'D:/ontario_inventory/imputation/distributions/height.png')

# Canopy Cover
ggplot() +
  geom_density(data = imp_full, mapping = aes(x = cc, color = 'FRI (% 2m Rtn)'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = cc, color = 'LiDAR (% 2m Rtn)'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = CC_o, color = 'FRI (CC)'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = CC, color = 'LiDAR Imp (CC)'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = CC, color = 'FRI Imp (CC)'), size = sz) +
  theme_bw() + 
  ggtitle("Canopy Cover Density Plot") +
  xlab("Cover (%)") + ylab("Density") + 
  scale_color_manual(name = "Dataset", values = c('FRI (% 2m Rtn)' = '#CCBB44', 'LiDAR (% 2m Rtn)' = '#AA3377',
                                                  'FRI (CC)' = '#4477AA', 'LiDAR Imp (CC)' = '#EE6677',
                                                  'FRI Imp (CC)' = '#228833'))
ggsave(filename = 'D:/ontario_inventory/imputation/distributions/canopy_cover.png')

# Basal Area
ggplot() +
  geom_density(data = imp_full, mapping = aes(x = ba, color = 'FRI (EFI gen)'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = ba, color = 'LiDAR (EFI gen)'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = BA_o, color = 'FRI (BA)'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = BA, color = 'LiDAR Imp (BA)'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = BA, color = 'FRI Imp (BA)'), size = sz) +
  theme_bw() + 
  ggtitle("Basal Area Density Plot") +
  xlab(expression("Basal Area (m"^2*"/ha)")) + ylab("Density") + 
  scale_color_manual(name = "Dataset", values = c('FRI (EFI gen)' = '#CCBB44', 'LiDAR (EFI gen)' = '#AA3377',
                                                  'FRI (BA)' = '#4477AA', 'LiDAR Imp (BA)' = '#EE6677',
                                                  'FRI Imp (BA)' = '#228833'))
ggsave(filename = 'D:/ontario_inventory/imputation/distributions/basal_area.png')

# Age
ggplot() +
  geom_density(data = imp_full, mapping = aes(x = AGE_o, color = 'FRI Polygons'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = AGE, color = 'LiDAR Imp Polygons'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = AGE, color = 'FRI Imp Polygons'), size = sz) +
  theme_bw() + 
  ggtitle("Interpreter Derived Age Density Plot") +
  xlab("Age (years)") + ylab("Density") + 
  colscale_fri
ggsave(filename = 'D:/ontario_inventory/imputation/distributions/age.png')

# Polytype
# better as freq table
# load full data categories
polyt_freq <- plyr::count(imp_full, 'POLYTYPE_o')

# change col names
polyt_freq <- rename(polyt_freq, POLYTYPE = POLYTYPE_o, fri_polygons = freq)

# remove WAT polytype since didn't use for imputation
polyt_freq <- filter(polyt_freq, POLYTYPE != "WAT")

# add other polys
add1 <- plyr::count(imp_lidar, 'POLYTYPE')
add1 <- drop_na(add1)
add2 <- plyr::count(imp_full, 'POLYTYPE')
add2 <- drop_na(add2)

polyt_freq$add1 <- 0
polyt_freq$add2 <- 0

polyt_freq$add1[polyt_freq$POLYTYPE %in% add1$POLYTYPE] <- add1$freq
polyt_freq$add2[polyt_freq$POLYTYPE %in% add2$POLYTYPE] <- add2$freq

# change col names
polyt_freq <- rename(polyt_freq, lidar_imp_polygons = add1, fri_imp_polygons = add2)

# add percentages
polyt_freq <- add_column(polyt_freq, fri_poly_perc = polyt_freq$fri_polygons/sum(polyt_freq$fri_polygons)*100,
                         lidar_imp_poly_perc = polyt_freq$lidar_imp_polygons/sum(polyt_freq$lidar_imp_polygons)*100,
                         fri_imp_poly_perc = polyt_freq$fri_imp_polygons/sum(polyt_freq$fri_imp_polygons)*100)

# round
polyt_freq[,(NCOL(polyt_freq)-2):NCOL(polyt_freq)] <- round(polyt_freq[,(NCOL(polyt_freq)-2):NCOL(polyt_freq)], 2)

# export as csv
write.csv(polyt_freq, 
          file = 'D:/ontario_inventory/imputation/distributions/polytype.csv', 
          row.names = F)

# Working group
# better as freq table
# load full data categories
# again don't include water polygons since skews the distribution
wg_freq <- plyr::count(imp_full %>% filter(POLYTYPE_o !='WAT'), 'WG_o')

# change col names
wg_freq <- rename(wg_freq, WG = WG_o, fri_polygons = freq)

# add other polys
add1 <- plyr::count(imp_lidar, 'WG')
add1 <- drop_na(add1)
add2 <- plyr::count(imp_full, 'WG')
add2 <- drop_na(add2)

wg_freq$add1 <- 0
wg_freq$add2 <- 0

wg_freq$add1[wg_freq$WG %in% add1$WG] <- add1$freq
wg_freq$add2[wg_freq$WG %in% add2$WG] <- add2$freq

# change col names
wg_freq <- rename(wg_freq, lidar_imp_polygons = add1, fri_imp_polygons = add2)

# add percentages
wg_freq <- add_column(wg_freq, fri_poly_perc = wg_freq$fri_polygons/sum(wg_freq$fri_polygons)*100,
                         lidar_imp_poly_perc = wg_freq$lidar_imp_polygons/sum(wg_freq$lidar_imp_polygons)*100,
                         fri_imp_poly_perc = wg_freq$fri_imp_polygons/sum(wg_freq$fri_imp_polygons)*100)

# round
wg_freq[,(NCOL(wg_freq)-2):NCOL(wg_freq)] <- round(wg_freq[,(NCOL(wg_freq)-2):NCOL(wg_freq)], 2)

# export as csv
write.csv(wg_freq, 
          file = 'D:/ontario_inventory/imputation/distributions/working_group.csv', 
          row.names = F)
