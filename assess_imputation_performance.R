# This code runs an imputation algorithm on LiDAR attributes connecting them to interpreter
# derived response attributes resulting in a model that can be applied to forest stands
# segmented from LiDAR rasters

# load packages
library(terra)
library(tidyverse)
library(yaImpute)
library(randomForest)

# example 1 from yaImpute user guide
data('iris')
set.seed(1)
refs <- sample(rownames(iris), 50)
x <- iris[, 1:3]
y <- iris[refs, 4:5]

ex <- yai(x = x, y = y, method = 'mahalanobis')
plot(ex)

head(impute(ex))

# example 2 from yaImpute user guide
data("MoscowMtStJoe")
x <- MoscowMtStJoe[, c("EASTING", "NORTHING", "ELEVMEAN",
                       "SLPMEAN", "ASPMEAN", "INTMEAN", "HTMEAN", "CCMEAN")]
x[, 5] <- (1 - cos((x[, 5] - 30) * pi/180))/2
names(x)[5] = "TrASP"
y <- MoscowMtStJoe[, c(1, 9, 12, 14, 18)]
mal <- yai(x = x, y = y, method = "mahalanobis")
msn <- yai(x = x, y = y, method = "msn")
gnn <- yai(x = x, y = y, method = "gnn")
ica <- yai(x = x, y = y, method = "ica")

y2 <- cbind(whatsMax(y[, 1:4]), y[, 5])
names(y2) <- c("MajorSpecies", "BasalAreaMajorSp", "TotalBA")
rf <- yai(x = x, y = y2, method = "randomForest")
head(y2)

plot(rf, vars = yvars(rf))

rfImp <- impute(rf)
rmsd <- compare.yai(mal, msn, gnn, rfImp, ica)
apply(rmsd, 2, mean, na.rm = TRUE)
plot(rmsd)

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

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# change all non-numeric variables to factor
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], 
                                       as.factor)

#load names of lidar columns added to data above
colnames(dat[,112:123])

# build dataframe of reference variables
ref_vars <- c('lor', 'ba', 'qmdbh', 'dens', 'agb', 'cc', 'cv')
x <- dat[, ref_vars]

# build dataframe of target variables
tar_vars <- c('HT', 'CC', 'BA', 'AGE', 'POLYTYPE')
y <- dat[, tar_vars]

# build dataframe of ancillary data
anci <- dat[, !(colnames(dat) %in% ref_vars)]

# run imputation algorithm
rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y))

# plot output
plot(rf, vars = yvars(rf))

# we want to create a system to compare the different outputs...
# R^2 values??

####################################################################
###LOAD LIDAR DERIVED POLYGONS AND POPULATE WITH LIDAR ATTRIBUTES###
####################################################################

# load lidar derived polygons
poly_lidar <- vect('D:/ontario_inventory/segmentation/ms_10_10_100_agg_na.shp')

# take a subset to test imputation
poly_lidar <- poly_lidar[1:1000,]

# convert to df
dat_lidar <- as.data.frame(poly_lidar)

# only keep label and nbPixels columns
dat_lidar <- dat_lidar %>% select(label, nbPixels)

# we need to set rownames of lidar polygons so they don't overlap reference data
# find final rowname from reference and add 1
start_r <- as.numeric(rownames(dat)[NROW(dat)]) + 1
end_r <- start_r + NROW(dat_lidar) - 1

# change row names so don't overlap with reference dataset
rownames(dat_lidar) <- start_r:end_r

# load LiDAR datasets we need to extract over polygons
# create named vector with variable names and data links

lidar <- c('cc' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif',
           'avg' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_avg.tif',
           'max' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_max.tif',
           'qav' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_qav.tif',
           'ske' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_ske.tif',
           'kur' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_kur.tif',
           'cv' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_cv.tif',
           'lor' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_lor.tif',
           'ba' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_ba_ha.tif',
           'qmdbh' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_qmdbh.tif',
           'dens' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_dens.tif',
           'agb' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_AGB_ha.tif')

# loop through LiDAR datasets and add to main dataframe
for(i in 1:length(lidar)){

  # load raster variable
  ras <- rast(lidar[i])

  # project poly to crs of raster
  poly_ras <- project(poly_lidar, ras)

  # extract median values within each polygon
  ras_med <- terra::extract(ras, poly_ras, fun = function(x){median(x, na.rm = T)})

  # add new column into dat
  dat_lidar <- dat_lidar %>% add_column(ras = ras_med[,2])

  # change column name
  colnames(dat_lidar)[NCOL(dat_lidar)] <- names(lidar[i])

  # clean up
  rm(ras, poly_ras, ras_med)
}

# clean up
rm(i, lidar)

# run imputation over new polygons
imp_lidar <- predict(object = rf,
                     newdata = dat_lidar,
                     ancillaryData = anci,
                     observed = F)

# remove columns from imputed dataset if they exist from lidar
imp_lidar <- imp_lidar[, !(colnames(imp_lidar) %in% colnames(dat_lidar))]

# add additional columns from lidar dataset that were not imputed
imp_lidar <- imp_lidar %>% add_column(dat_lidar)

# repopulate polygons with data
values(poly_lidar) <- as.data.frame(imp_lidar)
