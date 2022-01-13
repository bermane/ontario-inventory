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
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen_2_30_70.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# change missing WG values to "UCL" so they aren't blank
dat$WG[dat$WG == ""] <- "UCL"

# change all non-numeric variables to factor
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], 
                                       as.factor)

#load names of lidar columns added to data above
colnames(dat[,112:123])

# build dataframe of reference variables
ref_vars <- c('cc', 'p80', 'avg', 'qav', 'cv', 'kur', 'max', 'ske', 'agb', 'ba', 'dens', 'lor', 'qmdbh', 'top_height', 'v', 'v_merch')
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

# # load lidar derived polygons
# poly_lidar <- vect('D:/ontario_inventory/segmentation/ms_10_10_100_agg_na.shp')
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

# save loaded LiDAR attributes
#save(dat_lidar, file = 'D:/ontario_inventory/imputation/seg_df_ms_10_10_100_agg_na.RData')

# load lidar attributes
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

# remove columns from imputed dataset if they exist from lidar
imp_lidar <- imp_lidar[, !(colnames(imp_lidar) %in% colnames(dat_lidar))]

# add additional columns from lidar dataset that were not imputed
imp_lidar <- imp_lidar %>% add_column(dat_lidar)

# repopulate polygons with data
values(poly_lidar) <- as.data.frame(imp_lidar)

# do we need to re add big polygon of NAs?

# write polygons with imputed data
writeVector(poly_lidar, filename = 'D:/ontario_inventory/imputation/ms_10_10_100_agg_na_IMP.shp')

#################################
###COMPARISON WITH FRI DATASET###
#################################

# lets also apply the imputation to the entire dataset of FRI polygons
# since that will give a good indication of how some of the other attributes are comparing

# run imputation over new polygons
imp_fri <- impute(object = rf,
                  ancillaryData = anci,
                  observed = F)

# save dat and imp_fri
write.csv(dat, file = 'D:/ontario_inventory/imputation/test_imp_wg/dat_out.csv', row.names = F)
write.csv(imp_fri, file = 'D:/ontario_inventory/imputation/test_imp_wg/dat_out_imp.csv', row.names = F)

# compare input data and imputed data plots
plot(dat$HT, imp_fri$HT)
plot(dat$CC, imp_fri$CC)
plot(dat$BA, imp_fri$BA)
plot(dat$AGE, imp_fri$AGE)

# create df of POLYTYPE
polyt <- data.frame(obs = dat$POLYTYPE,
                    est = imp_fri$POLYTYPE)

# create column of match or not
polyt$match <- polyt$obs == polyt$est

# total percent of matching POLYTYPE
NROW(polyt[polyt$match == T,])/NROW(polyt)

# build accuracy table
accmat_pt <- table("pred" = polyt$est, "ref" = polyt$obs)

# UA
ua_pt <- diag(accmat_pt) / rowSums(accmat_pt) * 100

# PA
pa_pt <- diag(accmat_pt) / colSums(accmat_pt) * 100

# OA
oa_pt <- sum(diag(accmat_pt)) / sum(accmat_pt) * 100

# build confusion matrix
accmat_pt_ext <- addmargins(accmat_pt)
accmat_pt_ext <- rbind(accmat_pt_ext, "Users" = c(pa_pt, NA))
accmat_pt_ext <- cbind(accmat_pt_ext, "Producers" = c(ua_pt, NA, oa_pt))
#colnames(accmat_pt_ext) <- c(levels(as.factor(shp.train$classes)), "Sum", "PA")
#rownames(accmat_pt_ext) <- c(levels(as.factor(shp.train$classes)), "Sum", "UA")
accmat_pt_ext <- round(accmat_pt_ext, digits = 1)
dimnames(accmat_pt_ext) <- list("Prediction" = colnames(accmat_pt_ext),
                                "Reference" = rownames(accmat_pt_ext))
class(accmat_pt_ext) <- "table"
accmat_pt_ext

# write to csv
write.csv(accmat_pt_ext, file = 'D:/ontario_inventory/imputation/test_imp_wg/accmat_pt.csv')

# calculate kappa coefficient
kappa <- function(m) {
  N <- sum(m)
  No <- sum(diag(m))
  Ne <- 1 / N * sum(colSums(m) * rowSums(m))
  return( (No - Ne) / (N - Ne) )
}

kappa(accmat_pt)

# calculate rmse different variables
rmse <- function(obs, est) sqrt(mean((obs - est)^2))

rmse_df <- data.frame(AGE = c(rmse(dat$AGE, imp_fri$AGE),
                              rmse(dat$AGE, imp_fri$AGE)/sd(dat$AGE)),
                      HT = c(rmse(dat$HT, imp_fri$HT),
                             rmse(dat$HT, imp_fri$HT)/sd(dat$HT)),
                      BA = c(rmse(dat$BA, imp_fri$BA),
                             rmse(dat$BA, imp_fri$BA)/sd(dat$BA)),
                      CC = c(rmse(dat$CC, imp_fri$CC),
                             rmse(dat$CC, imp_fri$CC)/sd(dat$CC)))

rownames(rmse_df) <- c('raw_rmsd', 'scaled_rmsd')

# write table
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_wg/rmsd_dat_imp.csv')

# create df of WG
wg <- data.frame(obs = dat$WG,
                 est = imp_fri$WG)

# create column of match or not
wg$match <- wg$obs == wg$est

# total percent of matching WG
NROW(wg[wg$match == T,])/NROW(wg)

# check count of different working groups
plyr::count(wg, 'obs')
plyr::count(wg, 'est')

# build accuracy table
accmat_wg <- table("pred" = wg$est, "ref" = wg$obs)

# UA
ua_wg <- diag(accmat_wg) / rowSums(accmat_wg) * 100

# PA
pa_wg <- diag(accmat_wg) / colSums(accmat_wg) * 100

# OA
oa_wg <- sum(diag(accmat_wg)) / sum(accmat_wg) * 100

# build confusion matrix
accmat_wg_ext <- addmargins(accmat_wg)
accmat_wg_ext <- rbind(accmat_wg_ext, "Users" = c(pa_wg, NA))
accmat_wg_ext <- cbind(accmat_wg_ext, "Producers" = c(ua_wg, NA, oa_wg))
#colnames(accmat_wg_ext) <- c(levels(as.factor(shp.train$classes)), "Sum", "PA")
#rownames(accmat_wg_ext) <- c(levels(as.factor(shp.train$classes)), "Sum", "UA")
accmat_wg_ext <- round(accmat_wg_ext, digits = 1)
dimnames(accmat_wg_ext) <- list("Prediction" = colnames(accmat_wg_ext),
                             "Reference" = rownames(accmat_wg_ext))
class(accmat_wg_ext) <- "table"
accmat_wg_ext

# write to csv
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/test_imp_wg/accmat_wg.csv')

# calculate kappa coefficient
kappa(accmat_wg)
