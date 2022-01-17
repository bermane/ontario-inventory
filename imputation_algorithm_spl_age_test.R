# This code runs an imputation algorithm on LiDAR attributes connecting them to interpreter
# derived response attributes resulting in a model that can be applied to forest stands
# segmented from LiDAR rasters

# The model is applied to a model that now includes age to see if the results are any different from
# when age is not included in the model (see imputation_algorithm_spl_efi_test.R)

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

# save dat into different object
dat_full <- dat

# change missing WG values to "UCL" so they aren't blank
dat_full$WG[dat_full$WG == ""] <- "UCL"

# remove WAT POLYTYPE
dat_full <- dat_full[dat_full$POLYTYPE != "WAT",]

# change all non-numeric variables to factor
dat_full[sapply(dat_full, is.character)] <- lapply(dat_full[sapply(dat_full, is.character)], 
                                         as.factor)

# create ID col based on row names
dat_full$ID <- rownames(dat_full)

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen_2_30_70.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

################################################
###RUN IMPUTATION ALGORITHM SPL METRICS FIRST###
################################################

# change missing WG values to "UCL" so they aren't blank
dat$WG[dat$WG == ""] <- "UCL"

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

# build dataframe of ancillary data we actually want for the comparison
anci <- dat[, c('HT', 'CC', 'BA', 'AGE', 'POLYID', 'POLYTYPE', 'WG')]

# run imputation algorithm
rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y))

# plot output
plot(rf, vars = yvars(rf))

########################################
###COMPARISON WITH ENTIRE FRI DATASET###
########################################

# this will currently only compare with FRI polygons that were not used in the model
# I think that makes sense?

# run imputation over new polygons
# we want to keep observed cols to see relationship
imp_full <- predict(object = rf,
                     newdata = dat_full,
                     ancillaryData = anci,
                     observed = F)

# create ID col based on row names
imp_full$ID <- rownames(imp_full)

# add observed cols
imp_full <- imp_full %>% add_column(HT_o = dat_full$HT[dat_full$ID %in% imp_full$ID],
                                    CC_o = dat_full$CC[dat_full$ID %in% imp_full$ID],
                                    BA_o = dat_full$BA[dat_full$ID %in% imp_full$ID],
                                    AGE_o = dat_full$AGE[dat_full$ID %in% imp_full$ID],
                                    POLYID_o = dat_full$POLYID[dat_full$ID %in% imp_full$ID],
                                    POLYTYPE_o = dat_full$POLYTYPE[dat_full$ID %in% imp_full$ID],
                                    WG_o = dat_full$WG[dat_full$ID %in% imp_full$ID],
                                    ID_o = dat_full$ID[dat_full$ID %in% imp_full$ID])

# compare input data and imputed data plots

# load get_density function
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

imp_full$density <- get_density(imp_full$HT_o, imp_full$HT, n = 50)

png('D:/ontario_inventory/imputation/test_imp_spl_efi/ht_spl_full.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(HT_o, HT, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 35) + ylim(0, 35)
dev.off()

imp_full$density <- get_density(imp_full$CC_o, imp_full$CC, n = 50)

png('D:/ontario_inventory/imputation/test_imp_spl_efi/cc_spl_full.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(CC_o, CC, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 100) + ylim(0, 100)
dev.off()

imp_full$density <- get_density(imp_full$BA_o, imp_full$BA, n = 50)

png('D:/ontario_inventory/imputation/test_imp_spl_efi/ba_spl_full.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(BA_o, BA, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 60) + ylim(0, 60)
dev.off()

imp_full$density <- get_density(imp_full$AGE_o, imp_full$AGE, n = 50)

png('D:/ontario_inventory/imputation/test_imp_spl_efi/age_spl_full.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(AGE_o, AGE, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 175) + ylim(0, 175)
dev.off()

# calculate rmse different variables
rmse <- function(obs, est) sqrt(mean((obs - est)^2))

rmse_df <- data.frame(AGE = c(rmse(imp_full$AGE_o, imp_full$AGE),
                              rmse(imp_full$AGE_o, imp_full$AGE)/sd(imp_full$AGE_o)),
                      HT = c(rmse(imp_full$HT_o, imp_full$HT),
                             rmse(imp_full$HT_o, imp_full$HT)/sd(imp_full$HT_o)),
                      BA = c(rmse(imp_full$BA_o, imp_full$BA),
                             rmse(imp_full$BA_o, imp_full$BA)/sd(imp_full$BA_o)),
                      CC = c(rmse(imp_full$CC_o, imp_full$CC),
                             rmse(imp_full$CC_o, imp_full$CC)/sd(imp_full$CC_o)))

rownames(rmse_df) <- c('raw_rmsd', 'scaled_rmsd')

# write table
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/rmsd_spl_full.csv')

# create df of POLYTYPE
polyt <- data.frame(obs = imp_full$POLYTYPE_o,
                    est = imp_full$POLYTYPE)

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

# calculate kappa coefficient
kappa <- function(m) {
  N <- sum(m)
  No <- sum(diag(m))
  Ne <- 1 / N * sum(colSums(m) * rowSums(m))
  return( (No - Ne) / (N - Ne) )
}

kap <- kappa(accmat_pt) %>% round(2)

# add kappa to matrix
accmat_pt_ext <- rbind(accmat_pt_ext, c(kap, rep('', NCOL(accmat_pt_ext) - 1)))
rownames(accmat_pt_ext)[NROW(accmat_pt_ext)] <- 'Kappa'

# write to csv
write.csv(accmat_pt_ext, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/accmat_pt_spl_full.csv')

# create df of WG
wg <- data.frame(obs = imp_full$WG_o,
                 est = imp_full$WG)

# make estimate levels match obs levels
levels(wg$est) <- c(levels(wg$est), levels(wg$obs)[!(levels(wg$obs) %in% levels(wg$est))])
wg$est <- factor(wg$est, levels = levels(wg$obs))

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

# calculate kappa coefficient
kap <- kappa(accmat_wg) %>% round(2)

# add kappa to matrix
accmat_wg_ext <- rbind(accmat_wg_ext, c(kap, rep('', NCOL(accmat_wg_ext) - 1)))
rownames(accmat_wg_ext)[NROW(accmat_wg_ext)] <- 'Kappa'

# write to csv
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/accmat_wg_spl_full.csv')

######################################
###COMPARISON WITH SCREENED DATASET###
######################################

# lets also apply the imputation to dataset of FRI polygons used to build the model

# run imputation over new polygons
imp_fri <- impute(object = rf,
                  ancillaryData = anci,
                  observed = F)

# save dat and imp_fri
# write.csv(dat, file = 'D:/ontario_inventory/imputation/test_imp_wg/dat_out.csv', row.names = F)
# write.csv(imp_fri, file = 'D:/ontario_inventory/imputation/test_imp_wg/dat_out_imp.csv', row.names = F)

# compare input data and imputed data plots
png('D:/ontario_inventory/imputation/test_imp_spl_efi/ht_spl.png', width = 480)
plot(dat$HT, imp_fri$HT)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_spl_efi/cc_spl.png', width = 480)
plot(dat$CC, imp_fri$CC)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_spl_efi/ba_spl.png', width = 480)
plot(dat$BA, imp_fri$BA)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_spl_efi/age_spl.png', width = 480)
plot(dat$AGE, imp_fri$AGE)
dev.off()

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
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/rmsd_spl.csv')

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

# calculate kappa coefficient
kappa <- function(m) {
  N <- sum(m)
  No <- sum(diag(m))
  Ne <- 1 / N * sum(colSums(m) * rowSums(m))
  return( (No - Ne) / (N - Ne) )
}

kap <- kappa(accmat_pt) %>% round(2)

# add kappa to matrix
accmat_pt_ext <- rbind(accmat_pt_ext, c(kap, rep('', NCOL(accmat_pt_ext) - 1)))
rownames(accmat_pt_ext)[NROW(accmat_pt_ext)] <- 'Kappa'

# write to csv
write.csv(accmat_pt_ext, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/accmat_pt_spl.csv')

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

# calculate kappa coefficient
kap <- kappa(accmat_wg) %>% round(2)

# add kappa to matrix
accmat_wg_ext <- rbind(accmat_wg_ext, c(kap, rep('', NCOL(accmat_wg_ext) - 1)))
rownames(accmat_wg_ext)[NROW(accmat_wg_ext)] <- 'Kappa'

# write to csv
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/accmat_wg_spl.csv')

#################################################
###RUN IMPUTATION ALGORITHM EFI METRICS SECOND###
#################################################

# build dataframe of reference variables
ref_vars <- c('agb', 'ba', 'dens', 'lor', 'qmdbh', 'top_height', 'v', 'v_merch')
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

########################################
###COMPARISON WITH ENTIRE FRI DATASET###
########################################

# this will currently only compare with FRI polygons that were not used in the model
# I think that makes sense?

# run imputation over new polygons
# we want to keep observed cols to see relationship
imp_full <- predict(object = rf,
                    newdata = dat_full,
                    ancillaryData = anci,
                    observed = F)

# create ID col based on row names
imp_full$ID <- rownames(imp_full)

# add observed cols
imp_full <- imp_full %>% add_column(HT_o = dat_full$HT[dat_full$ID %in% imp_full$ID],
                                    CC_o = dat_full$CC[dat_full$ID %in% imp_full$ID],
                                    BA_o = dat_full$BA[dat_full$ID %in% imp_full$ID],
                                    AGE_o = dat_full$AGE[dat_full$ID %in% imp_full$ID],
                                    POLYID_o = dat_full$POLYID[dat_full$ID %in% imp_full$ID],
                                    POLYTYPE_o = dat_full$POLYTYPE[dat_full$ID %in% imp_full$ID],
                                    WG_o = dat_full$WG[dat_full$ID %in% imp_full$ID],
                                    ID_o = dat_full$ID[dat_full$ID %in% imp_full$ID])

# compare input data and imputed data plots

# load get_density function
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

imp_full$density <- get_density(imp_full$HT_o, imp_full$HT, n = 50)

png('D:/ontario_inventory/imputation/test_imp_spl_efi/ht_efi_full.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(HT_o, HT, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 35) + ylim(0, 35)
dev.off()

imp_full$density <- get_density(imp_full$CC_o, imp_full$CC, n = 50)

png('D:/ontario_inventory/imputation/test_imp_spl_efi/cc_efi_full.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(CC_o, CC, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 100) + ylim(0, 100)
dev.off()

imp_full$density <- get_density(imp_full$BA_o, imp_full$BA, n = 50)

png('D:/ontario_inventory/imputation/test_imp_spl_efi/ba_efi_full.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(BA_o, BA, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 60) + ylim(0, 60)
dev.off()

imp_full$density <- get_density(imp_full$AGE_o, imp_full$AGE, n = 50)

png('D:/ontario_inventory/imputation/test_imp_spl_efi/age_efi_full.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(AGE_o, AGE, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 175) + ylim(0, 175)
dev.off()

# calculate rmse different variables
rmse <- function(obs, est) sqrt(mean((obs - est)^2))

rmse_df <- data.frame(AGE = c(rmse(imp_full$AGE_o, imp_full$AGE),
                              rmse(imp_full$AGE_o, imp_full$AGE)/sd(imp_full$AGE_o)),
                      HT = c(rmse(imp_full$HT_o, imp_full$HT),
                             rmse(imp_full$HT_o, imp_full$HT)/sd(imp_full$HT_o)),
                      BA = c(rmse(imp_full$BA_o, imp_full$BA),
                             rmse(imp_full$BA_o, imp_full$BA)/sd(imp_full$BA_o)),
                      CC = c(rmse(imp_full$CC_o, imp_full$CC),
                             rmse(imp_full$CC_o, imp_full$CC)/sd(imp_full$CC_o)))

rownames(rmse_df) <- c('raw_rmsd', 'scaled_rmsd')

# write table
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/rmsd_efi_full.csv')

# create df of POLYTYPE
polyt <- data.frame(obs = imp_full$POLYTYPE_o,
                    est = imp_full$POLYTYPE)

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

# calculate kappa coefficient
kappa <- function(m) {
  N <- sum(m)
  No <- sum(diag(m))
  Ne <- 1 / N * sum(colSums(m) * rowSums(m))
  return( (No - Ne) / (N - Ne) )
}

kap <- kappa(accmat_pt) %>% round(2)

# add kappa to matrix
accmat_pt_ext <- rbind(accmat_pt_ext, c(kap, rep('', NCOL(accmat_pt_ext) - 1)))
rownames(accmat_pt_ext)[NROW(accmat_pt_ext)] <- 'Kappa'

# write to csv
write.csv(accmat_pt_ext, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/accmat_pt_efi_full.csv')

# create df of WG
wg <- data.frame(obs = imp_full$WG_o,
                 est = imp_full$WG)

# make estimate levels match obs levels
levels(wg$est) <- c(levels(wg$est), levels(wg$obs)[!(levels(wg$obs) %in% levels(wg$est))])
wg$est <- factor(wg$est, levels = levels(wg$obs))

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

# calculate kappa coefficient
kap <- kappa(accmat_wg) %>% round(2)

# add kappa to matrix
accmat_wg_ext <- rbind(accmat_wg_ext, c(kap, rep('', NCOL(accmat_wg_ext) - 1)))
rownames(accmat_wg_ext)[NROW(accmat_wg_ext)] <- 'Kappa'

# write to csv
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/accmat_wg_efi_full.csv')

######################################
###COMPARISON WITH SCREENED DATASET###
######################################

# lets also apply the imputation to dataset of FRI polygons used to build the model

# run imputation over new polygons
imp_fri <- impute(object = rf,
                  ancillaryData = anci,
                  observed = F)

# save dat and imp_fri
# write.csv(dat, file = 'D:/ontario_inventory/imputation/test_imp_wg/dat_out.csv', row.names = F)
# write.csv(imp_fri, file = 'D:/ontario_inventory/imputation/test_imp_wg/dat_out_imp.csv', row.names = F)

# compare input data and imputed data plots
png('D:/ontario_inventory/imputation/test_imp_spl_efi/ht_efi.png', width = 480)
plot(dat$HT, imp_fri$HT)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_spl_efi/cc_efi.png', width = 480)
plot(dat$CC, imp_fri$CC)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_spl_efi/ba_efi.png', width = 480)
plot(dat$BA, imp_fri$BA)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_spl_efi/age_efi.png', width = 480)
plot(dat$AGE, imp_fri$AGE)
dev.off()

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
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/rmsd_efi.csv')

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

# calculate kappa coefficient
kappa <- function(m) {
  N <- sum(m)
  No <- sum(diag(m))
  Ne <- 1 / N * sum(colSums(m) * rowSums(m))
  return( (No - Ne) / (N - Ne) )
}

kap <- kappa(accmat_pt) %>% round(2)

# add kappa to matrix
accmat_pt_ext <- rbind(accmat_pt_ext, c(kap, rep('', NCOL(accmat_pt_ext) - 1)))
rownames(accmat_pt_ext)[NROW(accmat_pt_ext)] <- 'Kappa'

# write to csv
write.csv(accmat_pt_ext, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/accmat_pt_efi.csv')

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

# calculate kappa coefficient
kap <- kappa(accmat_wg) %>% round(2)

# add kappa to matrix
accmat_wg_ext <- rbind(accmat_wg_ext, c(kap, rep('', NCOL(accmat_wg_ext) - 1)))
rownames(accmat_wg_ext)[NROW(accmat_wg_ext)] <- 'Kappa'

# write to csv
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/accmat_wg_efi.csv')

################################################
###RUN IMPUTATION ALGORITHM ALL METRICS THIRD###
################################################

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

########################################
###COMPARISON WITH ENTIRE FRI DATASET###
########################################

# this will currently only compare with FRI polygons that were not used in the model
# I think that makes sense?

# run imputation over new polygons
# we want to keep observed cols to see relationship
imp_full <- predict(object = rf,
                    newdata = dat_full,
                    ancillaryData = anci,
                    observed = F)

# create ID col based on row names
imp_full$ID <- rownames(imp_full)

# add observed cols
imp_full <- imp_full %>% add_column(HT_o = dat_full$HT[dat_full$ID %in% imp_full$ID],
                                    CC_o = dat_full$CC[dat_full$ID %in% imp_full$ID],
                                    BA_o = dat_full$BA[dat_full$ID %in% imp_full$ID],
                                    AGE_o = dat_full$AGE[dat_full$ID %in% imp_full$ID],
                                    POLYID_o = dat_full$POLYID[dat_full$ID %in% imp_full$ID],
                                    POLYTYPE_o = dat_full$POLYTYPE[dat_full$ID %in% imp_full$ID],
                                    WG_o = dat_full$WG[dat_full$ID %in% imp_full$ID],
                                    ID_o = dat_full$ID[dat_full$ID %in% imp_full$ID])

# compare input data and imputed data plots

# load get_density function
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

imp_full$density <- get_density(imp_full$HT_o, imp_full$HT, n = 50)

png('D:/ontario_inventory/imputation/test_imp_spl_efi/ht_all_full.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(HT_o, HT, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 35) + ylim(0, 35)
dev.off()

imp_full$density <- get_density(imp_full$CC_o, imp_full$CC, n = 50)

png('D:/ontario_inventory/imputation/test_imp_spl_efi/cc_all_full.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(CC_o, CC, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 100) + ylim(0, 100)
dev.off()

imp_full$density <- get_density(imp_full$BA_o, imp_full$BA, n = 50)

png('D:/ontario_inventory/imputation/test_imp_spl_efi/ba_all_full.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(BA_o, BA, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 60) + ylim(0, 60)
dev.off()

imp_full$density <- get_density(imp_full$AGE_o, imp_full$AGE, n = 50)

png('D:/ontario_inventory/imputation/test_imp_spl_efi/age_all_full.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(AGE_o, AGE, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 175) + ylim(0, 175)
dev.off()

# calculate rmse different variables
rmse <- function(obs, est) sqrt(mean((obs - est)^2))

rmse_df <- data.frame(AGE = c(rmse(imp_full$AGE_o, imp_full$AGE),
                              rmse(imp_full$AGE_o, imp_full$AGE)/sd(imp_full$AGE_o)),
                      HT = c(rmse(imp_full$HT_o, imp_full$HT),
                             rmse(imp_full$HT_o, imp_full$HT)/sd(imp_full$HT_o)),
                      BA = c(rmse(imp_full$BA_o, imp_full$BA),
                             rmse(imp_full$BA_o, imp_full$BA)/sd(imp_full$BA_o)),
                      CC = c(rmse(imp_full$CC_o, imp_full$CC),
                             rmse(imp_full$CC_o, imp_full$CC)/sd(imp_full$CC_o)))

rownames(rmse_df) <- c('raw_rmsd', 'scaled_rmsd')

# write table
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/rmsd_all_full.csv')

# create df of POLYTYPE
polyt <- data.frame(obs = imp_full$POLYTYPE_o,
                    est = imp_full$POLYTYPE)

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

# calculate kappa coefficient
kappa <- function(m) {
  N <- sum(m)
  No <- sum(diag(m))
  Ne <- 1 / N * sum(colSums(m) * rowSums(m))
  return( (No - Ne) / (N - Ne) )
}

kap <- kappa(accmat_pt) %>% round(2)

# add kappa to matrix
accmat_pt_ext <- rbind(accmat_pt_ext, c(kap, rep('', NCOL(accmat_pt_ext) - 1)))
rownames(accmat_pt_ext)[NROW(accmat_pt_ext)] <- 'Kappa'

# write to csv
write.csv(accmat_pt_ext, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/accmat_pt_all_full.csv')

# create df of WG
wg <- data.frame(obs = imp_full$WG_o,
                 est = imp_full$WG)

# make estimate levels match obs levels
levels(wg$est) <- c(levels(wg$est), levels(wg$obs)[!(levels(wg$obs) %in% levels(wg$est))])
wg$est <- factor(wg$est, levels = levels(wg$obs))

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

# calculate kappa coefficient
kap <- kappa(accmat_wg) %>% round(2)

# add kappa to matrix
accmat_wg_ext <- rbind(accmat_wg_ext, c(kap, rep('', NCOL(accmat_wg_ext) - 1)))
rownames(accmat_wg_ext)[NROW(accmat_wg_ext)] <- 'Kappa'

# write to csv
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/accmat_wg_all_full.csv')

######################################
###COMPARISON WITH SCREENED DATASET###
######################################

# lets also apply the imputation to dataset of FRI polygons used to build the model

# run imputation over new polygons
imp_fri <- impute(object = rf,
                  ancillaryData = anci,
                  observed = F)

# save dat and imp_fri
# write.csv(dat, file = 'D:/ontario_inventory/imputation/test_imp_wg/dat_out.csv', row.names = F)
# write.csv(imp_fri, file = 'D:/ontario_inventory/imputation/test_imp_wg/dat_out_imp.csv', row.names = F)

# compare input data and imputed data plots
png('D:/ontario_inventory/imputation/test_imp_spl_efi/ht_all.png', width = 480)
plot(dat$HT, imp_fri$HT)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_spl_efi/cc_all.png', width = 480)
plot(dat$CC, imp_fri$CC)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_spl_efi/ba_all.png', width = 480)
plot(dat$BA, imp_fri$BA)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_spl_efi/age_all.png', width = 480)
plot(dat$AGE, imp_fri$AGE)
dev.off()

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
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/rmsd_all.csv')

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

# calculate kappa coefficient
kappa <- function(m) {
  N <- sum(m)
  No <- sum(diag(m))
  Ne <- 1 / N * sum(colSums(m) * rowSums(m))
  return( (No - Ne) / (N - Ne) )
}

kap <- kappa(accmat_pt) %>% round(2)

# add kappa to matrix
accmat_pt_ext <- rbind(accmat_pt_ext, c(kap, rep('', NCOL(accmat_pt_ext) - 1)))
rownames(accmat_pt_ext)[NROW(accmat_pt_ext)] <- 'Kappa'

# write to csv
write.csv(accmat_pt_ext, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/accmat_pt_all.csv')

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

# calculate kappa coefficient
kap <- kappa(accmat_wg) %>% round(2)

# add kappa to matrix
accmat_wg_ext <- rbind(accmat_wg_ext, c(kap, rep('', NCOL(accmat_wg_ext) - 1)))
rownames(accmat_wg_ext)[NROW(accmat_wg_ext)] <- 'Kappa'

# write to csv
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/test_imp_spl_efi/accmat_wg_all.csv')