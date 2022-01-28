# This code runs an imputation algorithm on LiDAR attributes connecting them to interpreter
# derived response attributes resulting in a model that can be applied to forest stands
# segmented from LiDAR rasters

# The purpose is to compare the different data screening procedures, which have been
# updated to only use FOR polytype and only SPL variables for screening (p95 instead of lor)

# load packages
library(terra)
library(tidyverse)
library(yaImpute)
library(randomForest)

###################################
###SCREENED FOR 10-90% RESIDUALS###
###################################

###################################################
###LOAD POLYGON DATASET AND ADD LIDAR ATTRIBUTES###
###################################################

# load extracted data frame
load('D:/ontario_inventory/dat/dat_fri_100.RData')

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

# save dat into different object
dat_full <- dat

# change missing WG values to "UCL" so they aren't blank
dat_full$WG[dat_full$WG == ""] <- "UCL"

# create further separation of data to impute
# only FOR POLYTYPE
dat_imp <- dat_full[dat_full$POLYTYPE == "FOR",]

# change all non-numeric variables to factor
dat_imp[sapply(dat_imp, is.character)] <- lapply(dat_imp[sapply(dat_imp, is.character)], 
                                                 as.factor)

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen_1_2a_10perc_for.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# SPL vars only

# change missing WG values to "UCL" so they aren't blank
# dat$WG[dat$WG == ""] <- "UCL"

# change all non-numeric variables to factor
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], 
                                         as.factor)

# build dataframe of reference variables
ref_vars <- c('cc', 'p95', 'avg', 'qav', 'cv', 'kur', 'max', 'ske')
x <- dat[, ref_vars]

# build dataframe of target variables
tar_vars <- c('HT', 'CC', 'BA', 'WG')
y <- dat[, tar_vars]

# build dataframe of ancillary data
anci <- dat[, !(colnames(dat) %in% ref_vars)]

# run imputation algorithm
rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y))

# plot output
plot(rf, vars = yvars(rf))

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

png('D:/ontario_inventory/imputation/test_imp_for_only/ht_spl_full_10_90.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(HT_o, HT, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 35) + ylim(0, 35)
dev.off()

imp_full$density <- get_density(imp_full$CC_o, imp_full$CC, n = 50)

png('D:/ontario_inventory/imputation/test_imp_for_only/cc_spl_full_10_90.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(CC_o, CC, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 100) + ylim(0, 100)
dev.off()

imp_full$density <- get_density(imp_full$BA_o, imp_full$BA, n = 50)

png('D:/ontario_inventory/imputation/test_imp_for_only/ba_spl_full_10_90.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(BA_o, BA, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 60) + ylim(0, 60)
dev.off()

imp_full$density <- get_density(imp_full$AGE_o, imp_full$AGE, n = 50)

png('D:/ontario_inventory/imputation/test_imp_for_only/age_spl_full_10_90.png', width = 480)
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
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_for_only/rmsd_spl_full_10_90.csv')

###################################
###SCREENED FOR 30-70% RESIDUALS###
###################################

###################################################
###LOAD POLYGON DATASET AND ADD LIDAR ATTRIBUTES###
###################################################

# load extracted data frame
load('D:/ontario_inventory/dat/dat_fri_100.RData')

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

# save dat into different object
dat_full <- dat

# change missing WG values to "UCL" so they aren't blank
dat_full$WG[dat_full$WG == ""] <- "UCL"

# create further separation of data to impute
# only FOR POLYTYPE
dat_imp <- dat_full[dat_full$POLYTYPE == "FOR",]

# change all non-numeric variables to factor
dat_imp[sapply(dat_imp, is.character)] <- lapply(dat_imp[sapply(dat_imp, is.character)], 
                                                 as.factor)

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen_1_2a_30perc_for.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# SPL vars only

# change missing WG values to "UCL" so they aren't blank
# dat$WG[dat$WG == ""] <- "UCL"

# change all non-numeric variables to factor
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], 
                                         as.factor)

# build dataframe of reference variables
ref_vars <- c('cc', 'p95', 'avg', 'qav', 'cv', 'kur', 'max', 'ske')
x <- dat[, ref_vars]

# build dataframe of target variables
tar_vars <- c('HT', 'CC', 'BA', 'WG')
y <- dat[, tar_vars]

# build dataframe of ancillary data
anci <- dat[, !(colnames(dat) %in% ref_vars)]

# run imputation algorithm
rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y))

# plot output
plot(rf, vars = yvars(rf))

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

png('D:/ontario_inventory/imputation/test_imp_for_only/ht_spl_full_30_70.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(HT_o, HT, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 35) + ylim(0, 35)
dev.off()

imp_full$density <- get_density(imp_full$CC_o, imp_full$CC, n = 50)

png('D:/ontario_inventory/imputation/test_imp_for_only/cc_spl_full_30_70.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(CC_o, CC, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 100) + ylim(0, 100)
dev.off()

imp_full$density <- get_density(imp_full$BA_o, imp_full$BA, n = 50)

png('D:/ontario_inventory/imputation/test_imp_for_only/ba_spl_full_30_70.png', width = 480)
ggplot(imp_full) + 
  geom_point(aes(BA_o, BA, color = density)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 60) + ylim(0, 60)
dev.off()

imp_full$density <- get_density(imp_full$AGE_o, imp_full$AGE, n = 50)

png('D:/ontario_inventory/imputation/test_imp_for_only/age_spl_full_30_70.png', width = 480)
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
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_for_only/rmsd_spl_full_30_70.csv')

