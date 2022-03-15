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
library(viridis)

# load function to calculate kappa coefficient
kappa <- function(m) {
  N <- sum(m)
  No <- sum(diag(m))
  Ne <- 1 / N * sum(colSums(m) * rowSums(m))
  return( (No - Ne) / (N - Ne) )
}

###################################
###SCREENED FOR 10-90% RESIDUALS###
###################################

###################################################
###LOAD POLYGON DATASET AND ADD LIDAR ATTRIBUTES###
###################################################

# load extracted data frame
# dat_fri_100 is loaded with LiDAR arributes derived from 
# masking pixels that do not fall 100 % inside polygons
load('D:/ontario_inventory/dat/dat_fri_100.RData')

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen_1bc_2a_2pc_10perc_wat_ucl.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# SPL vars only

# set missing WG values to UCL
dat$WG[dat$WG == ""] <- "UCL"

# change all non-numeric variables to factor
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], 
                                         as.factor)

# build dataframe of reference variables
ref_vars <- c('cc', 'p95', 'avg', 'qav', 'cv', 'kur', 'max', 'ske')
x <- dat[, ref_vars]

# build dataframe of target variables
tar_vars <- c('HT', 'WG')
y <- dat[, tar_vars]

# build dataframe of ancillary data
anci <- dat[, !(colnames(dat) %in% ref_vars)]

# run imputation algorithm
rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y))

# plot output
plot(rf, vars = yvars(rf))

rmsd(rf)

#################################
###COMPARISON WITH FRI DATASET###
#################################

# apply the imputation algorithm to the training/forest data
# we are not comparing to the other data because it doesn't 
# meet our criteria for clean forested polygons

# run imputation over new polygons
# we want to keep observed cols to see relationship
imp <- impute(object = rf,
              ancillaryData = anci,
              observed = T)

# compare input data and imputed data plots

# load get_density function
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# calculate density values
imp$dens_ht <- get_density(imp$HT.o, imp$HT, n = 50)
imp$dens_cc <- get_density(imp$CC.o, imp$CC, n = 50)
imp$dens_ba <- get_density(imp$BA.o, imp$BA, n = 50)
imp$dens_age <- get_density(imp$AGE.o, imp$AGE, n = 50)

# plot
png('D:/ontario_inventory/imputation/test_imp_for_only/ht_10_90_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(HT.o, HT, color = dens_ht)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 35) + ylim(0, 35)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/cc_10_90_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(CC.o, CC, color = dens_cc)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 100) + ylim(0, 100)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/ba_10_90_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(BA.o, BA, color = dens_ba)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 60) + ylim(0, 60)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/age_10_90_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(AGE.o, AGE, color = dens_age)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 200) + ylim(0, 200)
dev.off()

# calculate rmse different variables
rmse <- function(obs, est) sqrt(mean((obs - est)^2))

rmse_df <- data.frame(AGE = c(rmse(imp$AGE.o, imp$AGE),
                              rmse(imp$AGE.o, imp$AGE)/sd(imp$AGE.o)),
                      HT = c(rmse(imp$HT.o, imp$HT),
                             rmse(imp$HT.o, imp$HT)/sd(imp$HT.o)),
                      BA = c(rmse(imp$BA.o, imp$BA),
                             rmse(imp$BA.o, imp$BA)/sd(imp$BA.o)),
                      CC = c(rmse(imp$CC.o, imp$CC),
                             rmse(imp$CC.o, imp$CC)/sd(imp$CC.o)))

rownames(rmse_df) <- c('raw_rmsd', 'scaled_rmsd')

# write table
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_for_only/rmsd_10_90.csv')

# create df of WG
wg <- data.frame(obs = imp$WG.o,
                 est = imp$WG)

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
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/test_imp_for_only/accmat_wg_10_90.csv')

# create df of DEVSTAGE
dev <- data.frame(obs = imp$DEVSTAGE.o,
                 est = imp$DEVSTAGE)

# make estimate levels match obs levels
levels(dev$est) <- c(levels(dev$est), levels(dev$obs)[!(levels(dev$obs) %in% levels(dev$est))])
dev$est <- factor(dev$est, levels = levels(dev$obs))

# create column of match or not
dev$match <- dev$obs == dev$est

# total percent of matching DEVSTAGE
NROW(dev[dev$match == T,])/NROW(dev)

# check count of different working groups
plyr::count(dev, 'obs')
plyr::count(dev, 'est')

# build accuracy table
accmat_dev <- table("pred" = dev$est, "ref" = dev$obs)

# UA
ua_dev <- diag(accmat_dev) / rowSums(accmat_dev) * 100

# PA
pa_dev <- diag(accmat_dev) / colSums(accmat_dev) * 100

# OA
oa_dev <- sum(diag(accmat_dev)) / sum(accmat_dev) * 100

# build confusion matrix
accmat_dev_ext <- addmargins(accmat_dev)
accmat_dev_ext <- rbind(accmat_dev_ext, "Users" = c(pa_dev, NA))
accmat_dev_ext <- cbind(accmat_dev_ext, "Producers" = c(ua_dev, NA, oa_dev))
#colnames(accmat_wg_ext) <- c(levels(as.factor(shp.train$classes)), "Sum", "PA")
#rownames(accmat_wg_ext) <- c(levels(as.factor(shp.train$classes)), "Sum", "UA")
accmat_dev_ext <- round(accmat_dev_ext, digits = 1)
dimnames(accmat_dev_ext) <- list("Prediction" = colnames(accmat_dev_ext),
                                "Reference" = rownames(accmat_dev_ext))
class(accmat_dev_ext) <- "table"
accmat_dev_ext

# calculate kappa coefficient
kap <- kappa(accmat_dev) %>% round(2)

# add kappa to matrix
accmat_dev_ext <- rbind(accmat_dev_ext, c(kap, rep('', NCOL(accmat_dev_ext) - 1)))
rownames(accmat_dev_ext)[NROW(accmat_dev_ext)] <- 'Kappa'

# write to csv
write.csv(accmat_dev_ext, file = 'D:/ontario_inventory/imputation/test_imp_for_only/accmat_dev_10_90_wg_devstage.csv')

###################################
###SCREENED FOR 30-70% RESIDUALS###
###################################

###################################################
###LOAD POLYGON DATASET AND ADD LIDAR ATTRIBUTES###
###################################################

# load extracted data frame
# dat_fri_100 is loaded with LiDAR arributes derived from 
# masking pixels that do not fall 100 % inside polygons
load('D:/ontario_inventory/dat/dat_fri_100.RData')

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen_1bc_2apc_30perc_for.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# SPL vars only

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

# apply the imputation algorithm to the training/forest data
# we are not comparing to the other data because it doesn't 
# meet our criteria for clean forested polygons

# run imputation over new polygons
# we want to keep observed cols to see relationship
imp <- impute(object = rf,
              ancillaryData = anci,
              observed = T)

# compare input data and imputed data plots

# load get_density function
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# calculate density values
imp$dens_ht <- get_density(imp$HT.o, imp$HT, n = 50)
imp$dens_cc <- get_density(imp$CC.o, imp$CC, n = 50)
imp$dens_ba <- get_density(imp$BA.o, imp$BA, n = 50)
imp$dens_age <- get_density(imp$AGE.o, imp$AGE, n = 50)

# plot
png('D:/ontario_inventory/imputation/test_imp_for_only/ht_30_70_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(HT.o, HT, color = dens_ht)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 35) + ylim(0, 35)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/cc_30_70_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(CC.o, CC, color = dens_cc)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 100) + ylim(0, 100)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/ba_30_70_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(BA.o, BA, color = dens_ba)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 60) + ylim(0, 60)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/age_30_70_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(AGE.o, AGE, color = dens_age)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 200) + ylim(0, 200)
dev.off()

# calculate rmse different variables
rmse <- function(obs, est) sqrt(mean((obs - est)^2))

rmse_df <- data.frame(AGE = c(rmse(imp$AGE.o, imp$AGE),
                              rmse(imp$AGE.o, imp$AGE)/sd(imp$AGE.o)),
                      HT = c(rmse(imp$HT.o, imp$HT),
                             rmse(imp$HT.o, imp$HT)/sd(imp$HT.o)),
                      BA = c(rmse(imp$BA.o, imp$BA),
                             rmse(imp$BA.o, imp$BA)/sd(imp$BA.o)),
                      CC = c(rmse(imp$CC.o, imp$CC),
                             rmse(imp$CC.o, imp$CC)/sd(imp$CC.o)))

rownames(rmse_df) <- c('raw_rmsd', 'scaled_rmsd')

# write table
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_for_only/rmsd_30_70.csv')

# create df of WG
wg <- data.frame(obs = imp$WG.o,
                 est = imp$WG)

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
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/test_imp_for_only/accmat_wg_30_70.csv')

#########################################
###SCREENED FOR 10-90% RESIDUALS NO BA###
#########################################

###################################################
###LOAD POLYGON DATASET AND ADD LIDAR ATTRIBUTES###
###################################################

# load extracted data frame
# dat_fri_100 is loaded with LiDAR arributes derived from 
# masking pixels that do not fall 100 % inside polygons
load('D:/ontario_inventory/dat/dat_fri_100.RData')

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen_1bc_2apc_10perc_for.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# SPL vars only

# change all non-numeric variables to factor
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], 
                                         as.factor)

# build dataframe of reference variables
ref_vars <- c('cc', 'p95', 'avg', 'qav', 'cv', 'kur', 'max', 'ske')
x <- dat[, ref_vars]

# build dataframe of target variables
tar_vars <- c('HT', 'CC', 'WG')
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

# apply the imputation algorithm to the training/forest data
# we are not comparing to the other data because it doesn't 
# meet our criteria for clean forested polygons

# run imputation over new polygons
# we want to keep observed cols to see relationship
imp <- impute(object = rf,
              ancillaryData = anci,
              observed = T)

# compare input data and imputed data plots

# load get_density function
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# calculate density values
imp$dens_ht <- get_density(imp$HT.o, imp$HT, n = 50)
imp$dens_cc <- get_density(imp$CC.o, imp$CC, n = 50)
imp$dens_ba <- get_density(imp$BA.o, imp$BA, n = 50)
imp$dens_age <- get_density(imp$AGE.o, imp$AGE, n = 50)

# plot
png('D:/ontario_inventory/imputation/test_imp_for_only/ht_10_90_no_ba_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(HT.o, HT, color = dens_ht)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 35) + ylim(0, 35)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/cc_10_90_no_ba_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(CC.o, CC, color = dens_cc)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 100) + ylim(0, 100)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/ba_10_90_no_ba_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(BA.o, BA, color = dens_ba)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 60) + ylim(0, 60)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/age_10_90_no_ba_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(AGE.o, AGE, color = dens_age)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 200) + ylim(0, 200)
dev.off()

# calculate rmse different variables
rmse <- function(obs, est) sqrt(mean((obs - est)^2))

rmse_df <- data.frame(AGE = c(rmse(imp$AGE.o, imp$AGE),
                              rmse(imp$AGE.o, imp$AGE)/sd(imp$AGE.o)),
                      HT = c(rmse(imp$HT.o, imp$HT),
                             rmse(imp$HT.o, imp$HT)/sd(imp$HT.o)),
                      BA = c(rmse(imp$BA.o, imp$BA),
                             rmse(imp$BA.o, imp$BA)/sd(imp$BA.o)),
                      CC = c(rmse(imp$CC.o, imp$CC),
                             rmse(imp$CC.o, imp$CC)/sd(imp$CC.o)))

rownames(rmse_df) <- c('raw_rmsd', 'scaled_rmsd')

# write table
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_for_only/rmsd_10_90_no_ba.csv')

# create df of WG
wg <- data.frame(obs = imp$WG.o,
                 est = imp$WG)

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
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/test_imp_for_only/accmat_wg_10_90_no_ba.csv')

#########################################
###SCREENED FOR 30-70% RESIDUALS NO BA###
#########################################

###################################################
###LOAD POLYGON DATASET AND ADD LIDAR ATTRIBUTES###
###################################################

# load extracted data frame
# dat_fri_100 is loaded with LiDAR arributes derived from 
# masking pixels that do not fall 100 % inside polygons
load('D:/ontario_inventory/dat/dat_fri_100.RData')

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen_1bc_2apc_30perc_for.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# SPL vars only

# change all non-numeric variables to factor
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], 
                                         as.factor)

# build dataframe of reference variables
ref_vars <- c('cc', 'p95', 'avg', 'qav', 'cv', 'kur', 'max', 'ske')
x <- dat[, ref_vars]

# build dataframe of target variables
tar_vars <- c('HT', 'CC', 'WG')
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

# apply the imputation algorithm to the training/forest data
# we are not comparing to the other data because it doesn't 
# meet our criteria for clean forested polygons

# run imputation over new polygons
# we want to keep observed cols to see relationship
imp <- impute(object = rf,
              ancillaryData = anci,
              observed = T)

# compare input data and imputed data plots

# load get_density function
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# calculate density values
imp$dens_ht <- get_density(imp$HT.o, imp$HT, n = 50)
imp$dens_cc <- get_density(imp$CC.o, imp$CC, n = 50)
imp$dens_ba <- get_density(imp$BA.o, imp$BA, n = 50)
imp$dens_age <- get_density(imp$AGE.o, imp$AGE, n = 50)

# plot
png('D:/ontario_inventory/imputation/test_imp_for_only/ht_30_70_no_ba_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(HT.o, HT, color = dens_ht)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 35) + ylim(0, 35)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/cc_30_70_no_ba_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(CC.o, CC, color = dens_cc)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 100) + ylim(0, 100)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/ba_30_70_no_ba_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(BA.o, BA, color = dens_ba)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 60) + ylim(0, 60)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/age_30_70_no_ba_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(AGE.o, AGE, color = dens_age)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 200) + ylim(0, 200)
dev.off()

# calculate rmse different variables
rmse <- function(obs, est) sqrt(mean((obs - est)^2))

rmse_df <- data.frame(AGE = c(rmse(imp$AGE.o, imp$AGE),
                              rmse(imp$AGE.o, imp$AGE)/sd(imp$AGE.o)),
                      HT = c(rmse(imp$HT.o, imp$HT),
                             rmse(imp$HT.o, imp$HT)/sd(imp$HT.o)),
                      BA = c(rmse(imp$BA.o, imp$BA),
                             rmse(imp$BA.o, imp$BA)/sd(imp$BA.o)),
                      CC = c(rmse(imp$CC.o, imp$CC),
                             rmse(imp$CC.o, imp$CC)/sd(imp$CC.o)))

rownames(rmse_df) <- c('raw_rmsd', 'scaled_rmsd')

# write table
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_for_only/rmsd_30_70_no_ba.csv')

# create df of WG
wg <- data.frame(obs = imp$WG.o,
                 est = imp$WG)

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
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/test_imp_for_only/accmat_wg_30_70_no_ba.csv')

###################################################
###SCREENED FOR 10-90% RESIDUALS NO BA 200 TREES###
###################################################

###################################################
###LOAD POLYGON DATASET AND ADD LIDAR ATTRIBUTES###
###################################################

# load extracted data frame
# dat_fri_100 is loaded with LiDAR arributes derived from 
# masking pixels that do not fall 100 % inside polygons
load('D:/ontario_inventory/dat/dat_fri_100.RData')

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen_1bc_2apc_10perc_for.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# SPL vars only

# change all non-numeric variables to factor
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], 
                                         as.factor)

# build dataframe of reference variables
ref_vars <- c('cc', 'p95', 'avg', 'qav', 'cv', 'kur', 'max', 'ske')
x <- dat[, ref_vars]

# build dataframe of target variables
tar_vars <- c('HT', 'CC', 'WG')
y <- dat[, tar_vars]

# build dataframe of ancillary data
anci <- dat[, !(colnames(dat) %in% ref_vars)]

# run imputation algorithm
rf <- yai(x = x, y = y, method = 'randomForest', ntree = 200*NCOL(y))

# plot output
plot(rf, vars = yvars(rf))

#################################
###COMPARISON WITH FRI DATASET###
#################################

# apply the imputation algorithm to the training/forest data
# we are not comparing to the other data because it doesn't 
# meet our criteria for clean forested polygons

# run imputation over new polygons
# we want to keep observed cols to see relationship
imp <- impute(object = rf,
              ancillaryData = anci,
              observed = T)

# compare input data and imputed data plots

# load get_density function
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# calculate density values
imp$dens_ht <- get_density(imp$HT.o, imp$HT, n = 50)
imp$dens_cc <- get_density(imp$CC.o, imp$CC, n = 50)
imp$dens_ba <- get_density(imp$BA.o, imp$BA, n = 50)
imp$dens_age <- get_density(imp$AGE.o, imp$AGE, n = 50)

# plot
png('D:/ontario_inventory/imputation/test_imp_for_only/ht_10_90_no_ba_200_trees_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(HT.o, HT, color = dens_ht)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 35) + ylim(0, 35)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/cc_10_90_no_ba_200_trees_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(CC.o, CC, color = dens_cc)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 100) + ylim(0, 100)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/ba_10_90_no_ba_200_trees_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(BA.o, BA, color = dens_ba)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 60) + ylim(0, 60)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/age_10_90_no_ba_200_trees_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(AGE.o, AGE, color = dens_age)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 200) + ylim(0, 200)
dev.off()

# calculate rmse different variables
rmse <- function(obs, est) sqrt(mean((obs - est)^2))

rmse_df <- data.frame(AGE = c(rmse(imp$AGE.o, imp$AGE),
                              rmse(imp$AGE.o, imp$AGE)/sd(imp$AGE.o)),
                      HT = c(rmse(imp$HT.o, imp$HT),
                             rmse(imp$HT.o, imp$HT)/sd(imp$HT.o)),
                      BA = c(rmse(imp$BA.o, imp$BA),
                             rmse(imp$BA.o, imp$BA)/sd(imp$BA.o)),
                      CC = c(rmse(imp$CC.o, imp$CC),
                             rmse(imp$CC.o, imp$CC)/sd(imp$CC.o)))

rownames(rmse_df) <- c('raw_rmsd', 'scaled_rmsd')

# write table
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_for_only/rmsd_10_90_no_ba_200_trees.csv')

# create df of WG
wg <- data.frame(obs = imp$WG.o,
                 est = imp$WG)

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
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/test_imp_for_only/accmat_wg_10_90_no_ba_200_trees.csv')

#########################################
###SCREENED FOR 10-90% RESIDUALS NO CC###
#########################################

###################################################
###LOAD POLYGON DATASET AND ADD LIDAR ATTRIBUTES###
###################################################

# load extracted data frame
# dat_fri_100 is loaded with LiDAR arributes derived from 
# masking pixels that do not fall 100 % inside polygons
load('D:/ontario_inventory/dat/dat_fri_100.RData')

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen_1bc_2apc_10perc_for.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# SPL vars only

# change all non-numeric variables to factor
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], 
                                         as.factor)

# build dataframe of reference variables
ref_vars <- c('cc', 'p95', 'avg', 'qav', 'cv', 'kur', 'max', 'ske')
x <- dat[, ref_vars]

# build dataframe of target variables
tar_vars <- c('HT', 'BA', 'WG')
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

# apply the imputation algorithm to the training/forest data
# we are not comparing to the other data because it doesn't 
# meet our criteria for clean forested polygons

# run imputation over new polygons
# we want to keep observed cols to see relationship
imp <- impute(object = rf,
              ancillaryData = anci,
              observed = T)

# compare input data and imputed data plots

# load get_density function
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# calculate density values
imp$dens_ht <- get_density(imp$HT.o, imp$HT, n = 50)
imp$dens_cc <- get_density(imp$CC.o, imp$CC, n = 50)
imp$dens_ba <- get_density(imp$BA.o, imp$BA, n = 50)
imp$dens_age <- get_density(imp$AGE.o, imp$AGE, n = 50)

# plot
png('D:/ontario_inventory/imputation/test_imp_for_only/ht_10_90_no_cc_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(HT.o, HT, color = dens_ht)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 35) + ylim(0, 35)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/cc_10_90_no_cc_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(CC.o, CC, color = dens_cc)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 100) + ylim(0, 100)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/ba_10_90_no_cc_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(BA.o, BA, color = dens_ba)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 60) + ylim(0, 60)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/age_10_90_no_cc_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(AGE.o, AGE, color = dens_age)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 200) + ylim(0, 200)
dev.off()

# calculate rmse different variables
rmse <- function(obs, est) sqrt(mean((obs - est)^2))

rmse_df <- data.frame(AGE = c(rmse(imp$AGE.o, imp$AGE),
                              rmse(imp$AGE.o, imp$AGE)/sd(imp$AGE.o)),
                      HT = c(rmse(imp$HT.o, imp$HT),
                             rmse(imp$HT.o, imp$HT)/sd(imp$HT.o)),
                      BA = c(rmse(imp$BA.o, imp$BA),
                             rmse(imp$BA.o, imp$BA)/sd(imp$BA.o)),
                      CC = c(rmse(imp$CC.o, imp$CC),
                             rmse(imp$CC.o, imp$CC)/sd(imp$CC.o)))

rownames(rmse_df) <- c('raw_rmsd', 'scaled_rmsd')

# write table
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_for_only/rmsd_10_90_no_cc.csv')

# create df of WG
wg <- data.frame(obs = imp$WG.o,
                 est = imp$WG)

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
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/test_imp_for_only/accmat_wg_10_90_no_cc.csv')

#########################################
###SCREENED FOR 30-70% RESIDUALS NO CC###
#########################################

###################################################
###LOAD POLYGON DATASET AND ADD LIDAR ATTRIBUTES###
###################################################

# load extracted data frame
# dat_fri_100 is loaded with LiDAR arributes derived from 
# masking pixels that do not fall 100 % inside polygons
load('D:/ontario_inventory/dat/dat_fri_100.RData')

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen_1bc_2apc_30perc_for.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# SPL vars only

# change all non-numeric variables to factor
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], 
                                         as.factor)

# build dataframe of reference variables
ref_vars <- c('cc', 'p95', 'avg', 'qav', 'cv', 'kur', 'max', 'ske')
x <- dat[, ref_vars]

# build dataframe of target variables
tar_vars <- c('HT', 'BA', 'WG')
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

# apply the imputation algorithm to the training/forest data
# we are not comparing to the other data because it doesn't 
# meet our criteria for clean forested polygons

# run imputation over new polygons
# we want to keep observed cols to see relationship
imp <- impute(object = rf,
              ancillaryData = anci,
              observed = T)

# compare input data and imputed data plots

# load get_density function
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# calculate density values
imp$dens_ht <- get_density(imp$HT.o, imp$HT, n = 50)
imp$dens_cc <- get_density(imp$CC.o, imp$CC, n = 50)
imp$dens_ba <- get_density(imp$BA.o, imp$BA, n = 50)
imp$dens_age <- get_density(imp$AGE.o, imp$AGE, n = 50)

# plot
png('D:/ontario_inventory/imputation/test_imp_for_only/ht_30_70_no_cc_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(HT.o, HT, color = dens_ht)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 35) + ylim(0, 35)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/cc_30_70_no_cc_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(CC.o, CC, color = dens_cc)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 100) + ylim(0, 100)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/ba_30_70_no_cc_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(BA.o, BA, color = dens_ba)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 60) + ylim(0, 60)
dev.off()

png('D:/ontario_inventory/imputation/test_imp_for_only/age_30_70_no_cc_plot.png', width = 480)
ggplot(imp) + 
  geom_point(aes(AGE.o, AGE, color = dens_age)) + 
  scale_color_viridis() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1) +
  xlim(0, 200) + ylim(0, 200)
dev.off()

# calculate rmse different variables
rmse <- function(obs, est) sqrt(mean((obs - est)^2))

rmse_df <- data.frame(AGE = c(rmse(imp$AGE.o, imp$AGE),
                              rmse(imp$AGE.o, imp$AGE)/sd(imp$AGE.o)),
                      HT = c(rmse(imp$HT.o, imp$HT),
                             rmse(imp$HT.o, imp$HT)/sd(imp$HT.o)),
                      BA = c(rmse(imp$BA.o, imp$BA),
                             rmse(imp$BA.o, imp$BA)/sd(imp$BA.o)),
                      CC = c(rmse(imp$CC.o, imp$CC),
                             rmse(imp$CC.o, imp$CC)/sd(imp$CC.o)))

rownames(rmse_df) <- c('raw_rmsd', 'scaled_rmsd')

# write table
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/test_imp_for_only/rmsd_30_70_no_cc.csv')

# create df of WG
wg <- data.frame(obs = imp$WG.o,
                 est = imp$WG)

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
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/test_imp_for_only/accmat_wg_30_70_no_cc.csv')


