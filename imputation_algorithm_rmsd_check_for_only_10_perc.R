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
tar_vars <- c('HT', 'BA', 'CC', 'WG')
y <- dat[, tar_vars]

# build list of target variable combinations to try
tar_list <- combn(tar_vars, 1, simplify = F)
tar_list <- append(tar_list, combn(tar_vars, 2, simplify = F))
tar_list <- append(tar_list, combn(tar_vars, 3, simplify = F))
tar_list <- append(tar_list, combn(tar_vars, 4, simplify = F))

# build dataframe of ancillary data
anci <- dat[, !(colnames(dat) %in% ref_vars)]

# apply imputation over list and generate dataframe of results
for(i in 1:length(tar_list)){
  y1 <- y %>% select(tar_list[[i]])
  
  # run imputation algorithm
  rf <- yai(x = x, y = y1, method = 'randomForest', ntree = 100*NCOL(y))
  
  # impute nearest neighbors
  imp <- impute(object = rf,
                ancillaryData = anci,
                observed = T)
  
  # calculate accuracy of WG
  # create df of WG
  wg <- data.frame(obs = imp$WG.o,
                   est = imp$WG)
  
  # make estimate levels match obs levels
  levels(wg$est) <- c(levels(wg$est), levels(wg$obs)[!(levels(wg$obs) %in% levels(wg$est))])
  wg$est <- factor(wg$est, levels = levels(wg$obs))
  
  # create column of match or not
  wg$match <- wg$obs == wg$est
  
  # total percent of matching WG
  wg_acc <- NROW(wg[wg$match == T,])/NROW(wg)
  
  # calculate rmse different variables
  rmse <- function(obs, est) sqrt(mean((obs - est)^2))
  
  if(i == 1){
    rmse_df <- data.frame(NAME = toString(tar_list[[i]]),
                          RMSD = c('raw', 'scaled'),
                          AGE = c(rmse(imp$AGE.o, imp$AGE),
                                  rmse(imp$AGE.o, imp$AGE)/sd(imp$AGE.o)),
                          HT = c(rmse(imp$HT.o, imp$HT),
                                 rmse(imp$HT.o, imp$HT)/sd(imp$HT.o)),
                          BA = c(rmse(imp$BA.o, imp$BA),
                                 rmse(imp$BA.o, imp$BA)/sd(imp$BA.o)),
                          CC = c(rmse(imp$CC.o, imp$CC),
                                 rmse(imp$CC.o, imp$CC)/sd(imp$CC.o)),
                          WG_ACC = wg_acc)
  }else{
    rmse_df <- rbind(rmse_df, data.frame(NAME = toString(tar_list[[i]]),
                          RMSD = c('raw', 'scaled'),
                          AGE = c(rmse(imp$AGE.o, imp$AGE),
                                  rmse(imp$AGE.o, imp$AGE)/sd(imp$AGE.o)),
                          HT = c(rmse(imp$HT.o, imp$HT),
                                 rmse(imp$HT.o, imp$HT)/sd(imp$HT.o)),
                          BA = c(rmse(imp$BA.o, imp$BA),
                                 rmse(imp$BA.o, imp$BA)/sd(imp$BA.o)),
                          CC = c(rmse(imp$CC.o, imp$CC),
                                 rmse(imp$CC.o, imp$CC)/sd(imp$CC.o)),
                          WG_ACC = wg_acc))
  }

}

# separate raw and scaled
raw <- rmse_df %>% filter(RMSD == 'raw')
scaled <- rmse_df %>% filter(RMSD == 'scaled')

# export datasets
# write to csv
write.csv(raw, file = 'D:/ontario_inventory/imputation/distributions_for_only/rmsd_compare_10_perc_raw.csv')
write.csv(scaled, file = 'D:/ontario_inventory/imputation/distributions_for_only/rmsd_compare_10_perc_scaled.csv')
