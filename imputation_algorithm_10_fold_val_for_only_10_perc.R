# This code runs an imputation algorithm on LiDAR attributes connecting them to interpreter
# derived response attributes resulting in a model that can be applied to forest stands
# segmented from LiDAR rasters

# The purpose is to compare the distributions of various attributes (HT, CC, BA, AGE, POLYTYPE, WG)
# from FRI polygons, LiDAR attributes in FRI polygons, imputed vars into FRI polys, imputed vars into LiDAR polys

# this new version only runs forested pixels since we are only really interested in modelling
# and imputing attributes in forests
# it uses data screened for 10 and 90th percentile of residuals

# This version of the code runs a 10 fold cross validation on the screened data 
# to better assess model performance using the best data

# load packages
library(terra)
library(tidyverse)
library(yaImpute)
library(randomForest)
library(viridis)

###################################################
###LOAD POLYGON DATASET AND ADD LIDAR ATTRIBUTES###
###################################################

# load extracted data frame
# dat_fri_100 is loaded with LiDAR arributes derived from 
# masking pixels that do not fall 100 % inside polygons
load('D:/ontario_inventory/dat/dat_fri_100.RData')

# save full original dataset
dat_orig <- dat

# only keep forested polygons
dat_orig <- filter(dat_orig, POLYTYPE == 'FOR')

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

# set seed so reproducible
set.seed(4)

# split into 10 folds
# randomly shuffle the data
dat_fold <- dat[sample(nrow(dat)),]

# create 10 equally size folds
folds <- cut(seq(1,nrow(dat_fold)),breaks=10,labels=FALSE)

# iterative over 10 folds
for(i in 1:10){
  
  # segment data by fold using the which() function 
  test_index <- which(folds==i,arr.ind=TRUE)
  test_data <- dat_fold[test_index, ]
  train_data <- dat_fold[-test_index, ]
  
  # change all non-numeric variables to factor
  train_data[sapply(train_data, is.character)] <- lapply(train_data[sapply(train_data, is.character)], 
                                           as.factor)
  
  # change all non-numeric variables to factor
  test_data[sapply(test_data, is.character)] <- lapply(test_data[sapply(test_data, is.character)], 
                                                         as.factor)
  
  # drop rows from test data that have levels not in train data
  test_data <- test_data[test_data$WG %in% levels(train_data$WG),]
  
  # make working group test levels match train levels
  levels(test_data$WG) <- c(levels(test_data$WG), levels(train_data$WG)[!(levels(train_data$WG) %in% levels(test_data$WG))])
  test_data$WG <- factor(test_data$WG, levels = levels(train_data$WG))
  
  # spl data only
  # build dataframe of reference variables
  ref_vars <- c('cc', 'p95', 'avg', 'qav', 'cv', 'kur', 'max', 'ske')
  x <- train_data[, ref_vars]
  
  # build dataframe of target variables
  tar_vars <- c('HT', 'WG')
  y <- train_data[, tar_vars]
  
  # build dataframe of ancillary data
  anci <- train_data[, c('HT', 'CC', 'BA', 'WG', 'AGE', 'p95', 'cc', 'cv', 'avg', 'max')]
  
  # run imputation algorithm
  rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y))
  
  # impute predictions over test data
  pred <- predict(object = rf,
                      newdata = test_data,
                  ancillaryData = anci,
                      observed = T)
  
  # not sure why getting NAs for observed data
  # but can fill from test_data df
  pred$HT.o <- test_data$HT
  pred$CC.o <- test_data$CC
  pred$BA.o <- test_data$BA
  pred$WG.o <- test_data$WG
  pred$AGE.o <- test_data$AGE
  pred$p95.o <- test_data$p95
  pred$cc.o <- test_data$cc
  pred$cv.o <- test_data$cv
  pred$avg.o <- test_data$avg
  pred$max.o <- test_data$max

  # calculate rmse different variables
  rmse <- function(obs, est) sqrt(mean((obs - est)^2))
  
  if(i == 1){
    
    # calculate rmse metrics
    perform_df <- data.frame(AGE_raw_rmse = rmse(pred$AGE.o, pred$AGE),
                        HT_raw_rmse = rmse(pred$HT.o, pred$HT),
                        BA_raw_rmse = rmse(pred$BA.o, pred$BA),
                        CC_raw_rmse = rmse(pred$CC.o, pred$CC),
                        p95_raw_rmse = rmse(pred$p95.o, pred$p95),
                        cc_raw_rmse = rmse(pred$cc.o, pred$cc),
                        cv_raw_rmse = rmse(pred$cv.o, pred$cv),
                        avg_raw_rmse = rmse(pred$avg.o, pred$avg),
                        max_raw_rmse = rmse(pred$max.o, pred$max),
                        AGE_scaled_rmse = rmse(pred$AGE.o, pred$AGE)/sd(pred$AGE.o),
                        HT_scaled_rmse = rmse(pred$HT.o, pred$HT)/sd(pred$HT.o),
                        BA_scaled_rmse = rmse(pred$BA.o, pred$BA)/sd(pred$BA.o),
                        CC_scaled_rmse = rmse(pred$CC.o, pred$CC)/sd(pred$CC.o),
                        p95_scaled_rmse = rmse(pred$p95.o, pred$p95)/sd(pred$p95.o),
                        cc_scaled_rmse = rmse(pred$cc.o, pred$cc)/sd(pred$cc.o),
                        cv_scaled_rmse = rmse(pred$cv.o, pred$cv)/sd(pred$cv.o),
                        avg_scaled_rmse = rmse(pred$avg.o, pred$avg)/sd(pred$avg.o),
                        max_scaled_rmse = rmse(pred$max.o, pred$max)/sd(pred$max.o))
    
    # calculate wg accuracy
    # create df of WG
    wg <- data.frame(obs = pred$WG.o,
                     est = pred$WG)
    
    # make estimate levels match obs levels
    levels(wg$est) <- c(levels(wg$est), levels(wg$obs)[!(levels(wg$obs) %in% levels(wg$est))])
    wg$est <- factor(wg$est, levels = levels(wg$obs))
    
    # create column of match or not
    wg$match <- wg$obs == wg$est
    
    # add total percent of matching WG to perform_df
    perform_df <- cbind(perform_df,
                        data.frame(wg_accuracy = NROW(wg[wg$match == T,])/NROW(wg)))
    
  }else{

    # calculate wg accuracy
    # create df of WG
    wg <- data.frame(obs = pred$WG.o,
                     est = pred$WG)
    
    # make estimate levels match obs levels
    levels(wg$est) <- c(levels(wg$est), levels(wg$obs)[!(levels(wg$obs) %in% levels(wg$est))])
    wg$est <- factor(wg$est, levels = levels(wg$obs))
    
    # create column of match or not
    wg$match <- wg$obs == wg$est
    
    # bind to perform_df
    perform_df <- rbind(perform_df,
                        data.frame(AGE_raw_rmse = rmse(pred$AGE.o, pred$AGE),
                                   HT_raw_rmse = rmse(pred$HT.o, pred$HT),
                                   BA_raw_rmse = rmse(pred$BA.o, pred$BA),
                                   CC_raw_rmse = rmse(pred$CC.o, pred$CC),
                                   p95_raw_rmse = rmse(pred$p95.o, pred$p95),
                                   cc_raw_rmse = rmse(pred$cc.o, pred$cc),
                                   cv_raw_rmse = rmse(pred$cv.o, pred$cv),
                                   avg_raw_rmse = rmse(pred$avg.o, pred$avg),
                                   max_raw_rmse = rmse(pred$max.o, pred$max),
                                   AGE_scaled_rmse = rmse(pred$AGE.o, pred$AGE)/sd(pred$AGE.o),
                                   HT_scaled_rmse = rmse(pred$HT.o, pred$HT)/sd(pred$HT.o),
                                   BA_scaled_rmse = rmse(pred$BA.o, pred$BA)/sd(pred$BA.o),
                                   CC_scaled_rmse = rmse(pred$CC.o, pred$CC)/sd(pred$CC.o),
                                   p95_scaled_rmse = rmse(pred$p95.o, pred$p95)/sd(pred$p95.o),
                                   cc_scaled_rmse = rmse(pred$cc.o, pred$cc)/sd(pred$cc.o),
                                   cv_scaled_rmse = rmse(pred$cv.o, pred$cv)/sd(pred$cv.o),
                                   avg_scaled_rmse = rmse(pred$avg.o, pred$avg)/sd(pred$avg.o),
                                   max_scaled_rmse = rmse(pred$max.o, pred$max)/sd(pred$max.o),
                                   wg_accuracy = NROW(wg[wg$match == T,])/NROW(wg)))
    }
  }

# calculate columnMeans
perform_avg <- colMeans(perform_df)

# write as csv
write.csv(perform_avg, 
          file = 'D:/ontario_inventory/imputation/distributions_for_only/rmsd_10_fold_ht_wg_10_perc.csv')
