# This code runs an imputation algorithm on LiDAR attributes connecting them to interpreter
# derived response attributes resulting in a model that can be applied to forest stands
# segmented from LiDAR rasters

# The purpose is to compare the distributions of various attributes (HT, CC, BA, AGE, POLYTYPE, WG)
# from FRI polygons, LiDAR attributes in FRI polygons, imputed vars into FRI polys, imputed vars into LiDAR polys

# this new version only runs forested pixels since we are only really interested in modelling
# and imputing attributes in forests
# it uses data screened for 10 and 90th percentile of residuals

# This version of the code runs a 10 fold cross validation on the screened data 
# to assess model performance

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

# # save full original dataset
# dat_orig <- dat
# 
# # only keep forested polygons
# dat_orig <- filter(dat_orig, POLYTYPE == 'FOR')

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/fri_data_screen_1bc_2a_2p95_2cc_10perc_2022_07_06.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# set seed so reproducible
set.seed(14)

# split into 10 folds
# randomly shuffle the data
dat_fold <- dat[sample(nrow(dat)), ]

# change all non-numeric variables to factor
dat_fold[sapply(dat_fold, is.character)] <-
  lapply(dat_fold[sapply(dat_fold, is.character)],
         as.factor)

# create 10 equally size folds
folds <- cut(seq(1, nrow(dat_fold)), breaks = 10, labels = FALSE)

# iterate over 10 folds
for(i in 1:10){
  # segment data by fold using the which() function
  test_index <- which(folds == i, arr.ind = TRUE)
  test_data <- dat_fold[test_index,]
  train_data <- dat_fold[-test_index,]
  
  # als data only
  # build dataframe of reference variables
  ref_vars <-
    c('cc', 'p95', 'avg', 'qav', 'cv', 'kur', 'max', 'ske')
  x <- train_data[, ref_vars]
  x <- rbind(x, test_data[, ref_vars])
  
  # build dataframe of target variables
  tar_vars <- c('HT', 'WG')
  y <- train_data[, tar_vars]
  
  # build dataframe of ancillary data
  anci <-
    train_data[, c('HT', 'CC', 'BA', 'WG', 'AGE', 'p95', 'cc', 'cv', 'avg', 'max')]
  
  # run imputation algorithms
  imp <-
    yai(
      x = x,
      y = y,
      method = 'randomForest',
      ntree = 100 * NCOL(y)
    )
  
  # regular RF data
  dat_rf <- train_data[, c(ref_vars, 'WG')]
  
  # run random forest
  rf <- randomForest(WG ~ ., data = dat_rf)
  
  # predict rf
  pred_rf <- predict(rf, newdata = test_data)
  
  # unsupervised
  # imp <-
  #   yai(
  #     x = x,
  #     method = 'randomForest',
  #     ntree = 100 * NCOL(x)
  #   )
  
  # impute predictions over test data
  # pred <- predict(
  #   object = imp,
  #   newdata = test_data,
  #   ancillaryData = anci,
  #   observed = T
  # )
  
  # impute predictions over test data
  pred <- impute(imp, ancillaryData = anci)
  pred <- pred[rownames(pred) %in% rownames(test_data),]
  
  # # predictions for multiple k
  # pred <- predict(
  #   object = imp,
  #   newdata = test_data,
  #   ancillaryData = anci,
  #   observed = T,
  #   method = 'median',
  #   method.factor = 'mean or median'
  # )
  
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
  rmse <- function(obs, est){
    sqrt(mean((obs - est) ^ 2))
  }
  
  if (i == 1) {
    # calculate rmse metrics
    perform_df <-
      data.frame(
        AGE_raw_rmse = rmse(pred$AGE.o, pred$AGE),
        HT_raw_rmse = rmse(pred$HT.o, pred$HT),
        BA_raw_rmse = rmse(pred$BA.o, pred$BA),
        CC_raw_rmse = rmse(pred$CC.o, pred$CC),
        p95_raw_rmse = rmse(pred$p95.o, pred$p95),
        cc_raw_rmse = rmse(pred$cc.o, pred$cc),
        cv_raw_rmse = rmse(pred$cv.o, pred$cv),
        avg_raw_rmse = rmse(pred$avg.o, pred$avg),
        max_raw_rmse = rmse(pred$max.o, pred$max),
        AGE_scaled_rmse = rmse(pred$AGE.o, pred$AGE) /
          sd(pred$AGE.o),
        HT_scaled_rmse = rmse(pred$HT.o, pred$HT) /
          sd(pred$HT.o),
        BA_scaled_rmse = rmse(pred$BA.o, pred$BA) /
          sd(pred$BA.o),
        CC_scaled_rmse = rmse(pred$CC.o, pred$CC) /
          sd(pred$CC.o),
        p95_scaled_rmse = rmse(pred$p95.o, pred$p95) /
          sd(pred$p95.o),
        cc_scaled_rmse = rmse(pred$cc.o, pred$cc) /
          sd(pred$cc.o),
        cv_scaled_rmse = rmse(pred$cv.o, pred$cv) /
          sd(pred$cv.o),
        avg_scaled_rmse = rmse(pred$avg.o, pred$avg) /
          sd(pred$avg.o),
        max_scaled_rmse = rmse(pred$max.o, pred$max) /
          sd(pred$max.o)
      )
    
    # calculate wg accuracy
    # create df of WG
    wg <- data.frame(obs = pred$WG.o,
                     est = pred$WG)
    
    # make estimate levels match obs levels
    levels(wg$est) <-
      c(levels(wg$est), levels(wg$obs)[!(levels(wg$obs) %in% levels(wg$est))])
    wg$est <- factor(wg$est, levels = levels(wg$obs))
    
    # create column of match or not
    wg$match <- wg$obs == wg$est
    
    # add total percent of matching WG to perform_df
    perform_df <- cbind(perform_df,
                        data.frame(wg_accuracy = NROW(wg[wg$match == T, ]) /
                                     NROW(wg)))
    
  } else{
    # calculate wg accuracy
    # create df of WG
    wg <- data.frame(obs = pred$WG.o,
                     est = pred$WG)
    
    # make estimate levels match obs levels
    levels(wg$est) <-
      c(levels(wg$est), levels(wg$obs)[!(levels(wg$obs) %in% levels(wg$est))])
    wg$est <- factor(wg$est, levels = levels(wg$obs))
    
    # create column of match or not
    wg$match <- wg$obs == wg$est
    
    # bind to perform_df
    perform_df <- rbind(
      perform_df,
      data.frame(
        AGE_raw_rmse = rmse(pred$AGE.o, pred$AGE),
        HT_raw_rmse = rmse(pred$HT.o, pred$HT),
        BA_raw_rmse = rmse(pred$BA.o, pred$BA),
        CC_raw_rmse = rmse(pred$CC.o, pred$CC),
        p95_raw_rmse = rmse(pred$p95.o, pred$p95),
        cc_raw_rmse = rmse(pred$cc.o, pred$cc),
        cv_raw_rmse = rmse(pred$cv.o, pred$cv),
        avg_raw_rmse = rmse(pred$avg.o, pred$avg),
        max_raw_rmse = rmse(pred$max.o, pred$max),
        AGE_scaled_rmse = rmse(pred$AGE.o, pred$AGE) /
          sd(pred$AGE.o),
        HT_scaled_rmse = rmse(pred$HT.o, pred$HT) /
          sd(pred$HT.o),
        BA_scaled_rmse = rmse(pred$BA.o, pred$BA) /
          sd(pred$BA.o),
        CC_scaled_rmse = rmse(pred$CC.o, pred$CC) /
          sd(pred$CC.o),
        p95_scaled_rmse = rmse(pred$p95.o, pred$p95) /
          sd(pred$p95.o),
        cc_scaled_rmse = rmse(pred$cc.o, pred$cc) /
          sd(pred$cc.o),
        cv_scaled_rmse = rmse(pred$cv.o, pred$cv) /
          sd(pred$cv.o),
        avg_scaled_rmse = rmse(pred$avg.o, pred$avg) /
          sd(pred$avg.o),
        max_scaled_rmse = rmse(pred$max.o, pred$max) /
          sd(pred$max.o),
        wg_accuracy = NROW(wg[wg$match == T, ]) /
          NROW(wg)
      )
    )
  }
}
  
# summarize performance
perform_avg <- tibble(colnames(perform_df), colMeans(perform_df))

# write as csv
# write.csv(perform_avg, wfile = 'D:/ontario_inventory/imputation/ten_fold_validation/rmsd_10_fold_age_wg.csv')
