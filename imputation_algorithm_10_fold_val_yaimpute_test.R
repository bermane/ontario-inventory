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

# create vector of imputation algorithms
methods <- c('euclidean',
             'raw',
             'mahalanobis',
             'ica',
             'msn',
             'msn2',
             'gnn',
             'randomForest')

# preallocate rmsd_df
rmsd_df <- data.frame(method = character(),
                      ht_raw = numeric(),
                      ht_scaled = numeric())

# iterate over methods
for(m in methods){
  # iterate over 10 folds
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
    ref_vars <- c('p95', 'cv')
    x <- train_data %>% select(p95, cv)
    
    # build dataframe of target variables
    tar_vars <- 'HT'
    y <- train_data %>% select('HT')
    
    # build dataframe of ancillary data
    anci <- train_data %>% select(HT)
    
    # run imputation algorithms
    if(m == 'randomForest') {
      imp <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y), k = 5)
    }else imp <- yai(x = x, y = y, method = m, k = 5)
    
    # impute predictions over test data
    pred <- predict(object = imp,
                    newdata = test_data,
                    ancillaryData = anci,
                    observed = T)
    
    # not sure why getting NAs for observed data
    # but can fill from test_data df
    pred$HT.o <- test_data$HT
    
    # calculate rmse different variables
    rmse <- function(obs, est) sqrt(mean((obs - est)^2))
    
    if(i == 1){
      
      # calculate rmse metrics
      perform_df <- data.frame(HT_raw_rmse = rmse(pred$HT.o, pred$HT),
                               HT_scaled_rmse = rmse(pred$HT.o, pred$HT)/sd(pred$HT.o))
      
      
      
    }else{
      
      # bind to perform_df
      perform_df <- rbind(perform_df,
                          data.frame(HT_raw_rmse = rmse(pred$HT.o, pred$HT),
                                     HT_scaled_rmse = rmse(pred$HT.o, pred$HT)/sd(pred$HT.o)))
    }
  }
  
  # add to final df
  rmsd_df <- rbind(rmsd_df,
                   data.frame(method = m,
                              ht_raw = colMeans(perform_df)[1],
                              ht_scaled = colMeans(perform_df)[2]))
}


# write as csv
write.csv(rmsd_df, 
          file = 'D:/ontario_inventory/imputation/distributions_for_only/rmsd_10_fold_yaimpute_p95_cv_ht.csv')
