# This code runs an imputation algorithm on LiDAR attributes connecting them to interpreter
# derived response attributes resulting in a model that can be applied to forest stands
# segmented from LiDAR rasters

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
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen_1b_2_10_90.csv')

# subset dat based on screening
dat <- dat[dat$POLYID %in% dat_screen$POLYID,]

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# change all non-numeric variables to factor
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], 
                                       as.factor)

#load names of lidar columns added to data above
colnames(dat[,112:127])

# PART A
# build list of reference variable combinations to try
ref_vars <- list(c('cc', 'p80', 'avg', 'qav', 'cv', 'kur', 'max', 'ske'),
                 c('cc', 'avg', 'qav', 'cv', 'max', 'ske', 'kur'),
                 c('cc', 'avg', 'cv', 'max'),
                 c('cc', 'cv', 'ske', 'max'),
                 c('cc', 'cv', 'kur'),
                 c('avg', 'cv', 'max', 'ske'),
                 c('cc', 'cv', 'max', 'ske'),
                 c('agb', 'ba', 'dens', 'lor', 'qmdbh', 'top_height', 'v', 'v_merch'),
                 c('ba', 'dens', 'lor', 'qmdbh'),
                 c('agb', 'ba', 'dens', 'lor', 'qmdbh'),
                 c('agb', 'ba', 'dens', 'lor', 'qmdbh', 'top_height'),
                 c('agb', 'ba', 'dens', 'lor', 'qmdbh', 'top_height', 'v'),
                 c('ba', 'lor', 'qmdbh'),
                 c('cc', 'p80', 'avg', 'qav', 'cv', 'kur', 'max', 'ske', 'agb', 'ba', 'dens', 'lor', 'qmdbh', 'top_height', 'v', 'v_merch'),
                 c('cc', 'cv', 'max', 'ske', 'ba', 'lor', 'qmdbh'),
                 c('cc', 'avg', 'cv', 'max', 'ba', 'dens', 'lor', 'qmdbh'),
                 c('cc', 'cv', 'max', 'ske', 'ba', 'dens', 'lor', 'qmdbh'),
                 c('cc', 'avg', 'cv', 'max', 'agb', 'ba', 'dens', 'lor', 'qmdbh', 'top_height'),
                 c('cc', 'avg', 'cv', 'max', 'agb', 'ba', 'dens', 'lor', 'qmdbh', 'top_height', 'v'),
                 c('cc', 'avg', 'qav', 'cv', 'max', 'ske', 'kur', 'agb', 'ba', 'dens', 'lor', 'qmdbh', 'top_height', 'v'),
                 c('cc', 'avg', 'qav', 'cv', 'max', 'ske', 'kur', 'agb', 'ba', 'dens', 'lor', 'qmdbh'))

# build list of target variable combinations to try
tar_vars <- list(c('HT', 'CC', 'BA', 'AGE', 'POLYTYPE'),
                 c('HT', 'CC', 'BA', 'POLYTYPE'),
                 c('HT', 'BA', 'POLYTYPE'),
                 c('HT', 'POLYTYPE'))

# allocate output lists
rmsd <- list()
rmsd_s <- list()

# start tick and timer
tick <- 1
tictoc::tic()

# run imputation algorithm over list of variable combinations
for(i in 1:length(ref_vars)){
  for(j in 1:length(tar_vars)){
    
    # set reference and target variables
    x <- dat[, ref_vars[[i]]]
    y <- dat[, tar_vars[[j]]]
    
    # run imputation algorithm
    rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y))
    
    # calculate RMSD for each variable
    rmsd[[tick]] <- rmsd(rf, vars = yvars(rf))
    
    # calculate scaled RMSD for model
    rmsd_s_hold <- rmsd(rf, vars = yvars(rf), scale = T)
    rmsd_s[[tick]] <- as.numeric(rmsd_s_hold[,1]) %>% mean(na.rm = T)
    
    # advance tick
    tick <- tick + 1
  }
}

# stop timer
tictoc::toc()

# put results into table
ref_df <- sapply(ref_vars, FUN = function(x){paste(x, collapse = ', ')}) %>% rep(each = 4)
tar_df <- sapply(tar_vars, FUN = function(x){paste(x, collapse = ', ')}) %>% rep(times = 21)

var_compare <- data.frame(ref_vars = ref_df,
                          tar_vars = tar_df,
                          rmsd_mean_scaled = sapply(rmsd_s, FUN = function(x) x))

# output csv
write.csv(var_compare, file = 'D:/ontario_inventory/imputation/test_results/ref_compare.csv', row.names = F)

# PART B
# now compare different combinations of target variables
ref_vars <- list(c('cc', 'p80', 'avg', 'qav', 'cv', 'kur', 'max', 'ske', 'agb', 'ba', 'dens', 'lor', 'qmdbh', 'top_height', 'v', 'v_merch'),
                 c('cc', 'avg', 'qav', 'cv', 'max', 'ske', 'kur', 'agb', 'ba', 'dens', 'lor', 'qmdbh'),
                 c('cc', 'avg', 'cv', 'max', 'agb', 'ba', 'dens', 'lor', 'qmdbh', 'top_height'))

# build list of target variable combinations to try
tar_vars <- list(c('HT', 'CC', 'BA', 'AGE', 'POLYTYPE'),
                 c('HT', 'CC', 'BA', 'AGE'),
                 c('HT', 'CC', 'BA', 'POLYTYPE'),
                 c('HT', 'CC', 'AGE', 'POLYTYPE'),
                 c('HT', 'BA', 'AGE', 'POLYTYPE'),
                 c('CC', 'BA', 'AGE', 'POLYTYPE'),
                 c('HT', 'CC', 'BA'),
                 c('HT', 'CC', 'POLYTYPE'),
                 c('HT', 'BA', 'POLYTYPE'),
                 c('CC', 'BA', 'POLYTYPE'),
                 c('HT', 'BA'),
                 c('HT', 'CC'),
                 'HT',
                 'CC',
                 'BA',
                 'AGE')

# allocate output lists
rmsd <- list()
rmsd_s <- list()

# start tick and timer
tick <- 1
tictoc::tic()

# run imputation algorithm over list of variable combinations
for(i in 1:length(ref_vars)){
  for(j in 1:length(tar_vars)){
    
    # set reference and target variables
    x <- dat[, ref_vars[[i]]]
    y <- dat[, tar_vars[[j]]]
    
    # run imputation algorithm
    rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y))
    
    # calculate RMSD for each variable
    rmsd[[tick]] <- rmsd(rf, vars = yvars(rf))
    
    # calculate scaled RMSD for model
    rmsd_s_hold <- rmsd(rf, vars = yvars(rf), scale = T)
    rmsd_s[[tick]] <- data.frame(rmsd_s = as.numeric(rmsd_s_hold[,1]) %>% mean(na.rm = T),
                                 ref_vars = paste(ref_vars[[i]], collapse = ', '),
                                 tar_vars = paste(tar_vars[[j]], collapse = ', '))
    
    # advance tick
    tick <- tick + 1
  }
}

# stop timer
tictoc::toc()

# put results into table
var_compare <- data.frame()

for(l in 1:length(rmsd_s)){
  var_compare <- rbind(var_compare, rmsd_s[[l]])
}

# output csv
write.csv(var_compare, file = 'D:/ontario_inventory/imputation/test_results/tar_compare.csv', row.names = F)

# PART C
# Compare imputation methods
ref_vars <- list(c('cc', 'p80', 'avg', 'qav', 'cv', 'kur', 'max', 'ske', 'agb', 'ba', 'dens', 'lor', 'qmdbh', 'top_height', 'v', 'v_merch'),
                 c('cc', 'avg', 'qav', 'cv', 'max', 'ske', 'kur', 'agb', 'ba', 'dens', 'lor', 'qmdbh'))

# build list of target variable combinations to try
tar_vars <- list(c('HT'),
                 c('HT', 'BA'),
                 c('HT', 'CC'),
                 c('HT', 'CC', 'BA'),
                 c('HT', 'CC', 'BA', 'AGE'))


# build list of imputation methods to try
imp_method <- c('euclidean',
                 'mahalanobis',
                 'msn',
                 'gnn',
                 'ica',
                 'randomForest')

# allocate output lists
rmsd <- list()
rmsd_s <- list()

# start tick and timer
tick <- 1
tictoc::tic()

# run imputation algorithm over list of variable combinations and imputation methods
for(meth in 1:length(imp_method)){
  for(i in 1:length(ref_vars)){
    for(j in 1:length(tar_vars)){
      
      # set reference and target variables
      x <- dat[, ref_vars[[i]]]
      y <- dat[, tar_vars[[j]]]
      
      # run imputation algorithm
      if(imp_method[meth] == 'randomForest'){
        rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y))
      }else{
        rf <- yai(x = x, y = y, method = imp_method[meth])
      }

      # calculate RMSD for each variable
      rmsd[[tick]] <- rmsd(rf, vars = yvars(rf))
      
      # calculate scaled RMSD for model
      rmsd_s_hold <- rmsd(rf, vars = yvars(rf), scale = T)
      rmsd_s[[tick]] <- data.frame(rmsd_s = as.numeric(rmsd_s_hold[,1]) %>% mean(na.rm = T),
                                   ref_vars = paste(ref_vars[[i]], collapse = ', '),
                                   tar_vars = paste(tar_vars[[j]], collapse = ', '),
                                   imp_method = imp_method[meth])
      
      # advance tick
      tick <- tick + 1
    }
  }
}

# stop timer
tictoc::toc()

# put results into table
var_compare <- data.frame()

for(l in 1:length(rmsd_s)){
  var_compare <- rbind(var_compare, rmsd_s[[l]])
}

# output csv
write.csv(var_compare, file = 'D:/ontario_inventory/imputation/test_results/imp_method_compare.csv', row.names = F)

# PART D
# Test different data screening procedures

# clear workspace
rm(list=ls())

# load extracted data frame
load('D:/ontario_inventory/imputation/imputation_df.RData')

# load final datasets after screening
dat_screen_10_90 <- read.csv('D:/ontario_inventory/imputation/dat_screen_2_10_90.csv')
dat_screen_20_80 <- read.csv('D:/ontario_inventory/imputation/dat_screen_2_20_80.csv')
dat_screen_30_70 <- read.csv('D:/ontario_inventory/imputation/dat_screen_2_30_70.csv')
dat_screen_25_75 <- read.csv('D:/ontario_inventory/imputation/dat_screen_2_25_75.csv')

# combine into list object
dat_screen_l <- list(dat_screen_10_90, dat_screen_20_80,
                     dat_screen_30_70, dat_screen_25_75)

# create var with screen names
dat_sc_name <- c('dat_screen_10_90', 'dat_screen_20_80',
                  'dat_screen_30_70', 'dat_screen_25_75')

# load ref vars
ref_vars <- list(c('cc', 'p80', 'avg', 'qav', 'cv', 'kur', 'max', 'ske', 'agb', 'ba', 'dens', 'lor', 'qmdbh', 'top_height', 'v', 'v_merch'),
                 c('cc', 'avg', 'qav', 'cv', 'max', 'ske', 'kur', 'agb', 'ba', 'dens', 'lor', 'qmdbh'))

# load tar vars
tar_vars <- list(c('HT'),
                 c('HT', 'BA'),
                 c('HT', 'CC'),
                 c('HT', 'CC', 'BA'),
                 c('HT', 'CC', 'BA', 'AGE'),
                 c('HT', 'POLYTYPE'),
                 c('HT', 'BA', 'POLYTYPE'),
                 c('HT', 'CC', 'POLYTYPE'),
                 c('HT', 'CC', 'BA', 'POLYTYPE'),
                 c('HT', 'CC', 'BA', 'AGE', 'POLYTYPE'))

# allocate output lists
rmsd <- list()
rmsd_s <- list()

# start tick and timer
tick <- 1
tictoc::tic()

# loop through list and run models
for(s in 1:length(dat_screen_l)){
  
  # subset dat based on screening
  dat_sc <- dat[dat$POLYID %in% dat_screen_l[[s]]$POLYID,]
  
  # change all non-numeric variables to factor
  dat_sc[sapply(dat_sc, is.character)] <- lapply(dat_sc[sapply(dat_sc, is.character)], 
                                           as.factor)
  
  # run imputation algorithm over list of variable combinations
  for(i in 1:length(ref_vars)){
    for(j in 1:length(tar_vars)){
      
      # set reference and target variables
      x <- dat_sc[, ref_vars[[i]]]
      y <- dat_sc[, tar_vars[[j]]]
      
      # run imputation algorithm
      rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y))
      
      # calculate RMSD for each variable
      rmsd[[tick]] <- rmsd(rf, vars = yvars(rf))
      
      # calculate scaled RMSD for model
      rmsd_s_hold <- rmsd(rf, vars = yvars(rf), scale = T)
      rmsd_s[[tick]] <- data.frame(rmsd_s = as.numeric(rmsd_s_hold[,1]) %>% mean(na.rm = T),
                                   ref_vars = paste(ref_vars[[i]], collapse = ', '),
                                   tar_vars = paste(tar_vars[[j]], collapse = ', '),
                                   dat_screen = dat_sc_name[s])
      
      # advance tick
      tick <- tick + 1
    }
  }
}

# stop timer
tictoc::toc()

# put results into table
var_compare <- data.frame()

for(l in 1:length(rmsd_s)){
  var_compare <- rbind(var_compare, rmsd_s[[l]])
}

# output csv
write.csv(var_compare, file = 'D:/ontario_inventory/imputation/test_results/dat_screen_compare.csv', row.names = F)

# PART E
# Test different number of trees in random forest algorithm

# clear workspace
rm(list=ls())

# load extracted data frame
load('D:/ontario_inventory/imputation/imputation_df.RData')

# load final datasets after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/dat_screen_2_30_70.csv')

# subset dat based on screening
dat_sc <- dat[dat$POLYID %in% dat_screen$POLYID,]

# change all non-numeric variables to factor
dat_sc[sapply(dat_sc, is.character)] <- lapply(dat_sc[sapply(dat_sc, is.character)], 
                                               as.factor)

# load ref vars
ref_vars <- list(c('cc', 'p80', 'avg', 'qav', 'cv', 'kur', 'max', 'ske', 'agb', 'ba', 'dens', 'lor', 'qmdbh', 'top_height', 'v', 'v_merch'))

# load tar vars
tar_vars <- list(c('HT', 'POLYTYPE'),
                 c('HT', 'BA', 'POLYTYPE'),
                 c('HT', 'CC', 'POLYTYPE'),
                 c('HT', 'CC', 'BA', 'POLYTYPE'))

# set number of different iterations of trees
tree_iter <- c('50%', '75%', '90%', '100%',
               '110%', '125%', '150%')

# allocate output lists
rmsd <- list()
rmsd_s <- list()

# start tick and timer
tick <- 1
tictoc::tic()
  
# run imputation algorithm over list of variable combinations
for(i in 1:length(ref_vars)){
  for(j in 1:length(tar_vars)){
    for(iter in 1:length(tree_iter)){
      
      # set reference and target variables
      x <- dat_sc[, ref_vars[[i]]]
      y <- dat_sc[, tar_vars[[j]]]
      
      # run imputation algorithm based on tree iter
      if(iter == 1) rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y)*0.5)
      if(iter == 2) rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y)*0.75)
      if(iter == 3) rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y)*0.9)
      if(iter == 4) rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y))
      if(iter == 5) rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y)*1.1)
      if(iter == 6) rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y)*1.25)
      if(iter == 7) rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y)*1.5)
      
      # calculate RMSD for each variable
      rmsd[[tick]] <- rmsd(rf, vars = yvars(rf))
      
      # calculate scaled RMSD for model
      rmsd_s_hold <- rmsd(rf, vars = yvars(rf), scale = T)
      rmsd_s[[tick]] <- data.frame(rmsd_s = as.numeric(rmsd_s_hold[,1]) %>% mean(na.rm = T),
                                   ref_vars = paste(ref_vars[[i]], collapse = ', '),
                                   tar_vars = paste(tar_vars[[j]], collapse = ', '),
                                   num_trees = tree_iter[iter])
      
      # advance tick
      tick <- tick + 1
    }
  }
}
    
# stop timer
tictoc::toc()

# put results into table
var_compare <- data.frame()

for(l in 1:length(rmsd_s)){
  var_compare <- rbind(var_compare, rmsd_s[[l]])
}

# output csv
write.csv(var_compare, file = 'D:/ontario_inventory/imputation/test_results/num_trees_compare.csv', row.names = F)
