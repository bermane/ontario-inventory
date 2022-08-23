# This code runs an imputation algorithm on LiDAR attributes connecting them to interpreter
# derived response attributes resulting in a model that can be applied to forest stands
# segmented from LiDAR rasters

# This version of the code uses a simple nearest neighbor approach to query
# the attributes of choice and find the sample that minimizes euclidean
# distance between the variables

# we can then impute all the attributes of the nearest neighbor into the 
# newly generated polygon

# load packages
library(terra)
library(tidyverse)
library(magrittr)
library(RANN)

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

# # load photo interpreted polygons
# poly <-
#   vect(
#     'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp'
#   )

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <-
  read.csv(
    'D:/ontario_inventory/imputation/fri_data_screen_1bc_2a_2p95_2cc_10perc_2022_07_06.csv'
  )

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID, ]

# remove any missing values
dat <- na.omit(dat)

#####################################
###RUN NEAREST NEIGHBOR IMPUTATION###
#####################################

#####################
###P95, CC, and CV###
#####################

# lets try using p95, cc and cv, just like the segmentation
dat_nn <- dat %>% select(p95, cc, cv)

# scale for nn computation
dat_nn_scaled <- dat_nn %>% scale

# run nearest neighbor
nn <- nn2(dat_nn_scaled, dat_nn_scaled, k = 2)

# create tibble of results
nn_tab <- tibble(p95 = dat_nn[,1],
                 cc = dat_nn[,2],
                 cv = dat_nn[,3],
                 age = dat$AGE,
                 p95_nn = dat_nn[nn[[1]][,2],1],
                 cc_nn = dat_nn[nn[[1]][,2],2],
                 cv_nn = dat_nn[nn[[1]][,2],3],
                 age_nn = dat$AGE[nn[[1]][,2]])

# calculate rmse different variables
rmse <- function(obs, est){
  sqrt(mean((obs - est) ^ 2))
}

# calculate rmse metrics
perform_df <-
  data.frame(
    p95_raw_rmse = rmse(nn_tab$p95, nn_tab$p95_nn),
    cc_raw_rmse = rmse(nn_tab$cc, nn_tab$cc_nn),
    cv_raw_rmse = rmse(nn_tab$cv, nn_tab$cv_nn),
    age_raw_rmse = rmse(nn_tab$age, nn_tab$age_nn),
    p95_scaled_rmse = rmse(nn_tab$p95, nn_tab$p95_nn) /
      sd(nn_tab$p95),
    cc_scaled_rmse = rmse(nn_tab$cc, nn_tab$cc_nn) /
      sd(nn_tab$cc),
    cv_scaled_rmse = rmse(nn_tab$cv, nn_tab$cv_nn) /
      sd(nn_tab$cv),
    age_scaled_rmse = rmse(nn_tab$age, nn_tab$age_nn) /
      sd(nn_tab$age)
  )

###########################
###ALS HT, CC, BA vs FRI###
###########################

# create df for als and fri metrics
dat_als <- dat %>% select(lor, cc, ba)
dat_fri <- dat %>% select(HT, CC, BA)

# change als colnames
colnames(dat_als) <- c('HT', 'CC', 'BA')

# scale for nn computation
# need to combine and scale all values together then separate again
dat_alsfri_scaled <- rbind(dat_als, dat_fri) %>% scale
dat_als_scaled <- dat_alsfri_scaled[1:NROW(dat),]
dat_fri_scaled <- dat_alsfri_scaled[(NROW(dat)+1):(NROW(dat)*2),]

# run nearest neighbor
nn <- nn2(dat_fri_scaled, dat_als_scaled, k = 2)

# create tibble of results
nn_tab <- tibble(ht = dat_als[,1],
                 cc = dat_als[,2],
                 ba = dat_als[,3],
                 age = dat$AGE,
                 ht_nn = dat_fri[nn[[1]][,1],1],
                 cc_nn = dat_fri[nn[[1]][,1],2],
                 ba_nn = dat_fri[nn[[1]][,1],3],
                 age_nn = dat$AGE[nn[[1]][,1]])

# calculate rmse metrics
perform_df <-
  data.frame(
    ht_raw_rmse = rmse(nn_tab$ht, nn_tab$ht_nn),
    cc_raw_rmse = rmse(nn_tab$cc, nn_tab$cc_nn),
    ba_raw_rmse = rmse(nn_tab$ba, nn_tab$ba_nn),
    age_raw_rmse = rmse(nn_tab$age, nn_tab$age_nn),
    ht_scaled_rmse = rmse(nn_tab$ht, nn_tab$ht_nn) /
      sd(nn_tab$ht),
    cc_scaled_rmse = rmse(nn_tab$cc, nn_tab$cc_nn) /
      sd(nn_tab$cc),
    ba_scaled_rmse = rmse(nn_tab$ba, nn_tab$ba_nn) /
      sd(nn_tab$ba),
    age_scaled_rmse = rmse(nn_tab$age, nn_tab$age_nn) /
      sd(nn_tab$age)
  )


library(RANN)

set.seed(20)
r <- seq(1,100,1)
q <- seq(1,100,1) %>% sample

nn <- nn2(r, q, k = 1)

check <- tibble(obs = q,
                est = r[nn[[1]][,1]])
