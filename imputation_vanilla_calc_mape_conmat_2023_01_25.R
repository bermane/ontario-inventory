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
library(corrplot)
library(reshape2)
library(gtools)
library(janitor)
library(doParallel)

###################################################
###LOAD POLYGON DATASET AND ADD LIDAR ATTRIBUTES###
###################################################

# load extracted data frame
load('D:/ontario_inventory/dat/dat_fri_extr.RData')

# load photo interpreted polygons
poly <-
  vect(
    'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp'
  )

# cbind centroids to dat
dat <- cbind(dat, centroids(poly) %>% crds)

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# # load final dataset after screening
# # old with HT LM
# dat_screen <-
#   read.csv(
#     'D:/ontario_inventory/imputation/fri_data_screen_1bc_2a_2p95_2cc_10perc_2022_07_06.csv'
#   )

# # load final dataset after screening
# # without HT lm and with cv 0.1
# dat_screen <-
#   read.csv(
#     'D:/ontario_inventory/imputation/dat_screen_1b_01_1c_2pc_for.csv'
#   )

# load final dataset after screening
# dom for, p95 > 5, cc > 50%
dat_screen <-
  read.csv(
    'D:/ontario_inventory/imputation/dat_screen_domfor_p95_cc.csv'
  )

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID, ]

# change age values to 2018 value
dat$AGE2018 <- 2018 - dat$YRORG

# select columns we need
dat %<>% select(POLYID, WG, SPCOMP,
                cc, avg, max, p95, 
                qav, ske, kur, cv, lor, 
                ba, qmdbh, dens, agb, 
                top_height, v, v_merch,
                B6, depth_q25, rumple, zentropy,
                zpcum8, zsd, x, y, AGE2018)

# remove any missing values
dat <- na.omit(dat)

####################
###EXPLORE SPCOMP###
####################

# look at first few rows
head(dat$SPCOMP)

# parse SPCOMP strings
sp <- str_split(dat$SPCOMP, pattern = "\\s{2}")

# add first species to dat
dat$SP1 <- sapply(sp, FUN = function(x){
  str <- x[1]
  str <- str_sub(str, start = 1, end = 2)
  return(str)
})

# add first species percent to dat
dat$SP1P <- sapply(sp, FUN = function(x){
  str <- x[2]
  if(is.na(str)){
    str <- 100
  } else{
    str <- str_sub(str, start = 1, end = 2)
  }
  return(str)
})

# add second species to dat
dat$SP2 <- sapply(sp, FUN = function(x){
  str <- x[2]
  if(is.na(str) == F){
    str <- str_sub(str, start = 3, end = 4)
  }
  return(str)
})

# add second species percent to dat
dat$SP2P <- sapply(sp, FUN = function(x){
  str <- x[3]
  if(is.na(str) == F){
    str <- str_sub(str, start = 1, end = 2)
  }
  return(str)
})

# load sp_group data
sp_group <- read.csv(
  'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest_SPGROUP.shp'
)

# change POLYID to numeric
dat %<>% mutate(POLYID = as.numeric(POLYID))
sp_group %<>% mutate(POLYID = as.numeric(POLYID))

# join to dat
dat <- left_join(dat, 
                 sp_group %>% select(POLYID, SpeciesGroup2, SpeciesGroup3), 
                 by = 'POLYID')

####################################################
###FUNCTIONS TO RUN K NEAREST NEIGHBOR IMPUTATION###
####################################################

# create mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# create rmse function
rmse <- function(obs, est){
  sqrt(mean((obs - est) ^ 2))
}

# create mae function
mae <- function(obs, est){
  sum(abs(obs-est)) / length(obs)
}

# create mape function
mape <- function(obs, est){
  sum(abs((obs - est) / obs)) / length(obs) * 100
}

# create knn function
run_knn <- function(dat, vars, k) {
  
  # subset data
  dat_nn <- dat %>% select(all_of(vars))
  
  # scale for nn computation
  dat_nn_scaled <- dat_nn %>% scale
  
  # run nearest neighbor
  nn <- nn2(dat_nn_scaled, dat_nn_scaled, k = k + 1)
  
  # get nn indices
  nni <- nn[[1]][, 2:(k + 1)]
  
  # add vars to tibble
  
  # take mean/mode if k > 1
  if(k > 1){
    for(i in seq_along(vars)){
      if(i == 1){
        nn_tab <- tibble(!!vars[i] := dat_nn[,i],
                         !!str_c(vars[i], '_nn') := apply(nni, MARGIN = 1, FUN = function(x){
                           mean(dat_nn[x, i])
                         }))
      }else{
        nn_tab %<>% mutate(!!vars[i] := dat_nn[,i],
                           !!str_c(vars[i], '_nn') := apply(nni, MARGIN = 1, FUN = function(x){
                             mean(dat_nn[x, i])
                           }))
      }
    }
    
    # add aux vars to tibble
    nn_tab %<>% mutate(age = dat$AGE,
                       wg = dat$WG,
                       sp1 = dat$SP1,
                       sp2 = dat$SP2,
                       group5 = dat$SpeciesGroup2,
                       group3 = dat$SpeciesGroup3,
                       age_nn = apply(nni, MARGIN = 1, FUN = function(x){
                         mean(dat$AGE[x])
                       }),
                       wg_nn = apply(nni, MARGIN = 1, FUN = function(x){
                         getmode(dat$WG[x])
                       }),
                       sp1_nn = apply(nni, MARGIN = 1, FUN = function(x){
                         getmode(dat$SP1[x])
                       }),
                       sp2_nn = apply(nni, MARGIN = 1, FUN = function(x){
                         getmode(dat$SP2[x])
                       }),
                       group5_nn = apply(nni, MARGIN = 1, FUN = function(x){
                         getmode(dat$SpeciesGroup2[x])
                       }),
                       group3_nn = apply(nni, MARGIN = 1, FUN = function(x){
                         getmode(dat$SpeciesGroup3[x])
                       }))
  }
  
  # take direct nn if k == 1
  if(k == 1){
    for(i in seq_along(vars)){
      if(i == 1){
        nn_tab <- tibble(!!vars[i] := dat_nn[,i],
                         !!str_c(vars[i], '_nn') := dat_nn[nn[[1]][,2],i])
      }else{
        nn_tab %<>% mutate(!!vars[i] := dat_nn[,i],
                           !!str_c(vars[i], '_nn') := dat_nn[nn[[1]][,2],i])
      }
    }
    
    # add aux vars to tibble
    nn_tab %<>% mutate(age = dat$AGE2018,
                       wg = dat$WG,
                       sp1 = dat$SP1,
                       sp2 = dat$SP2,
                       group5 = dat$SpeciesGroup2,
                       group3 = dat$SpeciesGroup3,
                       age_nn = dat$AGE2018[nn[[1]][,2]],
                       wg_nn = dat$WG[nn[[1]][,2]],
                       sp1_nn = dat$SP1[nn[[1]][,2]],
                       sp2_nn = dat$SP2[nn[[1]][,2]],
                       group5_nn = dat$SpeciesGroup2[nn[[1]][,2]],
                       group3_nn = dat$SpeciesGroup3[nn[[1]][,2]])
  }
  
  
  # calculate rmse for vars
  for(i in seq_along(vars)){
    if(i == 1){
      perform_df <- tibble(!!str_c(vars[i], '_raw_rmsd') := 
                             rmse(pull(nn_tab, vars[i]),
                                  pull(nn_tab, str_c(vars[i], '_nn'))),
                           !!str_c(vars[i], '_scaled_rmsd') := 
                             rmse(pull(nn_tab, vars[i]),
                                  pull(nn_tab, str_c(vars[i], '_nn'))) / 
                             sd(pull(nn_tab, vars[i])),
                           !!str_c(vars[i], '_mae') := 
                             mae(pull(nn_tab, vars[i]),
                                 pull(nn_tab, str_c(vars[i], '_nn'))),
                           !!str_c(vars[i], '_mape') := 
                             mape(pull(nn_tab, vars[i]),
                                  pull(nn_tab, str_c(vars[i], '_nn'))))
    }else{
      perform_df %<>% mutate(!!str_c(vars[i], '_raw_rmsd') := 
                               rmse(pull(nn_tab, vars[i]),
                                    pull(nn_tab, str_c(vars[i], '_nn'))),
                             !!str_c(vars[i], '_scaled_rmsd') := 
                               rmse(pull(nn_tab, vars[i]),
                                    pull(nn_tab, str_c(vars[i], '_nn'))) / 
                               sd(pull(nn_tab, vars[i])),
                             !!str_c(vars[i], '_mae') := 
                               mae(pull(nn_tab, vars[i]),
                                   pull(nn_tab, str_c(vars[i], '_nn'))),
                             !!str_c(vars[i], '_mape') := 
                               mape(pull(nn_tab, vars[i]),
                                    pull(nn_tab, str_c(vars[i], '_nn'))))
    }
  }
  
  # calculate rmse for aux vars
  perform_df %<>% mutate(age_raw_rmse = rmse(nn_tab$age, nn_tab$age_nn),
                         age_scaled_rmse = rmse(nn_tab$age, nn_tab$age_nn) /
                           sd(nn_tab$age),
                         age_mae = mae(nn_tab$age, nn_tab$age_nn),
                         age_mape = mape(nn_tab$age, nn_tab$age_nn))
  
  # calculate wg accuracy
  # create df of WG
  wg <- data.frame(obs = nn_tab$wg,
                   est = nn_tab$wg_nn)
  
  # create column of match or not
  wg$match <- wg$obs == wg$est
  
  # add total percent of matching WG to perform_df
  perform_df <- cbind(perform_df,
                      data.frame(wg_accuracy = NROW(wg[wg$match == T,]) /
                                   NROW(wg)))
  
  # calculate SP1 accuracy
  # create df of SP1
  sp1 <- data.frame(obs = nn_tab$sp1,
                    est = nn_tab$sp1_nn)
  
  # create column of match or not
  sp1$match <- sp1$obs == sp1$est
  
  # add total percent of matching SP1 to perform_df
  perform_df <- cbind(perform_df,
                      data.frame(sp1_accuracy = NROW(sp1[sp1$match == T,]) /
                                   NROW(sp1)))
  
  # calculate SP2 accuracy
  # create df of SP2
  sp2 <- data.frame(obs = nn_tab$sp2,
                    est = nn_tab$sp2_nn)
  
  # create column of match or not
  sp2$match <- sp2$obs == sp2$est
  
  # add total percent of matching SP2 to perform_df
  perform_df <- cbind(perform_df,
                      data.frame(sp2_accuracy = NROW(sp2[sp2$match == T,]) /
                                   NROW(sp2)))
  
  # calculate GROUP3 accuracy
  # create df of GROUP3
  group3 <- data.frame(obs = nn_tab$group3,
                       est = nn_tab$group3_nn)
  
  # create column of match or not
  group3$match <- group3$obs == group3$est
  
  # add total percent of matching SP2 to perform_df
  perform_df <- cbind(perform_df,
                      data.frame(group3_accuracy = NROW(group3[group3$match == T,]) /
                                   NROW(group3)))
  
  # calculate GROUP5 accuracy
  # create df of GROUP5
  group5 <- data.frame(obs = nn_tab$group5,
                       est = nn_tab$group5_nn)
  
  # create column of match or not
  group5$match <- group5$obs == group5$est
  
  # add total percent of matching SP2 to perform_df
  perform_df <- cbind(perform_df,
                      data.frame(group5_accuracy = NROW(group5[group5$match == T,]) /
                                   NROW(group5)))
  
  # melt df
  perform_df <- melt(perform_df)
  
  # return df
  return(perform_df)
}

##############################################
### RUN KNN BEST COMBINATIONS TO CALC MAPE ###
##############################################

# best als metrics only

# select variables
vars <- c('B6', 'cc',
          'rumple',
          'x',
          'y',
          'zpcum8',
          'zsd')

# run knn
df <- run_knn(dat, vars, k = 5)

# filter specific columns
df %<>% filter(str_detect(variable, 'mape')) %>%
  mutate(value = round(value, 2))

# best efi attributes only

# select variables
vars <- c('B6', 'dens',
          'depth_q25',
          'lor',
          'top_height',
          'x',
          'y')

# run knn
df <- run_knn(dat, vars, k = 5)

# filter specific columns
df %<>% filter(str_detect(variable, 'mape')) %>%
  mutate(value = round(value, 2))

      
