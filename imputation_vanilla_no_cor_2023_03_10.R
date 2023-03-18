# This code runs an imputation algorithm on LiDAR attributes connecting them to interpreter
# derived response attributes resulting in a model that can be applied to forest stands
# segmented from LiDAR rasters

# This version of the code uses a simple nearest neighbor approach to query
# the attributes of choice and find the sample that minimizes euclidean
# distance between the variables

# we can then impute all the attributes of the nearest neighbor into the 
# newly generated polygon

# Updated with new performance metrics:
# MD and RMSD (Absolute and relative)

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

# create rmsd function
rmsd <- function(obs, est){
  sqrt(mean((est - obs) ^ 2))
}

# create rrmsd function
rrmsd <- function(obs, est){
  sqrt(mean((est - obs) ^ 2)) / mean(obs) * 100
}

# create md function
md <- function(obs, est){
  mean(est - obs)
}

# create rmd function
rmd <- function(obs, est){
  mean(est - obs) / mean(obs) * 100
}

# create mae function
mae <- function(obs, est){
  mean(abs(est - obs))
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
  
  
  # calculate performance metrics for vars
  for(i in seq_along(vars)){
    if(i == 1){
      perform_df <- tibble(!!str_c(vars[i], '_rmsd') := 
                             rmsd(pull(nn_tab, vars[i]),
                                  pull(nn_tab, str_c(vars[i], '_nn'))),
                           !!str_c(vars[i], '_rrmsd') := 
                             rrmsd(pull(nn_tab, vars[i]),
                                   pull(nn_tab, str_c(vars[i], '_nn'))),
                           !!str_c(vars[i], '_md') := 
                             md(pull(nn_tab, vars[i]),
                                pull(nn_tab, str_c(vars[i], '_nn'))),
                           !!str_c(vars[i], '_rmd') := 
                             rmd(pull(nn_tab, vars[i]),
                                 pull(nn_tab, str_c(vars[i], '_nn'))),
                           !!str_c(vars[i], '_mae') := 
                             mae(pull(nn_tab, vars[i]),
                                 pull(nn_tab, str_c(vars[i], '_nn'))))
    }else{
      perform_df %<>% mutate(!!str_c(vars[i], '_rmsd') := 
                               rmsd(pull(nn_tab, vars[i]),
                                    pull(nn_tab, str_c(vars[i], '_nn'))),
                             !!str_c(vars[i], '_rrmsd') := 
                               rrmsd(pull(nn_tab, vars[i]),
                                     pull(nn_tab, str_c(vars[i], '_nn'))),
                             !!str_c(vars[i], '_md') := 
                               md(pull(nn_tab, vars[i]),
                                  pull(nn_tab, str_c(vars[i], '_nn'))),
                             !!str_c(vars[i], '_rmd') := 
                               rmd(pull(nn_tab, vars[i]),
                                   pull(nn_tab, str_c(vars[i], '_nn'))),
                             !!str_c(vars[i], '_mae') := 
                               mae(pull(nn_tab, vars[i]),
                                   pull(nn_tab, str_c(vars[i], '_nn'))))
    }
  }
  
  # calculate perf metrics for aux vars
  perform_df %<>% mutate(age_rmsd = rmsd(nn_tab$age, nn_tab$age_nn),
                         age_rrmsd = rrmsd(nn_tab$age, nn_tab$age_nn),
                         age_md = md(nn_tab$age, nn_tab$age_nn),
                         age_rmd = rmd(nn_tab$age, nn_tab$age_nn),
                         age_mae = mae(nn_tab$age, nn_tab$age_nn))
  
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

################################################################
### RUN KNN IMPUTATION OF MODELLED VARIABLES + B6 (group 1b) ###
################################################################

# select variables
vars <- c('lor', 'ba', 'qmdbh', 'dens',
          'agb', 'top_height', 'v', 
          'v_merch', 
          'x', 'y', 'B6')

# check correlation
corrplot(cor(dat %>% select(all_of(vars))), method = 'number')

# allocate list of combinations
comb <- list()

for(i in c(3, 5, 7)){
  
  # generate combinations
  c <- combinations(n = length(vars), r = i, v = vars)
  
  # morph into list
  d <- apply(c, 1, function(x){
    list(x)
  })
  
  # fix list nesting
  d <- lapply(d, function(x){
    x[[1]]
  })
  
  # append comb
  comb <- c(comb, d)
}

# split comb into 10 chunks
comb <- split(comb,
              ceiling(seq_along(comb) / (length(comb)/10)))

# register parallel backend
cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)

# loop through chunks
foreach::foreach(i = seq_along(comb)) %dopar% {
  
  # load packages
  library(terra)
  library(tidyverse)
  library(magrittr)
  library(RANN)
  library(corrplot)
  library(reshape2)
  library(gtools)
  library(janitor)
  
  # get list of variables to run inside this chunk
  var_list <- comb[[i]]
  
  # loop over all possible combinations
  for(x in seq_along(var_list)){
    
    # loop over different values of k
    for(k in c(1, 3, 5)){
      # run knn
      df <- suppressMessages(run_knn(dat, var_list[[x]], k = k))
      
      # filter specific columns
      df %<>% filter(variable %in% c('age_rmsd', 'age_rrmsd',
                                     'age_md', 'age_rmd',
                                     'age_mae', 'sp1_accuracy',
                                     'sp2_accuracy', 'group3_accuracy',
                                     'group5_accuracy'))
      
      # add number of vars and k
      df %<>% add_row(variable = 'n', value = length(var_list[[x]])) %>%
        add_row(variable = 'k', value = k)
      
      # transpose
      df %<>% t %>% as_tibble
      
      # set col names
      colnames(df) <- c('age_rmsd', 'age_rrmsd',
                        'age_md', 'age_rmd',
                        'age_mae', 'sp1_accuracy',
                        'sp2_accuracy', 'group3_accuracy',
                        'group5_accuracy',
                        'n', 'k')
      
      # drop variable row
      df %<>% filter(!row_number() == 1)
      
      # add vars attribute
      df %<>% add_column(vars = paste(var_list[[x]], collapse = " "))
      
      # save output
      if(x == 1 & k == 1){
        out <- df
      }else {
        out <- rbind(out, df)
      }
    }
  }
  
  # write out and read back in
  write.csv(out, str_c('D:/ontario_inventory/imputation/vanilla/group1b_mod_vars_chunk_', i, '2023_02_28.csv'))
}

# stop parallel cluster
parallel::stopCluster(cl)

#################################################
### RUN KNN IMPUTATION OF ALS + RAW VARIABLES ###
#################################################

# select variables
vars <- c('cc', 'avg', 'max',
          'p95', 'qav', 'ske',
          'kur', 'cv', 'B6', 'rumple',
          'zentropy', 'zpcum8',
          'depth_q25',
          'zsd', 'x', 'y')

# check correlation
corrplot(cor(dat %>% select(all_of(vars))), method = 'number')

# allocate list of combinations
comb <- list()

for(i in c(3, 5, 7)){
  
  # generate combinations
  c <- combinations(n = length(vars), r = i, v = vars)
  
  # morph into list
  d <- apply(c, 1, function(x){
    list(x)
  })
  
  # fix list nesting
  d <- lapply(d, function(x){
    x[[1]]
  })
  
  # append comb
  comb <- c(comb, d)
}

# split comb into 10 chunks
comb <- split(comb,
              ceiling(seq_along(comb) / (length(comb)/10)))

# register parallel backend
cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)

# loop through chunks
foreach::foreach(i = seq_along(comb)) %dopar% {
  
  # load packages
  library(terra)
  library(tidyverse)
  library(magrittr)
  library(RANN)
  library(corrplot)
  library(reshape2)
  library(gtools)
  library(janitor)
  
  # get list of variables to run inside this chunk
  var_list <- comb[[i]]
  
  # loop over all possible combinations
  for(x in seq_along(var_list)){
    
    # loop over different values of k
    for(k in c(1, 3, 5)){
      # run knn
      df <- suppressMessages(run_knn(dat, var_list[[x]], k = k))
      
      # filter specific columns
      df %<>% filter(variable %in% c('age_rmsd', 'age_rrmsd',
                                     'age_md', 'age_rmd',
                                     'age_mae', 'sp1_accuracy',
                                     'sp2_accuracy', 'group3_accuracy',
                                     'group5_accuracy'))
      
      # add number of vars and k
      df %<>% add_row(variable = 'n', value = length(var_list[[x]])) %>%
        add_row(variable = 'k', value = k)
      
      # transpose
      df %<>% t %>% as_tibble
      
      # set col names
      colnames(df) <- c('age_rmsd', 'age_rrmsd',
                        'age_md', 'age_rmd',
                        'age_mae', 'sp1_accuracy',
                        'sp2_accuracy', 'group3_accuracy',
                        'group5_accuracy',
                        'n', 'k')
      
      # drop variable row
      df %<>% filter(!row_number() == 1)
      
      # add vars attribute
      df %<>% add_column(vars = paste(var_list[[x]], collapse = " "))
      
      # save output
      if(x == 1 & k == 1){
        out <- df
      }else {
        out <- rbind(out, df)
      }
    }
  }
  
  # write out and read back in
  write.csv(out, str_c('D:/ontario_inventory/imputation/vanilla/group2_raw_vars_chunk_', i, '2023_02_28.csv'))
}

# stop parallel cluster
parallel::stopCluster(cl)

#############################################
### RUN KNN IMPUTATION ON JOINT VARIABLES ###
#############################################

# # select variables
# vars <- c('cc', 'avg', 'max',
#           'p95', 'qav', 'B6', 'rumple', 'zpcum8',
#           'zsd', 'x', 'y',
#           'depth_q25', 'agb',
#           'ba', 'lor', 'v', 'qmdbh',
#           'dens', 'top_height')

vars <- c('cc', 'avg', 'max',
          'p95', 'qav', 'ske',
          'kur', 'cv', 'B6', 'rumple',
          'zentropy', 'zpcum8',
          'zsd', 'x', 'y', 
          'lor', 'ba', 'qmdbh', 'dens',
          'agb', 'top_height', 'v', 
          'v_merch', 'depth_q25')

# allocate list of combinations
comb <- list()

for(i in c(3, 5, 7)){
  
  # generate combinations
  c <- combinations(n = length(vars), r = i, v = vars)
  
  # morph into list
  d <- apply(c, 1, function(x){
    list(x)
  })
  
  # fix list nesting
  d <- lapply(d, function(x){
    x[[1]]
  })
  
  # append comb
  comb <- c(comb, d)
}

# split comb into 20 chunks
comb <- split(comb,
              ceiling(seq_along(comb) / (length(comb)/20)))

# register parallel backend
cl <- parallel::makeCluster(20)
doParallel::registerDoParallel(cl)

# loop through chunks
foreach::foreach(i = seq_along(comb)) %dopar% {
  
  # load packages
  library(terra)
  library(tidyverse)
  library(magrittr)
  library(RANN)
  library(corrplot)
  library(reshape2)
  library(gtools)
  library(janitor)
  
  # get list of variables to run inside this chunk
  var_list <- comb[[i]]
  
  # loop over all possible combinations
  for(x in seq_along(var_list)){
    
    # loop over different values of k
    for(k in c(1, 3, 5)){
      # run knn
      df <- suppressMessages(run_knn(dat, var_list[[x]], k = k))
      
      # filter specific columns
      df %<>% filter(variable %in% c('age_rmsd', 'age_rrmsd',
                                     'age_md', 'age_rmd',
                                     'age_mae', 'sp1_accuracy',
                                     'sp2_accuracy', 'group3_accuracy',
                                     'group5_accuracy'))
      
      # add number of vars and k
      df %<>% add_row(variable = 'n', value = length(var_list[[x]])) %>%
        add_row(variable = 'k', value = k)
      
      # transpose
      df %<>% t %>% as_tibble
      
      # set col names
      colnames(df) <- c('age_rmsd', 'age_rrmsd',
                        'age_md', 'age_rmd',
                        'age_mae', 'sp1_accuracy',
                        'sp2_accuracy', 'group3_accuracy',
                        'group5_accuracy',
                        'n', 'k')
      
      # drop variable row
      df %<>% filter(!row_number() == 1)
      
      # add vars attribute
      df %<>% add_column(vars = paste(var_list[[x]], collapse = " "))
      
      # save output
      if(x == 1 & k == 1){
        out <- df
      }else {
        out <- rbind(out, df)
      }
    }
  }
  
  # write out and read back in
  write.csv(out, str_c('D:/ontario_inventory/imputation/vanilla/group3_comb_vars_chunk_', i, '2023_02_28.csv'))
}

# stop parallel cluster
parallel::stopCluster(cl)

#################################
### ORGANIZE OUTPUTS GROUP 1b ###
#################################

# load filenames to read back in
gr1_files <- list.files('D:/ontario_inventory/imputation/vanilla',
                        pattern = glob2rx('group1b_mod_vars_chunk_*2023_02_28.csv'),
                        full.names = T)

# read back in
gr1 <- do.call(rbind,lapply(gr1_files,read.csv))

# remove row names
gr1 %<>% select(-X)

# create function to rescale values from 0 to 100
scale_100 <- function(x){
  x <- (x-min(x))/(max(x) - min(x))
  return(x)
}

# lets scale age rmsd
gr1$age_rmsd_scaled <- scale(gr1$age_rmsd)

# take combined average of species comp acc
gr1$sp_acc <- rowMeans(gr1 %>% select(sp1_accuracy, sp2_accuracy, group3_accuracy, group5_accuracy))

# take error of combined species accuracy
gr1$sp_error <- 1 - gr1$sp_acc

# scale combined species error
gr1$sp_error_scaled <- scale(gr1$sp_error)

# take combined error mean age and species
gr1$error_mean <- rowMeans(gr1 %>% select(age_rmsd_scaled, sp_error_scaled))

#############################
### CALCULATE BEST MODELS ###
#############################

# best model
gr1_top <- gr1 %>% 
  filter(error_mean %in% sort(error_mean)[1]) %>% 
  select(vars, n, k, error_mean, age_rmsd, 
         age_mae, age_md, group3_accuracy, 
         group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
  arrange(error_mean) %>%
  add_column(result = 'top', .before = 1)

# best model k = 1
gr1_top %<>% add_row(gr1 %>% 
                       filter(k == 1) %>%
                       filter(error_mean %in% sort(error_mean)[1]) %>% 
                       select(vars, n, k, error_mean, age_rmsd, 
                              age_mae, age_md, group3_accuracy, 
                              group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                       arrange(error_mean) %>%
                       add_column(result = 'k = 1', .before = 1))

# best model n = 5
gr1_top %<>% add_row(gr1 %>% 
                       filter(n == 5) %>%
                       filter(error_mean %in% sort(error_mean)[1]) %>% 
                       select(vars, n, k, error_mean, age_rmsd, 
                              age_mae, age_md, group3_accuracy, 
                              group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                       arrange(error_mean) %>%
                       add_column(result = 'n = 5', .before = 1))

# best model n = 3
gr1_top %<>% add_row(gr1 %>% 
                       filter(n == 3) %>%
                       filter(error_mean %in% sort(error_mean)[1]) %>% 
                       select(vars, n, k, error_mean, age_rmsd, 
                              age_mae, age_md, group3_accuracy, 
                              group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                       arrange(error_mean) %>%
                       add_column(result = 'n = 3', .before = 1))

# best model k = 1, n = 5
gr1_top %<>% add_row(gr1 %>% 
                       filter(k == 1) %>%
                       filter(n == 5) %>%
                       filter(error_mean %in% sort(error_mean)[1]) %>% 
                       select(vars, n, k, error_mean, age_rmsd, 
                              age_mae, age_md, group3_accuracy, 
                              group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                       arrange(error_mean) %>%
                       add_column(result = 'k = 1 n = 5', .before = 1))

# best model k = 3, n = 7
gr1_top %<>% add_row(gr1 %>% 
                       filter(k == 3) %>%
                       filter(n == 7) %>%
                       filter(error_mean %in% sort(error_mean)[1]) %>% 
                       select(vars, n, k, error_mean, age_rmsd, 
                              age_mae, age_md, group3_accuracy, 
                              group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                       arrange(error_mean) %>%
                       add_column(result = 'k = 3 n = 7', .before = 1))

# best model k = 3, n = 5
gr1_top %<>% add_row(gr1 %>% 
                       filter(k == 3) %>%
                       filter(n == 5) %>%
                       filter(error_mean %in% sort(error_mean)[1]) %>% 
                       select(vars, n, k, error_mean, age_rmsd, 
                              age_mae, age_md, group3_accuracy, 
                              group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                       arrange(error_mean) %>%
                       add_column(result = 'k = 3 n = 5', .before = 1))

# write out
write.csv(gr1_top, 'D:/ontario_inventory/imputation/vanilla/gr1b_top_2023_02_28.csv')

##############################
### CALCULATE TOP 5 MODELS ###
##############################

# best models
gr1_top_5 <- gr1 %>% 
  filter(error_mean %in% sort(error_mean)[1:5]) %>% 
  select(vars, n, k, error_mean, age_rmsd, 
         age_mae, age_md, group3_accuracy, 
         group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
  arrange(error_mean) %>%
  add_column(result = 'top 5', .before = 1)

# best models k = 1
gr1_top_5 %<>% add_row(gr1 %>% 
                         filter(k == 1) %>%
                         filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                         select(vars, n, k, error_mean, age_rmsd, 
                                age_mae, age_md, group3_accuracy, 
                                group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                         arrange(error_mean) %>%
                         add_column(result = 'top 5 k = 1', .before = 1))

# best models n = 5
gr1_top_5 %<>% add_row(gr1 %>% 
                         filter(n == 5) %>%
                         filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                         select(vars, n, k, error_mean, age_rmsd, 
                                age_mae, age_md, group3_accuracy, 
                                group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                         arrange(error_mean) %>%
                         add_column(result = 'top 5 n = 5', .before = 1))

# best models n = 3
gr1_top_5 %<>% add_row(gr1 %>% 
                         filter(n == 3) %>%
                         filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                         select(vars, n, k, error_mean, age_rmsd, 
                                age_mae, age_md, group3_accuracy, 
                                group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                         arrange(error_mean) %>%
                         add_column(result = 'top 5 n = 3', .before = 1))

# best models k = 1, n = 5
gr1_top_5 %<>% add_row(gr1 %>% 
                         filter(k == 1) %>%
                         filter(n == 5) %>%
                         filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                         select(vars, n, k, error_mean, age_rmsd, 
                                age_mae, age_md, group3_accuracy, 
                                group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                         arrange(error_mean) %>%
                         add_column(result = 'top 5 k = 1 n = 5', .before = 1))

# best models k = 3, n = 7
gr1_top_5 %<>% add_row(gr1 %>% 
                         filter(k == 3) %>%
                         filter(n == 7) %>%
                         filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                         select(vars, n, k, error_mean, age_rmsd, 
                                age_mae, age_md, group3_accuracy, 
                                group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                         arrange(error_mean) %>%
                         add_column(result = 'top 5 k = 3 n = 5', .before = 1))

# best models k = 3, n = 5
gr1_top_5 %<>% add_row(gr1 %>% 
                         filter(k == 3) %>%
                         filter(n == 5) %>%
                         filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                         select(vars, n, k, error_mean, age_rmsd, 
                                age_mae, age_md, group3_accuracy, 
                                group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                         arrange(error_mean) %>%
                         add_column(result = 'top 5 k = 3 n = 5', .before = 1))

# write out
write.csv(gr1_top_5, 'D:/ontario_inventory/imputation/vanilla/gr1b_top_5_2023_02_28.csv')

# select variables
vars <- c('lor', 'ba', 'qmdbh', 'dens',
          'agb', 'top_height', 'v', 
          'v_merch', 'x', 'y')

# look at correlation of variables
cor <- cor(dat[, vars], use = 'complete.obs')
corrplot(cor, method = 'number')

################################
### ORGANIZE OUTPUTS GROUP 2 ###
################################

# load filenames to read back in
gr_files <- list.files('D:/ontario_inventory/imputation/vanilla',
                       pattern = glob2rx('group2_raw_vars_chunk_*2023_02_28.csv'),
                       full.names = T)

# read back in
gr <- do.call(rbind,lapply(gr_files,read.csv))

# remove row names
gr %<>% select(-X)

# create function to rescale values from 0 to 100
scale_100 <- function(x){
  x <- (x-min(x))/(max(x) - min(x))
  return(x)
}

# lets scale age rmsd
gr$age_rmsd_scaled <- scale(gr$age_rmsd)

# take combined average of species comp acc
gr$sp_acc <- rowMeans(gr %>% select(sp1_accuracy, sp2_accuracy, group3_accuracy, group5_accuracy))

# take error of combined species accuracy
gr$sp_error <- 1 - gr$sp_acc

# scale combined species error
gr$sp_error_scaled <- scale(gr$sp_error)

# take combined error mean age and species
gr$error_mean <- rowMeans(gr %>% select(age_rmsd_scaled, sp_error_scaled))

#############################
### CALCULATE BEST MODELS ###
#############################

# best model
gr_top <- gr %>% 
  filter(error_mean %in% sort(error_mean)[1]) %>% 
  select(vars, n, k, error_mean, age_rmsd, 
         age_mae, age_md, group3_accuracy, 
         group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
  arrange(error_mean) %>%
  add_column(result = 'top', .before = 1)

# best model k = 1
gr_top %<>% add_row(gr %>% 
                      filter(k == 1) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 1', .before = 1))

# best model n = 5
gr_top %<>% add_row(gr %>% 
                      filter(n == 5) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'n = 5', .before = 1))

# best model n = 3
gr_top %<>% add_row(gr %>% 
                      filter(n == 3) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'n = 3', .before = 1))

# best model k = 1, n = 5
gr_top %<>% add_row(gr %>% 
                      filter(k == 1) %>%
                      filter(n == 5) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 1 n = 5', .before = 1))

# best model k = 3, n = 7
gr_top %<>% add_row(gr %>% 
                      filter(k == 3) %>%
                      filter(n == 7) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 3 n = 7', .before = 1))

# best model k = 3, n = 5
gr_top %<>% add_row(gr %>% 
                      filter(k == 3) %>%
                      filter(n == 5) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 3 n = 5', .before = 1))

# write out
write.csv(gr_top, 'D:/ontario_inventory/imputation/vanilla/gr2_top_2023_02_28.csv')

##############################
### CALCULATE TOP 5 MODELS ###
##############################

# best models
gr_top_5 <- gr %>% 
  filter(error_mean %in% sort(error_mean)[1:5]) %>% 
  select(vars, n, k, error_mean, age_rmsd, 
         age_mae, age_md, group3_accuracy, 
         group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
  arrange(error_mean) %>%
  add_column(result = 'top 5', .before = 1)

# best models k = 1
gr_top_5 %<>% add_row(gr %>% 
                        filter(k == 1) %>%
                        filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                        select(vars, n, k, error_mean, age_rmsd, 
                               age_mae, age_md, group3_accuracy, 
                               group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                        arrange(error_mean) %>%
                        add_column(result = 'top 5 k = 1', .before = 1))

# best models n = 5
gr_top_5 %<>% add_row(gr %>% 
                        filter(n == 5) %>%
                        filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                        select(vars, n, k, error_mean, age_rmsd, 
                               age_mae, age_md, group3_accuracy, 
                               group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                        arrange(error_mean) %>%
                        add_column(result = 'top 5 n = 5', .before = 1))

# best models n = 3
gr_top_5 %<>% add_row(gr %>% 
                        filter(n == 3) %>%
                        filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                        select(vars, n, k, error_mean, age_rmsd, 
                               age_mae, age_md, group3_accuracy, 
                               group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                        arrange(error_mean) %>%
                        add_column(result = 'top 5 n = 3', .before = 1))

# best models k = 1, n = 5
gr_top_5 %<>% add_row(gr %>% 
                        filter(k == 1) %>%
                        filter(n == 5) %>%
                        filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                        select(vars, n, k, error_mean, age_rmsd, 
                               age_mae, age_md, group3_accuracy, 
                               group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                        arrange(error_mean) %>%
                        add_column(result = 'top 5 k = 1 n = 5', .before = 1))

# best models k = 3, n = 7
gr_top_5 %<>% add_row(gr %>% 
                        filter(k == 3) %>%
                        filter(n == 7) %>%
                        filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                        select(vars, n, k, error_mean, age_rmsd, 
                               age_mae, age_md, group3_accuracy, 
                               group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                        arrange(error_mean) %>%
                        add_column(result = 'top 5 k = 3 n = 5', .before = 1))

# best models k = 3, n = 5
gr_top_5 %<>% add_row(gr %>% 
                        filter(k == 3) %>%
                        filter(n == 5) %>%
                        filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                        select(vars, n, k, error_mean, age_rmsd, 
                               age_mae, age_md, group3_accuracy, 
                               group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                        arrange(error_mean) %>%
                        add_column(result = 'top 5 k = 3 n = 5', .before = 1))

# write out
write.csv(gr_top_5, 'D:/ontario_inventory/imputation/vanilla/gr2_top_5_2023_02_28.csv')

# select variables
vars <- c('cc', 'avg', 'max',
          'p95', 'qav', 'ske',
          'kur', 'cv', 'B6', 'rumple',
          'zentropy', 'zpcum8',
          'depth_q25',
          'zsd', 'x', 'y')

# look at correlation of variables
cor <- cor(dat[, vars], use = 'complete.obs')
corrplot(cor, method = 'number')

################################
### ORGANIZE OUTPUTS GROUP 3 ###
################################

# load filenames to read back in
gr_files <- list.files('D:/ontario_inventory/imputation/vanilla',
                       pattern = 'group3_comb_vars_chunk_*2023_02_28.csv',
                       full.names = T)

# read back in
gr <- do.call(rbind,lapply(gr_files,read.csv))

# remove row names
gr %<>% select(-X)

# lets scale age rmsd
gr$age_rmsd_scaled <- scale(gr$age_rmsd)

# take combined average of species comp acc
gr$sp_acc <- rowMeans(gr %>% select(sp1_accuracy, sp2_accuracy, group3_accuracy, group5_accuracy))

# take error of combined species accuracy
gr$sp_error <- 1 - gr$sp_acc

# scale combined species error
gr$sp_error_scaled <- scale(gr$sp_error)

# take combined error mean age and species
gr$error_mean <- rowMeans(gr %>% select(age_rmsd_scaled, sp_error_scaled))

#############################
### CALCULATE BEST MODELS ###
#############################

# best model
gr_top <- gr %>% 
  filter(error_mean %in% sort(error_mean)[1]) %>% 
  select(vars, n, k, error_mean, age_rmsd, 
         age_mae, age_md, group3_accuracy, 
         group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
  arrange(error_mean) %>%
  add_column(result = 'top', .before = 1)

# best model k = 1
gr_top %<>% add_row(gr %>% 
                      filter(k == 1) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 1', .before = 1))

# best model n = 5
gr_top %<>% add_row(gr %>% 
                      filter(n == 5) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'n = 5', .before = 1))

# best model n = 3
gr_top %<>% add_row(gr %>% 
                      filter(n == 3) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'n = 3', .before = 1))

# best model k = 1, n = 5
gr_top %<>% add_row(gr %>% 
                      filter(k == 1) %>%
                      filter(n == 5) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 1 n = 5', .before = 1))

# best model k = 3, n = 7
gr_top %<>% add_row(gr %>% 
                      filter(k == 3) %>%
                      filter(n == 7) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 3 n = 7', .before = 1))

# best model k = 3, n = 5
gr_top %<>% add_row(gr %>% 
                      filter(k == 3) %>%
                      filter(n == 5) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 3 n = 5', .before = 1))

# write out
write.csv(gr_top, 'D:/ontario_inventory/imputation/vanilla/gr3_top_2023_01_25.csv')

##############################
### CALCULATE TOP 5 MODELS ###
##############################

# best models
gr_top_5 <- gr %>% 
  filter(error_mean %in% sort(error_mean)[1:5]) %>% 
  select(vars, n, k, error_mean, age_rmsd, 
         age_mae, age_md, group3_accuracy, 
         group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
  arrange(error_mean) %>%
  add_column(result = 'top 5', .before = 1)

# best models k = 1
gr_top_5 %<>% add_row(gr %>% 
                        filter(k == 1) %>%
                        filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                        select(vars, n, k, error_mean, age_rmsd, 
                               age_mae, age_md, group3_accuracy, 
                               group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                        arrange(error_mean) %>%
                        add_column(result = 'top 5 k = 1', .before = 1))

# best models n = 5
gr_top_5 %<>% add_row(gr %>% 
                        filter(n == 5) %>%
                        filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                        select(vars, n, k, error_mean, age_rmsd, 
                               age_mae, age_md, group3_accuracy, 
                               group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                        arrange(error_mean) %>%
                        add_column(result = 'top 5 n = 5', .before = 1))

# best models n = 3
gr_top_5 %<>% add_row(gr %>% 
                        filter(n == 3) %>%
                        filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                        select(vars, n, k, error_mean, age_rmsd, 
                               age_mae, age_md, group3_accuracy, 
                               group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                        arrange(error_mean) %>%
                        add_column(result = 'top 5 n = 3', .before = 1))

# best models k = 1, n = 5
gr_top_5 %<>% add_row(gr %>% 
                        filter(k == 1) %>%
                        filter(n == 5) %>%
                        filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                        select(vars, n, k, error_mean, age_rmsd, 
                               age_mae, age_md, group3_accuracy, 
                               group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                        arrange(error_mean) %>%
                        add_column(result = 'top 5 k = 1 n = 5', .before = 1))

# best models k = 3, n = 7
gr_top_5 %<>% add_row(gr %>% 
                        filter(k == 3) %>%
                        filter(n == 7) %>%
                        filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                        select(vars, n, k, error_mean, age_rmsd, 
                               age_mae, age_md, group3_accuracy, 
                               group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                        arrange(error_mean) %>%
                        add_column(result = 'top 5 k = 3 n = 5', .before = 1))

# best models k = 3, n = 5
gr_top_5 %<>% add_row(gr %>% 
                        filter(k == 3) %>%
                        filter(n == 5) %>%
                        filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                        select(vars, n, k, error_mean, age_rmsd, 
                               age_mae, age_md, group3_accuracy, 
                               group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                        arrange(error_mean) %>%
                        add_column(result = 'top 5 k = 3 n = 5', .before = 1))

# write out
write.csv(gr_top_5, 'D:/ontario_inventory/imputation/vanilla/gr3_top_5_2023_02_28.csv')

# select variables
vars <- c('cc', 'avg', 'max',
          'p95', 'qav', 'ske',
          'kur', 'cv', 'B6', 'rumple',
          'zentropy', 'zpcum8',
          'zsd', 'x', 'y')

# look at correlation of variables
cor <- cor(dat[, vars], use = 'complete.obs')
corrplot(cor, method = 'number')

# ########################################################
# ###RUN NEAREST NEIGHBOR IMPUTATION USING FRI AS QUERY###
# ########################################################
# 
# # create vector of variables
# # vars <- c('cc', 'p95', 'qav', 'cv', 'dens', 'agb', 'top_height',
# #           'B6', 'B11', 'PZABOVEMEAN', 'PZABOVE2', 'depth_mean',
# #           'depth_q25', 'depth_q50', 'depth_q75',
# #           'zpcum8', 'zpcum9', 'zq85', 'x', 'y',
# #           'zsd', 'zskew', 'rumple')
# 
# # create smaller vector of variables
# vars <- c('cc', 'cv',
#           'B6', 'B11', 'PZABOVEMEAN',
#           'depth_q25',
#           'zpcum8', 'zq85', 'x', 'y',
#           'zsd', 'rumple')
# 
# # allocate list of combinations
# comb <- list()
# 
# # loop over different lengths of r
# for(i in 2:5){
#   
#   # generate combinations
#   c <- combinations(n = length(vars), r = i, v = vars)
#   
#   # morph into list
#   d <- apply(c, 1, function(x){
#     list(x)
#   })
#   
#   # fix list nesting
#   d <- lapply(d, function(x){
#     x[[1]]
#   })
#   
#   # append comb
#   comb <- c(comb, d)
# }
# 
# # loop over all possible combinations
# for(x in seq_along(comb)){
#   # run knn
#   df <- suppressMessages(run_knn(dat, comb[[x]], k = 10))
#   
#   # filter specific columns
#   df %<>% filter(variable %in% c('age_raw_rmse', 'age_scaled_rmse',
#                                  'wg_accuracy', 'sp1_accuracy',
#                                  'sp2_accuracy', 'group3_accuracy',
#                                  'group5_accuracy'))
#   
#   # save output
#   if(x == 1){
#     out <- data.frame(df$value)
#   }else {
#     out <- cbind(out, data.frame(df$value))
#   }
# }
# 
# # write out and read back in
# write.csv(out, 'D:/ontario_inventory/imputation/vanilla/out.csv')
# out <- read.csv('D:/ontario_inventory/imputation/vanilla/out.csv')
# 
# # remove row names
# out %<>% select(-X)
# 
# # transpose
# out %<>% t
# 
# # set colnames
# colnames(out) <- c('age_raw_rmse', 'age_scaled_rmse',
#                    'wg_accuracy', 'sp1_accuracy',
#                    'sp2_accuracy', 'group3_accuracy',
#                    'group5_accuracy')
# 
# # create vector of imputation attributes
# for(i in seq_along(comb)){
#   if(i == 1){
#     attr <- paste(comb[[i]], collapse = " ")
#   } else{
#     attr <- append(attr, paste(comb[[i]], collapse = " "))
#   }
# }
# 
# # add to out df
# out <- cbind(out, tibble(attr = attr))
# 
# # check results
# out$attr[out$age_raw_rmse == min(out$age_raw_rmse)]
# out$attr[out$group3_accuracy == max(out$group3_accuracy)]
# out$attr[out$group5_accuracy == max(out$group5_accuracy)]
# 
# # depth_q25, zpcum8, red_edge_2 (B6)
# df <- run_knn(dat, c('depth_q25', 'zpcum8', 'B6'), k = 10)
# 
# df <- run_knn(dat, c('depth_q25', 'zpcum8', 'B6', 'B11'), k = 10)
# 
# df <- run_knn(dat, c('depth_q25', 'zpcum8', 'zentropy',
#                      'zq85', 'B6', 'B11'), k = 10)
# 
# df <- run_knn(dat, c('top_height', 'cc', 'B6'), k = 10)
# 
# df <- run_knn(dat, c('B6', 'rumple', 'zsd',
#                      'zq85', 'cc', 'zpcum8'), k = 10)
# 
# df2 <- run_knn(dat, c('B6', 'rumple', 'zsd',
#                       'p95', 'cc', 'zpcum8',
#                       'depth_q25'), k = 10)
# 
# # best?
# df3 <- run_knn(dat, c('p95', 'rumple',
#                       'zpcum8', 'depth_q25',
#                       'x', 'y', 'B6'), k = 10)
# 
# df4 <- run_knn(dat, c('B6', 'rumple',
#                       'p95', 'zpcum8',
#                       'depth_q25', 'x', 'y',
#                       'cc', 'zsd'), k = 10)
# 
# cor(dat %>% select('B6', 'rumple',
#                    'p95', 'zpcum8',
#                    'depth_q25', 'x', 'y',
#                    'cc', 'zsd'))
# 
# #################################################
# ###top_height, cc, red_edge_2, SWIR_2 (K = 10)###
# #################################################
# 
# # red edge 2 is B6
# # SWIR_2 is B11
# 
# # create mode function
# getmode <- function(v) {
#   uniqv <- unique(v)
#   uniqv[which.max(tabulate(match(v, uniqv)))]
# }
# 
# # subset data
# dat_nn <- dat %>% select(top_height, cc, B6, B11)
# 
# # scale for nn computation
# dat_nn_scaled <- dat_nn %>% scale
# 
# # run nearest neighbor
# nn <- nn2(dat_nn_scaled, dat_nn_scaled, k = 11)
# 
# # get nn indices
# nni <- nn[[1]][, 2:11]
# 
# # create tibble of results
# nn_tab <- tibble(top_height = dat_nn[,1],
#                  cc = dat_nn[,2],
#                  red_edge_2 = dat_nn[,3],
#                  swir_2 = dat_nn[,4],
#                  age = dat$AGE,
#                  wg = dat$WG,
#                  sp1 = dat$SP1,
#                  sp2 = dat$SP2,
#                  top_height_nn = apply(nni, MARGIN = 1, FUN = function(x){
#                    mean(dat_nn[x, 1])
#                  }),
#                  cc_nn = apply(nni, MARGIN = 1, FUN = function(x){
#                    mean(dat_nn[x, 2])
#                  }),
#                  red_edge_2_nn = apply(nni, MARGIN = 1, FUN = function(x){
#                    mean(dat_nn[x, 3])
#                  }),
#                  swir_2_nn = apply(nni, MARGIN = 1, FUN = function(x){
#                    mean(dat_nn[x, 4])
#                  }),
#                  age_nn = apply(nni, MARGIN = 1, FUN = function(x){
#                    mean(dat$AGE[x])
#                  }),
#                  wg_nn = apply(nni, MARGIN = 1, FUN = function(x){
#                    getmode(dat$WG[x])
#                  }),
#                  sp1_nn = apply(nni, MARGIN = 1, FUN = function(x){
#                    getmode(dat$SP1[x])
#                  }),
#                  sp2_nn = apply(nni, MARGIN = 1, FUN = function(x){
#                    getmode(dat$SP2[x])
#                  }),
#                  group5 = dat$SpeciesGroup2,
#                  group3 = dat$SpeciesGroup3,
#                  group5_nn = apply(nni, MARGIN = 1, FUN = function(x){
#                    getmode(dat$SpeciesGroup2[x])
#                  }),
#                  group3_nn = apply(nni, MARGIN = 1, FUN = function(x){
#                    getmode(dat$SpeciesGroup3[x])
#                  }))
# 
# # calculate rmse different variables
# rmse <- function(obs, est){
#   sqrt(mean((obs - est) ^ 2))
# }
# 
# # calculate rmse metrics
# perform_df <-
#   data.frame(
#     top_height_raw_rmse = rmse(nn_tab$top_height, nn_tab$top_height_nn),
#     cc_raw_rmse = rmse(nn_tab$cc, nn_tab$cc_nn),
#     red_edge_2_rmse = rmse(nn_tab$red_edge_2, nn_tab$red_edge_2_nn),
#     swir_2_rmse = rmse(nn_tab$swir_2, nn_tab$swir_2_nn),
#     age_raw_rmse = rmse(nn_tab$age, nn_tab$age_nn),
#     top_height_scaled_rmse = rmse(nn_tab$top_height, nn_tab$top_height_nn) /
#       sd(nn_tab$top_height),
#     cc_scaled_rmse = rmse(nn_tab$cc, nn_tab$cc_nn) /
#       sd(nn_tab$cc),
#     red_edge_2_scaled_rmse = rmse(nn_tab$red_edge_2, nn_tab$red_edge_2_nn) / 
#       sd(nn_tab$red_edge_2),
#     swir_2_scaled_rmse = rmse(nn_tab$swir_2, nn_tab$swir_2_nn) / 
#       sd(nn_tab$swir_2),
#     age_scaled_rmse = rmse(nn_tab$age, nn_tab$age_nn) /
#       sd(nn_tab$age),
#     top_height_mae = mae(nn_tab$top_height, nn_tab$top_height_nn),
#     age_mae = mae(nn_tab$age, nn_tab$age_nn)
#   )
# 
# # calculate wg accuracy
# # create df of WG
# wg <- data.frame(obs = nn_tab$wg,
#                  est = nn_tab$wg_nn)
# 
# # create column of match or not
# wg$match <- wg$obs == wg$est
# 
# # add total percent of matching WG to perform_df
# perform_df <- cbind(perform_df,
#                     data.frame(wg_accuracy = NROW(wg[wg$match == T,]) /
#                                  NROW(wg)))
# 
# # calculate SP1 accuracy
# # create df of SP1
# sp1 <- data.frame(obs = nn_tab$sp1,
#                   est = nn_tab$sp1_nn)
# 
# # create column of match or not
# sp1$match <- sp1$obs == sp1$est
# 
# # add total percent of matching SP1 to perform_df
# perform_df <- cbind(perform_df,
#                     data.frame(sp1_accuracy = NROW(sp1[sp1$match == T,]) /
#                                  NROW(sp1)))
# 
# # calculate SP2 accuracy
# # create df of SP2
# sp2 <- data.frame(obs = nn_tab$sp2,
#                   est = nn_tab$sp2_nn)
# 
# # create column of match or not
# sp2$match <- sp2$obs == sp2$est
# 
# # add total percent of matching SP2 to perform_df
# perform_df <- cbind(perform_df,
#                     data.frame(sp2_accuracy = NROW(sp2[sp2$match == T,]) /
#                                  NROW(sp2)))
# 
# # calculate GROUP3 accuracy
# # create df of GROUP3
# group3 <- data.frame(obs = nn_tab$group3,
#                      est = nn_tab$group3_nn)
# 
# # create column of match or not
# group3$match <- group3$obs == group3$est
# 
# # add total percent of matching SP2 to perform_df
# perform_df <- cbind(perform_df,
#                     data.frame(group3_accuracy = NROW(group3[group3$match == T,]) /
#                                  NROW(group3)))
# 
# # calculate GROUP5 accuracy
# # create df of GROUP5
# group5 <- data.frame(obs = nn_tab$group5,
#                      est = nn_tab$group5_nn)
# 
# # create column of match or not
# group5$match <- group5$obs == group5$est
# 
# # add total percent of matching SP2 to perform_df
# perform_df <- cbind(perform_df,
#                     data.frame(group5_accuracy = NROW(group5[group5$match == T,]) /
#                                  NROW(group5)))
# 
# ############################
# ###top_height, ba, cc, cv###
# ############################
# 
# dat_nn <- dat %>% select(top_height, ba, cc, cv)
# 
# # scale for nn computation
# dat_nn_scaled <- dat_nn %>% scale
# 
# # run nearest neighbor
# nn <- nn2(dat_nn_scaled, dat_nn_scaled, k = 2)
# 
# # create tibble of results
# nn_tab <- tibble(top_height = dat_nn[,1],
#                  ba = dat_nn[,2],
#                  cc = dat_nn[,3],
#                  cv = dat_nn[,4],
#                  age = dat$AGE,
#                  wg = dat$WG,
#                  sp1 = dat$SP1,
#                  sp2 = dat$SP2,
#                  group5 = dat$SpeciesGroup2,
#                  group3 = dat$SpeciesGroup3,
#                  top_height_nn = dat_nn[nn[[1]][,2],1],
#                  ba_nn = dat_nn[nn[[1]][,2],2],
#                  cc_nn = dat_nn[nn[[1]][,2],3],
#                  cv_nn = dat_nn[nn[[1]][,2],4],
#                  age_nn = dat$AGE[nn[[1]][,2]],
#                  wg_nn = dat$WG[nn[[1]][,2]],
#                  sp1_nn = dat$SP1[nn[[1]][,2]],
#                  sp2_nn = dat$SP2[nn[[1]][,2]],
#                  group5_nn = dat$SpeciesGroup2[nn[[1]][,2]],
#                  group3_nn = dat$SpeciesGroup3[nn[[1]][,2]])
# 
# # calculate rmse different variables
# rmse <- function(obs, est){
#   sqrt(mean((obs - est) ^ 2))
# }
# 
# # calculate rmse metrics
# perform_df <-
#   data.frame(
#     top_height_raw_rmse = rmse(nn_tab$top_height, nn_tab$top_height_nn),
#     ba_raw_rmse = rmse(nn_tab$ba, nn_tab$ba_nn),
#     cc_raw_rmse = rmse(nn_tab$cc, nn_tab$cc_nn),
#     cv_raw_rmse = rmse(nn_tab$cv, nn_tab$cv_nn),
#     age_raw_rmse = rmse(nn_tab$age, nn_tab$age_nn),
#     top_height_scaled_rmse = rmse(nn_tab$top_height, nn_tab$top_height_nn) /
#       sd(nn_tab$top_height),
#     ba_scaled_rmse = rmse(nn_tab$ba, nn_tab$ba_nn) /
#       sd(nn_tab$ba),
#     cc_scaled_rmse = rmse(nn_tab$cc, nn_tab$cc_nn) /
#       sd(nn_tab$cc),
#     cv_scaled_rmse = rmse(nn_tab$cv, nn_tab$cv_nn) /
#       sd(nn_tab$cv),
#     age_scaled_rmse = rmse(nn_tab$age, nn_tab$age_nn) /
#       sd(nn_tab$age)
#   )
# 
# 
# # calculate wg accuracy
# # create df of WG
# wg <- data.frame(obs = nn_tab$wg,
#                  est = nn_tab$wg_nn)
# 
# # create column of match or not
# wg$match <- wg$obs == wg$est
# 
# # add total percent of matching WG to perform_df
# perform_df <- cbind(perform_df,
#                     data.frame(wg_accuracy = NROW(wg[wg$match == T,]) /
#                                  NROW(wg)))
# 
# # calculate SP1 accuracy
# # create df of SP1
# sp1 <- data.frame(obs = nn_tab$sp1,
#                   est = nn_tab$sp1_nn)
# 
# # create column of match or not
# sp1$match <- sp1$obs == sp1$est
# 
# # add total percent of matching SP1 to perform_df
# perform_df <- cbind(perform_df,
#                     data.frame(sp1_accuracy = NROW(sp1[sp1$match == T,]) /
#                                  NROW(sp1)))
# 
# # calculate SP2 accuracy
# # create df of SP2
# sp2 <- data.frame(obs = nn_tab$sp2,
#                   est = nn_tab$sp2_nn)
# 
# # create column of match or not
# sp2$match <- sp2$obs == sp2$est
# 
# # add total percent of matching SP2 to perform_df
# perform_df <- cbind(perform_df,
#                     data.frame(sp2_accuracy = NROW(sp2[sp2$match == T,]) /
#                                  NROW(sp2)))
# 
# # calculate GROUP3 accuracy
# # create df of GROUP3
# group3 <- data.frame(obs = nn_tab$group3,
#                      est = nn_tab$group3_nn)
# 
# # create column of match or not
# group3$match <- group3$obs == group3$est
# 
# # add total percent of matching SP2 to perform_df
# perform_df <- cbind(perform_df,
#                     data.frame(group3_accuracy = NROW(group3[group3$match == T,]) /
#                                  NROW(group3)))
# 
# # calculate GROUP5 accuracy
# # create df of GROUP5
# group5 <- data.frame(obs = nn_tab$group5,
#                      est = nn_tab$group5_nn)
# 
# # create column of match or not
# group5$match <- group5$obs == group5$est
# 
# # add total percent of matching SP2 to perform_df
# perform_df <- cbind(perform_df,
#                     data.frame(group5_accuracy = NROW(group5[group5$match == T,]) /
#                                  NROW(group5)))
# 
# ####################################
# ###RUN IMPUTATION OVER GRM OUTPUT###
# ####################################
# 
# # load output GRM polygon dataset
# load('D:/ontario_inventory/segmentation/grm/lidar_extracted/grm_10_01_05_ext_100.RData')
# 
# # subset forested polygons
# dat_lidar_for <- dat_lidar %>% filter(dom_for == 'Yes')
# 
# ##################################################
# ###RUN CORRELATION ANALYSIS ON LIDAR ATTRIBUTES###
# ##################################################
# 
# # both within FRI polygons and within LiDAR segmented polygons
# 
# # FRI polygons
# # cols 112:127 in dat
# cor <- cor(dat[, 112:127], use = 'complete.obs')
# corrplot(cor, method = 'number')
# 
# # GRM polygons
# cor <- cor(dat_lidar_for[,9:24], use = 'complete.obs')
# corrplot(cor, method = 'number')
# 
# # all polygons together
# cor_mat <- rbind(dat[,112:127], dat_lidar_for[, 9:24])
# cor <- cor(cor_mat, use = 'complete.obs')
# corrplot(cor, method = 'number')
# 
# #####################
# ###P95, CC, and CV###
# #####################
# 
# # create df for spl and fri metrics
# dat_spl <- dat_lidar_for %>% select(p95, cc, cv) %>% na.omit
# dat_fri <- dat %>% select(p95, cc, cv) %>% na.omit
# 
# # scale for nn computation
# # need to combine and scale all values together then separate again
# dat_splfri_scaled <- rbind(dat_spl, dat_fri) %>% scale
# dat_spl_scaled <- dat_splfri_scaled[1:NROW(dat_spl),]
# dat_fri_scaled <- dat_splfri_scaled[(NROW(dat_spl)+1):(NROW(dat_spl)+NROW(dat_fri)),]
# 
# # run nearest neighbor
# nn <- nn2(dat_fri_scaled, dat_spl_scaled, k = 1)
# 
# # create tibble of results
# nn_tab <- tibble(p95 = dat_spl$p95,
#                  cc = dat_spl$cc,
#                  cv = dat_spl$cv,
#                  p95_nn = dat_fri$p95[nn[[1]]],
#                  cc_nn = dat_fri$cc[nn[[1]]],
#                  cv_nn = dat_fri$cv[nn[[1]]])
# 
# # calculate rmse metrics
# perform_df <-
#   data.frame(
#     p95_raw_rmse = rmse(nn_tab$p95, nn_tab$p95_nn),
#     cc_raw_rmse = rmse(nn_tab$cc, nn_tab$cc_nn),
#     cv_raw_rmse = rmse(nn_tab$cv, nn_tab$cv_nn),
#     p95_scaled_rmse = rmse(nn_tab$p95, nn_tab$p95_nn) /
#       sd(nn_tab$p95),
#     cc_scaled_rmse = rmse(nn_tab$cc, nn_tab$cc_nn) /
#       sd(nn_tab$cc),
#     cv_scaled_rmse = rmse(nn_tab$cv, nn_tab$cv_nn) /
#       sd(nn_tab$cv)
#   )
# 
# # write.csv(perform_df, 'D:/ontario_inventory/imputation/vanilla/grm_als_als_performance.csv')
# 
# ###########################
# ###ALS HT, CC, BA vs FRI###
# ###########################
# 
# # create df for als and fri metrics
# dat_spl <- dat_lidar_for %>% select(lor, cc, ba) %>% na.omit
# dat_fri <- dat %>% select(HT, CC, BA) %>% na.omit
# 
# # change als colnames
# colnames(dat_spl) <- c('HT', 'CC', 'BA')
# 
# # scale for nn computation
# # need to combine and scale all values together then separate again
# dat_splfri_scaled <- rbind(dat_spl, dat_fri) %>% scale
# dat_spl_scaled <- dat_splfri_scaled[1:NROW(dat_spl),]
# dat_fri_scaled <- dat_splfri_scaled[(NROW(dat_spl)+1):(NROW(dat_spl)+NROW(dat_fri)),]
# 
# # run nearest neighbor
# nn <- nn2(dat_fri_scaled, dat_spl_scaled, k = 1)
# 
# # create tibble of results
# nn_tab <- tibble(ht = dat_spl$HT,
#                  cc = dat_spl$CC,
#                  ba = dat_spl$BA,
#                  ht_nn = dat_fri$HT[nn[[1]]],
#                  cc_nn = dat_fri$CC[nn[[1]]],
#                  ba_nn = dat_fri$BA[nn[[1]]])
# 
# # calculate rmse metrics
# perform_df <-
#   data.frame(
#     ht_raw_rmse = rmse(nn_tab$ht, nn_tab$ht_nn),
#     cc_raw_rmse = rmse(nn_tab$cc, nn_tab$cc_nn),
#     ba_raw_rmse = rmse(nn_tab$ba, nn_tab$ba_nn),
#     ht_scaled_rmse = rmse(nn_tab$ht, nn_tab$ht_nn) /
#       sd(nn_tab$ht),
#     cc_scaled_rmse = rmse(nn_tab$cc, nn_tab$cc_nn) /
#       sd(nn_tab$cc),
#     ba_scaled_rmse = rmse(nn_tab$ba, nn_tab$ba_nn) /
#       sd(nn_tab$ba)
#   )
# 
# #  write.csv(perform_df, 'D:/ontario_inventory/imputation/vanilla/grm_als_fri_performance.csv')
# 
# ############################
# ###top_height, ba, cc, cv###
# ############################
# 
# # create df for spl and fri metrics
# dat_spl <- dat_lidar_for %>% select(top_height, ba, cc, cv) %>% na.omit
# dat_fri <- dat %>% select(top_height, ba, cc, cv) %>% na.omit
# 
# # scale for nn computation
# # need to combine and scale all values together then separate again
# dat_splfri_scaled <- rbind(dat_spl, dat_fri) %>% scale
# dat_spl_scaled <- dat_splfri_scaled[1:NROW(dat_spl),]
# dat_fri_scaled <- dat_splfri_scaled[(NROW(dat_spl)+1):(NROW(dat_spl)+NROW(dat_fri)),]
# 
# # run nearest neighbor
# nn <- nn2(dat_fri_scaled, dat_spl_scaled, k = 1)
# 
# # create tibble of results
# nn_tab <- tibble(top_height = dat_spl$top_height,
#                  ba = dat_spl$ba,
#                  cc = dat_spl$cc,
#                  cv = dat_spl$cv,
#                  top_height_nn = dat_fri$top_height[nn[[1]]],
#                  ba_nn = dat_fri$ba[nn[[1]]],
#                  cc_nn = dat_fri$cc[nn[[1]]],
#                  cv_nn = dat_fri$cv[nn[[1]]])
# 
# # calculate rmse metrics
# perform_df <-
#   data.frame(
#     top_height_raw_rmse = rmse(nn_tab$top_height, nn_tab$top_height_nn),
#     ba_raw_rmse = rmse(nn_tab$ba, nn_tab$ba_nn),
#     cc_raw_rmse = rmse(nn_tab$cc, nn_tab$cc_nn),
#     cv_raw_rmse = rmse(nn_tab$cv, nn_tab$cv_nn),
#     top_height_scaled_rmse = rmse(nn_tab$top_height, nn_tab$top_height_nn) /
#       sd(nn_tab$top_height),
#     ba_scaled_rmse = rmse(nn_tab$ba, nn_tab$ba_nn) /
#       sd(nn_tab$ba),
#     cc_scaled_rmse = rmse(nn_tab$cc, nn_tab$cc_nn) /
#       sd(nn_tab$cc),
#     cv_scaled_rmse = rmse(nn_tab$cv, nn_tab$cv_nn) /
#       sd(nn_tab$cv)
#   )
