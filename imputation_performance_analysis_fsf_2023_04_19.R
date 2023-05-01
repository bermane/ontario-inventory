# This code runs an imputation algorithm on LiDAR attributes
# connecting them to interpreter derived response attributes resulting in a 
# model that can be applied to forest stands segmented from LiDAR rasters

# This version of the code runs all possible iterations of model X-variables
# to assess which variables are optimal for final imputation model

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
load('D:/ontario_inventory/dat/dat_fri_fsf.RData')

# load photo interpreted polygons
poly <-
  vect(
    'D:/ontario_inventory/FSF/FRI/FSF_opi_polygon_CSRS_NAD83_17.shp'
  )

# cbind centroids to dat
dat <- cbind(dat, centroids(poly) %>% crds)

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
# dom for, p95 > 5, cc > 50%
dat_screen <-
  read.csv(
    'D:/ontario_inventory/imputation/fsf/dat_screen_fri_fsf.csv'
  )

# subset dat based on screening
dat <- dat[dat$OPI_ID %in% dat_screen$OPI_ID,]

# change age values to 2018 value
dat$AGE2018 <- 2018 - dat$YRORG

# select columns we need
dat %<>% select(OPI_ID, OLEADSPC, SPCOMP,
                cc, avg, max, p95, 
                ske, kur, cv, lor, 
                ba, qmdbh, stems, agb, 
                top_height, gtv, gmvwl,
                b6, rumple, zentropy,
                zpcum8, sd, x, y, AGE2018)

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
  'D:/ontario_inventory/imputation/fsf/fsf_SPGROUP.csv'
)

# change POLYID to numeric
dat %<>% mutate(OPI_ID = as.numeric(OPI_ID))
sp_group %<>% mutate(OPI_ID = as.numeric(OPI_ID))

# join to dat
dat <- left_join(dat, 
                 sp_group %>% select(OPI_ID, SpeciesGroup2, SpeciesGroup3), 
                 by = 'OPI_ID')

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

###################################################
### RUN KNN IMPUTATION ON FULL SET OF VARIABLES ###
###################################################

# select variables

vars <- c('cc', 'avg', 'max',
          'p95', 'ske',
          'kur', 'cv', 'b6', 'rumple',
          'zentropy', 'zpcum8',
          'sd', 'x', 'y', 
          'lor', 'ba', 'qmdbh',
          'agb', 'top_height', 'gtv', 
          'gmvwl', 'stems')

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
  
  # write out
  write.csv(out, str_c('D:/ontario_inventory/imputation/fsf/imp_perf_comb_vars_chunk_', i, '.csv'))
}

# stop parallel cluster
parallel::stopCluster(cl)

#################################
### ORGANIZE COMBINED OUTPUTS ###
#################################

# load filenames to read back in
gr_files <- list.files('D:/ontario_inventory/imputation/fsf',
                       pattern = glob2rx('imp_perf_comb_vars_chunk_*.csv'),
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

# best model k = 3
gr_top %<>% add_row(gr %>% 
                      filter(k == 3) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 3', .before = 1))

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
write.csv(gr_top, 'D:/ontario_inventory/imputation/fsf/gr3_top_2023_05_01.csv')

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

# best models k = 3
gr_top_5 %<>% add_row(gr %>% 
                        filter(k == 3) %>%
                        filter(error_mean %in% sort(error_mean)[1:5]) %>% 
                        select(vars, n, k, error_mean, age_rmsd, 
                               age_mae, age_md, group3_accuracy, 
                               group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                        arrange(error_mean) %>%
                        add_column(result = 'top 5 k = 3', .before = 1))

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
write.csv(gr_top_5, 'D:/ontario_inventory/imputation/fsf/gr3_top_5_2023_05_01.csv')

########################################################
### CALCULATE BEST MODELS ALS ONLY FROM COMB RESULTS ###
########################################################

# filter rows for ALS vars (OR for n=3 B6, x and y only)
# but since none of the other top runs included only those 3
# i think we're OK to not keep that specific selection

# best way is to remove rows with EFI attributes
# can't remove v since single letter
efi_vars <- c('lor', 'ba', 'qmdbh', 'agb',
          'top_height', 'gtv', 'gmvwl', 'stems')

gr_als <- gr %>% filter(str_detect(vars, paste(efi_vars, collapse = "|")) == F)


# best model
gr_top <- gr_als %>% 
  filter(error_mean %in% sort(error_mean)[1]) %>% 
  select(vars, n, k, error_mean, age_rmsd, 
         age_mae, age_md, group3_accuracy, 
         group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
  arrange(error_mean) %>%
  add_column(result = 'top', .before = 1)

# best model k = 1
gr_top %<>% add_row(gr_als %>% 
                      filter(k == 1) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 1', .before = 1))

# best model n = 5
gr_top %<>% add_row(gr_als %>% 
                      filter(n == 5) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'n = 5', .before = 1))

# best model n = 3
gr_top %<>% add_row(gr_als %>% 
                      filter(n == 3) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'n = 3', .before = 1))

# best model k = 1, n = 5
gr_top %<>% add_row(gr_als %>% 
                      filter(k == 1) %>%
                      filter(n == 5) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 1 n = 5', .before = 1))

# best model k = 3
gr_top %<>% add_row(gr_als %>% 
                      filter(k == 3) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 3', .before = 1))

# best model k = 3, n = 5
gr_top %<>% add_row(gr_als %>% 
                      filter(k == 3) %>%
                      filter(n == 5) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 3 n = 5', .before = 1))

# write out
write.csv(gr_top, 'D:/ontario_inventory/imputation/fsf/gr3_top_als_2023_05_01.csv')

########################################################
### CALCULATE BEST MODELS EFI ONLY FROM COMB RESULTS ###
########################################################

# filter rows for EFI vars (OR for n=3 B6, x and y only)
# but since none of the other top runs included only those 3
# i think we're OK to not keep that specific selection

# best way is to remove rows with ALS attributes
als_vars <- c('cc', 'avg', 'max',
              'p95', 'ske', 'kur',
              'cv', 'rumple',
              'zentropy', 'zpcum8',
              'sd')

gr_efi <- gr %>% filter(str_detect(vars, paste(als_vars, collapse = "|")) == F)


# best model
gr_top <- gr_efi %>% 
  filter(error_mean %in% sort(error_mean)[1]) %>% 
  select(vars, n, k, error_mean, age_rmsd, 
         age_mae, age_md, group3_accuracy, 
         group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
  arrange(error_mean) %>%
  add_column(result = 'top', .before = 1)

# best model k = 1
gr_top %<>% add_row(gr_efi %>% 
                      filter(k == 1) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 1', .before = 1))

# best model n = 5
gr_top %<>% add_row(gr_efi %>% 
                      filter(n == 5) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'n = 5', .before = 1))

# best model n = 3
gr_top %<>% add_row(gr_efi %>% 
                      filter(n == 3) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'n = 3', .before = 1))

# best model k = 1, n = 5
gr_top %<>% add_row(gr_efi %>% 
                      filter(k == 1) %>%
                      filter(n == 5) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 1 n = 5', .before = 1))

# best model k = 3
gr_top %<>% add_row(gr_efi %>% 
                      filter(k == 3) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 3', .before = 1))

# best model k = 3, n = 5
gr_top %<>% add_row(gr_efi %>% 
                      filter(k == 3) %>%
                      filter(n == 5) %>%
                      filter(error_mean %in% sort(error_mean)[1]) %>% 
                      select(vars, n, k, error_mean, age_rmsd, 
                             age_mae, age_md, group3_accuracy, 
                             group5_accuracy, sp1_accuracy, sp2_accuracy) %>%
                      arrange(error_mean) %>%
                      add_column(result = 'k = 3 n = 5', .before = 1))

# write out
write.csv(gr_top, 'D:/ontario_inventory/imputation/fsf/gr3_top_efi_2023_05_01.csv')
