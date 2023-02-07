# Imputation {-#imputation}

## Intro {-#intro}

This page provides an example of the imputation approach used to estimate age and species composition in newly-generated forest stands. The imputation is based on a k-nearest neighbor algorithm. "Link" variables between the photo-interpreted FRI and Generic Region Merging (GRM) generated forest stands are identified and for each GRM forest polygon, the k-nearest neighbor FRI polygons that minimize the Euclidean distance of the link variables are identified and are used to impute age and species composition for the new polygon. If k = 1 the age and species variables are imputed directly from the best matching polygon. If k > 1 age is imputed as the mean age of k-matching polygons and species variables are imputed as the mode of k-matching polygons. The basic workflow is as follows:

0. Run GRM Segmentation to derive new polygon dataset (see previous posts)
1. Extract LiDAR attributes in FRI polygons and GRM segmented polygons
2. Screen FRI polygons to curate an optimal dataset to use for imputation
3. Run imputation on FRI polygons ONLY to assess imputation performance
4. Run imputation on GRM polygon dataset to derive estimates of age and species composition

## 1. Extract LiDAR attributes in FRI polygons and newly derived polygons {-#imputationp1}

Before running imputation we need to extract LiDAR attributes in all FRI polygons and GRM segmented polygons. For this example, we only extract the attributes we will use in the imputation (listed below) as well as the attributes needed for FRI polygon data screening in section 2 of this walk-through. The value extracted for each polygon is the median cell value, weighted by the fraction of each cell that is covered by the polygon.

The variables used in the imputation algorithm are as follows, and were selected with guidance from previous works as well as a sensitivity analysis:

* avg: Average height of returns (> 1.3 m classified as vegetation ??)
* zsd: Standard deviation of returns height (> 1.3 m classified as vegetation ??)
* rumple: Ratio of canopy outer surface area to ground surface area
* zpcum8: Cumulative percentage of LiDAR returns found in 80% percentile of LiDAR height
* x and y: Coordinates of polygon centroid
* red_edge_2: Sentinel 2 surface reflectance band (740 nm)

* p95: 95th percentile of LiDAR height returns > 1.3 m classified as vegetation
* depth_q25: 25th percentile of signal attenuation depth of all returns

## 1a. Extract attributes in FRI polygons {-#imputationp1a}




```r
# load packages
library(terra)
library(tidyverse)
library(exactextractr)
library(sf)
library(magrittr)
library(gridExtra)

# load FRI polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

# convert to df
dat_fri <- as.data.frame(poly)

# cbind centroids to dat
dat_fri <- cbind(dat_fri, centroids(poly) %>% crds)

# load LiDAR datasets we need to extract over polygons
# create named vector with variable names and data links

lidar <- c('p95' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_p95.tif',
           'avg' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_avg.tif',
           'zsd' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_std.tif',
           'rumple' = 'D:/ontario_inventory/romeo/SPL metrics/Z_METRICS_MOSAIC/individual/RMF_RUMPLE_MOSAIC_r_rumple.tif',
           'zpcum8' = 'D:/ontario_inventory/romeo/SPL metrics/Z_METRICS_MOSAIC/individual/RMF_Z_METRICS_MOSAIC_zpcum8.tif',
           #'depth_q25' = 'D:/ontario_inventory/romeo/SPL metrics/Z_METRICS_MOSAIC/individual/RMF_SAD_METRICS_MOSAIC_depth_q25.tif',
           'red_edge_2' = 'D:/ontario_inventory/romeo/Sentinel/red_edge_2.tif',
           'cc' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif')

# loop through LiDAR attributes to extract values
for (i in seq_along(lidar)) {
  # load LiDAR rasters as raster stack
  lidar_ras <- rast(lidar[i])
  
  # project poly to crs of raster
  poly_ras <- project(poly, lidar_ras)
  
  # convert to sf
  poly_ras <- st_as_sf(poly_ras)
  
  #extract median values
  vec <-
    exact_extract(lidar_ras, poly_ras, 'median')
  
  # aggregate into data frame
  if(i == 1){
    vec_df <- as.data.frame(vec)
  } else{
    vec_df <- cbind(vec_df, as.data.frame(vec))
  }
  
}

# change column names of extracted attribute data frame
colnames(vec_df) <- names(lidar)

# add LiDAR attributes to FRI polygon data frame
dat_fri <- cbind(dat_fri, vec_df)

# add 2018 age values
dat_fri$AGE2018 <- 2018 - dat_fri$YRORG

# save extracted dataframe for fast rebooting
save(dat_fri, file = 'D:/ontario_inventory/imputation/example/dat_fri_extr.RData')
```

## 1b. Extract attributes in newly segmented polygons {-#imputationp1b}


```r
# load GRM segmented polygons
poly <- vect('D:/ontario_inventory/segmentation/grm/shp/grm_10_01_05.shp')

# convert to df
dat_grm <- as.data.frame(poly)

# cbind centroids to dat
dat_grm <- cbind(dat_grm, centroids(poly) %>% crds)

# remove p95 and cc as we'll re calculate here for consistency
dat_grm %<>% select(-c(p95, cc))

# load LiDAR datasets we need to extract over polygons
# create named vector with variable names and data links

lidar <- c('p95' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_p95.tif',
           'avg' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_avg.tif',
           'zsd' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_std.tif',
           'rumple' = 'D:/ontario_inventory/romeo/SPL metrics/Z_METRICS_MOSAIC/individual/RMF_RUMPLE_MOSAIC_r_rumple.tif',
           'zpcum8' = 'D:/ontario_inventory/romeo/SPL metrics/Z_METRICS_MOSAIC/individual/RMF_Z_METRICS_MOSAIC_zpcum8.tif',
           #'depth_q25' = 'D:/ontario_inventory/romeo/SPL metrics/Z_METRICS_MOSAIC/individual/RMF_SAD_METRICS_MOSAIC_depth_q25.tif',
           'red_edge_2' = 'D:/ontario_inventory/romeo/Sentinel/red_edge_2.tif',
           'cc' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif')

# loop through LiDAR attributes to extract values
for (i in seq_along(lidar)) {
  # load LiDAR rasters as raster stack
  lidar_ras <- rast(lidar[i])
  
  # project poly to crs of raster
  poly_ras <- project(poly, lidar_ras)
  
  # convert to sf
  poly_ras <- st_as_sf(poly_ras)
  
  #extract median values
  vec <-
    exact_extract(lidar_ras, poly_ras, 'median')
  
  # aggregate into data frame
  if(i == 1){
    vec_df <- as.data.frame(vec)
  } else{
    vec_df <- cbind(vec_df, as.data.frame(vec))
  }
  
}

# change column names of extracted attribute data frame
colnames(vec_df) <- names(lidar)

# add LiDAR attributes to FRI polygon data frame
dat_grm <- cbind(dat_grm, vec_df)

# save extracted data frame for fast rebooting
save(dat_grm, file = 'D:/ontario_inventory/imputation/example/grm_10_01_05_extr.RData')

# clear workspace
rm(list=ls())
```

## 2. Screen FRI polygons to curate an optimal dataset to use for imputation {-#imputationp2}

In order to ensure the best possible imputation results, it is important to screen the FRI dataset and remove polygons that do not fit certain data quality criteria, or are not representative of forest stands.

The current criteria being used are:
a) POLYTYPE == 'FOR'
b) polygon >= 50% forested landcover
c) p95 (95th percentile of LiDAR height returns) >= 5 meters (definition of 'forest')
d) Canopy cover (% of returns > 2 m) >= 50%


```r
# load FRI polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

# load FRI polygon data frame
load('D:/ontario_inventory/imputation/example/dat_fri_extr.RData')

##################################
### SCREEN FOR POLYTYPE == FOR ###
##################################

# filter POLYTYPE == 'FOR'
dat_fri <- filter(dat_fri, POLYTYPE == 'FOR')
poly_fri <- poly[poly$POLYTYPE == 'FOR']

####################################################
### SCREEN FOR POLYGON >= 50% FORESTED LANDCOVER ###
####################################################

# load VLCE 2.0 landcover dataset from 2018
lc <- rast('D:/ontario_inventory/VLCE/CA_forest_VLCE2_2018_CLIPPED.tif')

# project poly to crs of raster
poly_lc <- project(poly_fri, lc)

# convert to sf
poly_lcsf <- st_as_sf(poly_lc)

# extract landcover values
lc_poly <- exact_extract(lc, poly_lcsf)

# set landcover class key
lc_key <- c(`20` = 'Water',
            `31` = 'Snow/Ice',
            `32` = 'Rock/Rubble',
            `33` = 'Exposed/Barren Land',
            `40` = 'Bryoids',
            `50` = 'Shrubland',
            `80` = 'Wetland',
            `81` = 'Wetland-Treed',
            `100` = 'Herbs',
            `210` = 'Coniferous',
            `220` = 'Broadleaf',
            `230` = 'Mixed Wood')

# # find number of unique lc types in each polygon
# # apply over list
# lc_uni <- sapply(lc_poly, function(x){
#   x$value <- recode(x$value, !!!lc_key)
#   x %<>% filter(coverage_fraction >= 0.5)
#   return(length(unique(x$value)))
# })

# set landcover class key with single forested class
lc_key_for <- c(`20` = 'Water',
                `31` = 'Snow/Ice',
                `32` = 'Rock/Rubble',
                `33` = 'Exposed/Barren Land',
                `40` = 'Bryoids',
                `50` = 'Shrubland',
                `80` = 'Wetland',
                `81` = 'Forest',
                `100` = 'Herbs',
                `210` = 'Forest',
                `220` = 'Forest',
                `230` = 'Forest')

# find pixels with forest at least 50% of pixel
# apply over list
lc_dom_for <- sapply(lc_poly, function(x){
  x$value <- recode(x$value, !!!lc_key_for)
  x <- x %>% group_by(value) %>% summarize(sum = sum(coverage_fraction))
  m <- x$value[which(x$sum == max(x$sum))]
  if((length(m) == 1) & (m == 'Forest')[1]){
    if(x$sum[x$value == m]/sum(x$sum) >= 0.5){
      return('Yes')
    }else{return('No')}
  }else{return('No')}
})

# add into FRI data frame
dat_fri <- dat_fri %>% add_column(dom_for = lc_dom_for)

# subset FRI data frame based on whether polygon dominated by forest
dat_dom_for <- dat_fri[dat_fri$dom_for == 'Yes',]

##################################
### SCREEN FOR P95 >= 5 METERS ###
##################################

# subset FRI data frame
dat_p95 <- dat_fri[dat_fri$p95 >= 5,]

############################
### SCREEN FOR CC >= 50% ###
############################

# subset FRI data frame
dat_cc <- dat_fri[dat_fri$cc >= 50,]

############################################
### COMBINE INTERSECTION OF DATA SCREENS ###
############################################

# first combine only by POLYID
dat_fri_scr <- intersect(dat_dom_for %>% subset(select = POLYID),
                        dat_p95 %>% subset(select = POLYID)) %>%
  intersect(., dat_cc %>% subset(select = POLYID))

# load full data frame attributes
dat_fri_scr <- dat_fri[dat_fri$POLYID %in% dat_fri_scr$POLYID,]

# save extracted data frame for fast rebooting
save(dat_fri_scr, file = 'D:/ontario_inventory/imputation/example/dat_fri_scr.RData')

# clear workspace
rm(list=ls())
```

## 3. Run imputation on FRI polygons ONLY to assess imputation performance {-#imputationp3}

The goal of this imputation procedure is to estimate age and species composition in newly segmented forest polygons. But it is also important to assess the performance of the algorithm. To do this, we can conduct the imputation over the FRI dataset ONLY. For each FRI polygon, we find the k-nearest neighbor matches, and calculate the age and species composition to impute. We can then compare the observed age and species composition of the polygon to the imputed values. Note we have to do this calculation on the FRI dataset alone because we do not have observed age and species composition values for the GRM segmented polygons. Also note that ALL attributes are imputed from the same k-nearest neighbors. The algorithm is not run for individual attributes.

We present the performance of numeric attributes (age as well as the variables used to conduct the imputation algorithm) as Mean Absolute Error (MAE) and Mean Absolute Percentage Error (MAPE). MAE gives an indication of performance in the units of each attribute (such as Age) while MAPE gives a percentage error that can be compared across variables.

For species composition, we report the percent of observed and imputed values that match (Accuracy). Species composition is broken down into several distinct attributes:

a) Working group
b) First dominant species (from FRI SPCOMP attribute)
c) Second dominant species (from FRI SPCOMP attribute)
d) 3-Group species classification (softwood, mixedwood, hardwood)
e) 5-Group species classification (jack pine dominated, black spruce dominated, mixed conifer, mixedwood, hardwood)

## 3a. Create functions needed {-#imputationp3a}


```r
####################################################
###FUNCTIONS TO RUN K NEAREST NEIGHBOR IMPUTATION###
####################################################

# load packages
library(RANN)
library(reshape2)

# create mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# create mae function
mae <- function(obs, est){
  sum(abs(obs - est)) / length(obs)
}

# create mape function
mape <- function(obs, est){
  sum(abs((obs - est) / obs)) / length(obs) * 100
}

# create knn function
run_knn_fri <- function(dat, vars, k) {
  
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
    nn_tab %<>% mutate(age = dat$AGE2018,
                       wg = dat$WG,
                       sp1 = dat$SP1,
                       sp2 = dat$SP2,
                       group5 = dat$SpeciesGroup2,
                       group3 = dat$SpeciesGroup3,
                       age_nn = apply(nni, MARGIN = 1, FUN = function(x){
                         mean(dat$AGE2018[x])
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
  
  
  # calculate fit metrics for vars
  for(i in seq_along(vars)){
    if(i == 1){
      perform_df <- tibble(!!str_c(vars[i], '_mae') := 
                             mae(pull(nn_tab, vars[i]),
                                  pull(nn_tab, str_c(vars[i], '_nn'))),
                           !!str_c(vars[i], '_mape') := 
                             mape(pull(nn_tab, vars[i]),
                                  pull(nn_tab, str_c(vars[i], '_nn'))))
    }else{
      perform_df %<>% mutate(!!str_c(vars[i], '_mae') := 
                             mae(pull(nn_tab, vars[i]),
                                  pull(nn_tab, str_c(vars[i], '_nn'))),
                           !!str_c(vars[i], '_mape') := 
                             mape(pull(nn_tab, vars[i]),
                                  pull(nn_tab, str_c(vars[i], '_nn'))))
    }
  }
  
  # calculate rmse for aux vars
  perform_df %<>% mutate(age_mae = mae(nn_tab$age, nn_tab$age_nn),
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
```

## 3b. Run the imputation and assess performance {-#imputationp3b}


```r
###########################
### LOAD FRI DATA FRAME ###
###########################

# load from part 2 above
load('D:/ontario_inventory/imputation/example/dat_fri_scr.RData')

# subset only the attributes we need
dat_fri_scr %<>% select(POLYID, AGE2018, SPCOMP, WG, 
                        avg, zsd,
                        rumple, zpcum8,
                        x, y, red_edge_2)

# remove any polygons with missing values
dat_fri_scr <- na.omit(dat_fri_scr)

###########################################
### CALCULATE AND ADD SPCOMP ATTRIBUTES ###
###########################################

# parse SPCOMP strings
sp <- str_split(dat_fri_scr$SPCOMP, pattern = "\\s{2}")

# add first species to dat
dat_fri_scr$SP1 <- sapply(sp, FUN = function(x){
  str <- x[1]
  str <- str_sub(str, start = 1, end = 2)
  return(str)
})

# add first species percent to dat
dat_fri_scr$SP1P <- sapply(sp, FUN = function(x){
  str <- x[2]
  if(is.na(str)){
    str <- 100
  } else{
    str <- str_sub(str, start = 1, end = 2)
  }
  return(str)
})

# add second species to dat
dat_fri_scr$SP2 <- sapply(sp, FUN = function(x){
  str <- x[2]
  if(is.na(str) == F){
    str <- str_sub(str, start = 3, end = 4)
  }
  return(str)
})

# add second species percent to dat
dat_fri_scr$SP2P <- sapply(sp, FUN = function(x){
  str <- x[3]
  if(is.na(str) == F){
    str <- str_sub(str, start = 1, end = 2)
  }
  return(str)
})

# change second species missing values
dat_fri_scr$SP2[is.na(dat_fri_scr$SP2)] <- 'MIS'
dat_fri_scr$SP2P[is.na(dat_fri_scr$SP2P)] <- 0

# load species group data -- calculated in separate code (can provide details)
sp_group <- read.csv(
  'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest_SPGROUP.shp'
)

# change POLYID to numeric
dat_fri_scr %<>% mutate(POLYID = as.numeric(POLYID))
sp_group %<>% mutate(POLYID = as.numeric(POLYID))

# join to dat
dat_fri_scr <- left_join(dat_fri_scr, 
                 sp_group %>% select(POLYID, SpeciesGroup2, SpeciesGroup3), 
                 by = 'POLYID')

##############################################
### RUN NEAREST NEIGHBOR IMPUTATION K = 5 ###
##############################################

# run_knn function already created
# we use top_height, cc, red_edge_2 in the algorithm
perf <- run_knn_fri(dat_fri_scr, c('avg', 'zsd', 'rumple', 
                                   'zpcum8', 'x', 
                                   'y', 'red_edge_2'), k = 5)

# round values
perf %<>% mutate(value = round(value, 2))

# display results
knitr::kable(perf, caption = "Imputation Performance of FRI Forest Stand Polygons", label = NA)
```

EDIT!

*** Mean Absolute Percent Error (MAPE) of the imputation attributes (p95, rumple, zpcum8, depth_q25, x, y, red_edge_2) is below 5% for all attributes. These are low values, which demonstrate that the imputation algorithm is finding optimal matches within the database of available FRI polygons. The Mean Absolute Error of Age is 16.27, and the MAPE is 38.71 percent. Accuracy of first species classification (sp1) is 67%, and a much lower 36% for second species. Species classification with 3 and 5 classes have respective accuracies of 74% and 63%. ***

## 4. Run imputation on GRM polygon dataset to derive estimates of age and species composition {-#imputationp4}

The last step is to run the imputation between the screened FRI data frame and the GRM segmented polygons. For each GRM segmented polygon, the k-nearest neighbors in the FRI data are found and used to impute age and species composition. Note we do not conduct imputation on all the GRM segmented polygons, but only the polygons that have >= 50% forested landcover, p95 >= 5 meters, and canopy cover >= 50% (same criteria as FRI data screening conducted above in step 2).

Although we cannot assess performance of age and species composition when imputing into the GRM segmented polygons, we can still assess the fit of the variables used in the imputation: avg, zsd, rumple, zpcum8, x, y, and Sentinel red_edge_2.

We can also review maps and distributions comparing FRI age/species composition against the same attributes imputed into GRM segmented polygons.


```r
# load additional packages
library(viridis)
library(scales)
library(janitor)

# load GRM polygon data frame
load('D:/ontario_inventory/imputation/example/grm_10_01_05_extr.RData')

# screen for forested polygons
dat_grm_scr <- dat_grm %>% filter(dom_for == 'Yes',
                                  p95 >= 5,
                                  cc >= 50)

# create data frame for grm and fri metrics used in imputation
dat_grm_imp <- dat_grm_scr %>% select(id, avg, zsd, 
                                      rumple, zpcum8,
                                      x, y, red_edge_2) %>% na.omit
dat_fri_imp <- dat_fri_scr %>% select(avg, zsd, 
                                      rumple, zpcum8,
                                      x, y, red_edge_2) %>% na.omit

# need to combine and scale all values together then separate again
dat_comb_scaled <- rbind(dat_grm_imp %>% select(-id), 
                         dat_fri_imp) %>% scale
dat_grm_scaled <- dat_comb_scaled[1:NROW(dat_grm_imp),]
dat_fri_scaled <- dat_comb_scaled[(NROW(dat_grm_imp)+1):(NROW(dat_grm_imp)+NROW(dat_fri_imp)),]

# run nearest neighbor imputation k = 5
nn <- nn2(dat_fri_scaled, dat_grm_scaled, k = 5)

# get nn indices
nni <- nn[[1]]

# add imputed attributes into GRM imputation data frame
dat_grm_imp %<>% add_column(
  avg_imp = apply(
    nni,
    MARGIN = 1,
    FUN = function(x) {
      mean(dat_fri_imp[x, 'avg'])
    }
  ),
  zsd_imp = apply(
    nni,
    MARGIN = 1,
    FUN = function(x) {
      mean(dat_fri_imp[x, 'zsd'])
    }
  ),
  rumple_imp = apply(
    nni,
    MARGIN = 1,
    FUN = function(x) {
      mean(dat_fri_imp[x, 'rumple'])
    }
  ),
  zpcum8_imp = apply(
    nni,
    MARGIN = 1,
    FUN = function(x) {
      mean(dat_fri_imp[x, 'zpcum8'])
    }
  ),
  x_imp = apply(
    nni,
    MARGIN = 1,
    FUN = function(x) {
      mean(dat_fri_imp[x, 'x'])
    }
  ),
  y_imp = apply(
    nni,
    MARGIN = 1,
    FUN = function(x) {
      mean(dat_fri_imp[x, 'y'])
    }
  ),
  red_edge_2_imp = apply(
    nni,
    MARGIN = 1,
    FUN = function(x) {
      mean(dat_fri_imp[x, 'red_edge_2'])
    }
  ),
  age = apply(
    nni,
    MARGIN = 1,
    FUN = function(x) {
      mean(dat_fri_scr[x, 'AGE2018'])
    }
  ),
  wg = apply(
    nni,
    MARGIN = 1,
    FUN = function(x) {
      getmode(dat_fri_scr[x, 'WG'])
    }
  ),
  sp1 = apply(
    nni,
    MARGIN = 1,
    FUN = function(x) {
      getmode(dat_fri_scr[x, 'SP1'])
    }
  ),
  sp2 = apply(
    nni,
    MARGIN = 1,
    FUN = function(x) {
      getmode(dat_fri_scr[x, 'SP2'])
    }
  ),
  class3 = apply(
    nni,
    MARGIN = 1,
    FUN = function(x) {
      getmode(dat_fri_scr[x, 'SpeciesGroup3'])
    }
  ),
  class5 = apply(
    nni,
    MARGIN = 1,
    FUN = function(x) {
      getmode(dat_fri_scr[x, 'SpeciesGroup2'])
    }
  )
)

# add values back into main GRM data frame (missing values for polygons not
# included in the imputation)
dat_grm <- left_join(dat_grm, dat_grm_imp)

# calculate performance across imputation attributes
perf <- tibble(avg_mae = mae(dat_grm_imp$avg, dat_grm_imp$avg_imp),
               avg_mape = mape(dat_grm_imp$avg, dat_grm_imp$avg_imp),
               zsd_mae = mae(dat_grm_imp$zsd, dat_grm_imp$zsd_imp),
               zsd_mape = mape(dat_grm_imp$zsd, dat_grm_imp$zsd_imp),
               rumple_mae = mae(dat_grm_imp$rumple, dat_grm_imp$rumple_imp),
               rumple_mape = mape(dat_grm_imp$rumple, dat_grm_imp$rumple_imp),
               zpcum8_mae = mae(dat_grm_imp$zpcum8, dat_grm_imp$zpcum8_imp),
               zpcum8_mape = mape(dat_grm_imp$zpcum8, dat_grm_imp$zpcum8_imp),
               x_mae = mae(dat_grm_imp$x, dat_grm_imp$x_imp),
               x_mape = mape(dat_grm_imp$x, dat_grm_imp$x_imp),
               y_mae = mae(dat_grm_imp$y, dat_grm_imp$y_imp),
               y_mape = mape(dat_grm_imp$y, dat_grm_imp$y_imp),
               red_edge_2_mae = mae(dat_grm_imp$red_edge_2, dat_grm_imp$red_edge_2_imp),
               red_edge_2_mape = mape(dat_grm_imp$red_edge_2, dat_grm_imp$red_edge_2_imp)
               ) %>% melt

# round to two decimal places
perf %<>% mutate(value = round(value, 2))

# display results of imputation rmsd
knitr::kable(perf, caption = "Imputation Performance between FRI and GRM Forest Stand Polygons", label = NA)
```

## 4a. Figures of Forest Stand Age {-#imputationp4a}


```r
# load GRM polygons
poly_grm <- vect('D:/ontario_inventory/segmentation/grm/shp/grm_10_01_05.shp')

# add new data frame to polygons
values(poly_grm) <- dat_grm

# save grm polygon output
writeVector(poly_grm, 'D:/ontario_inventory/imputation/example/grm_10_01_05_imp.shp', overwrite = T)

# load FRI polygons
poly_fri <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

# set age == 0 to NA
poly_fri$AGE[poly_fri$AGE == 0] <- NA

# load species group data -- calculated in separate code (can provide details)
sp_group <- read.csv(
  'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest_SPGROUP.shp'
)

# load poly_Fri dataframe
dat_fri <- as.data.frame(poly_fri)

# change POLYID to numeric
dat_fri %<>% mutate(POLYID = as.numeric(POLYID))
sp_group %<>% mutate(POLYID = as.numeric(POLYID))

# rename species groups
sp_group %<>% rename(class3 = SpeciesGroup3, class5 = SpeciesGroup2)

# join to dat
dat_fri <- left_join(dat_fri, 
                 sp_group %>% select(POLYID, class3, class5), 
                 by = 'POLYID')

# update age to 2018
dat_fri$AGE2018 <- 2018 - dat_fri$YRORG

# re-input attributes into FRI polygons
values(poly_fri) <- dat_fri
rm(dat_fri)

# since we only have values of age and species comp for dom_fom == yes
# and p95 >= 5 m we should only compare the screened FRI polygons that 
# contain those same attributes
poly_fri$AGE2018[!(poly_fri$POLYID %in% dat_fri_scr$POLYID)] <- NA
poly_fri$class3[!(poly_fri$POLYID %in% dat_fri_scr$POLYID)] <- NA
poly_fri$class5[!(poly_fri$POLYID %in% dat_fri_scr$POLYID)] <- NA

# plot age
par(mfrow=c(2,1))

plot(poly_grm, 'age', col = viridis(14), type = 'interval', breaks = seq(0, 140, 10),
     plg = list(x = 'topright', cex = 2, title = 'Age'),
     main = '')
title(main = 'Imputed Age of GRM Segmented Forest Stands', cex.main = 3)

plot(poly_fri, 'AGE2018', col = viridis(14), type = 'interval', breaks = seq(0, 140, 10),
     plg = list(x = 'topright', cex = 2, title = 'Age'),
     main = '')
title(main = 'Age of FRI Forest Stands', cex.main = 3)
```

<img src="05-Imputation_files/figure-html/imputationp4a-1.png" width="1728" />

## 4b. Figures of Forest Stand 3-Species Classification {-#imputationp4b}


```r
# plot group of 3 species
par(mfrow=c(2,1))

plot(poly_grm, 'class3', col = viridis(3), type = 'classes',
     plg = list(x = 'topright', cex = 2, title = 'Species Class'),
     main = '')
title(main = 'Imputed Species Class (3) of GRM Segmented Forest Stands', cex.main = 3)

plot(poly_fri, 'class3', col = viridis(3), type = 'classes',
     plg = list(x = 'topright', cex = 2, title = 'Species Class'),
     main = '')
title(main = 'Species Class (3) FRI Forest Stands', cex.main = 3)
```

<img src="05-Imputation_files/figure-html/imputationp4b-1.png" width="1728" />

## 4c. Figures of Forest Stand 5-Species Classification {-#imputationp4c}


```r
# plot group of 5 species
par(mfrow=c(2,1))

plot(poly_grm, 'class5', col = viridis(5), type = 'classes',
     plg = list(x = 'topright', cex = 2, title = 'Species Class'),
     main = '')
title(main = 'Imputed Species Class (5) of GRM Segmented Forest Stands', cex.main = 3)

plot(poly_fri, 'class5', col = viridis(5), type = 'classes',
     plg = list(x = 'topright', cex = 2, title = 'Species Class'),
     main = '')
title(main = 'Species Class (5) FRI Forest Stands', cex.main = 3)
```

<img src="05-Imputation_files/figure-html/imputationp4c-1.png" width="1728" />

## 4d. Density Plots of Forest Stand Age {-#imputationp4d}


```r
# density plots of age
p1 <- ggplot(dat_grm, aes(x = age)) +
  geom_density(fill = 'grey') +
  geom_vline(aes(xintercept = median(age, na.rm = T)),
             linetype = "dashed",
             size = 0.6) +
  xlim(c(0,200)) +
  ylim(c(0, 0.02)) +
  theme_bw() +
  xlab('Age') +
  ylab('Density') +
  ggtitle('Imputed Age of GRM Segmented Forest Stands') +
  theme(text = element_text(size = 25),
        plot.title = element_text(size=30))

p2 <- ggplot(as.data.frame(poly_fri), aes(x = AGE2018)) +
  geom_density(fill = 'grey') +
  geom_vline(aes(xintercept = median(AGE, na.rm = T)),
             linetype = "dashed",
             size = 0.6) +
  xlim(c(0, 200)) +
  ylim(c(0, 0.02)) +
  theme_bw() +
  xlab('Age') +
  ylab('Density') +
  ggtitle('Age of FRI Forest Stands') +
  theme(text = element_text(size = 25),
        plot.title = element_text(size=30))

grid.arrange(p1, p2, ncol = 2)
```

<img src="05-Imputation_files/figure-html/imputationp4d-1.png" width="1728" />

## 4e. Distribution of 3-Species Classification {-#imputationp4e}


```r
# distribution of 3 species classes
# create data frame for GRM 3 classes
dat_grm_c3 <- dat_grm %>% 
  tabyl(class3) %>%  
  filter(is.na(class3) == F) %>% 
  arrange(desc(class3)) %>%
  mutate(prop = n / sum(.$n)*100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop) %>%
  mutate(lbl = round(prop))

# create data frame for FRI 3 classes
dat_fri_c3 <- poly_fri %>% as.data.frame %>%
  tabyl(class3) %>%  
  filter(is.na(class3) == F) %>% 
  arrange(desc(class3)) %>%
  mutate(prop = n / sum(.$n)*100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop) %>%
  mutate(lbl = round(prop))

# plot
p1 <- ggplot(dat_grm_c3, aes(x = "", y = prop, fill = class3)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(y = ypos, label = str_c(lbl, "%")), size = 15) +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=30)) +
  labs(fill = "3 Species Classification") +
  ggtitle("Imputed 3 Species Classification \nof GRM Segmented Forest Stands")

p2 <- ggplot(dat_fri_c3, aes(x = "", y = prop, fill = class3)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(y = ypos, label = str_c(lbl, "%")), size = 15) +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=30)) +
  labs(fill = "3 Species Classification") +
  ggtitle("3 Species Classification \nof FRI Forest Stands")

grid.arrange(p1, p2, ncol = 2)
```

<img src="05-Imputation_files/figure-html/imputationp4e-1.png" width="1728" />

## 4f. Distribution of 5-Species Classification {-#imputationp4f}


```r
# distribution of 5 species classes
# create data frame for GRM 5 classes
dat_grm_c5 <- dat_grm %>% 
  tabyl(class5) %>%  
  filter(is.na(class5) == F) %>% 
  arrange(desc(class5)) %>%
  mutate(prop = n / sum(.$n)*100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop) %>%
  mutate(lbl = round(prop))

# create data frame for FRI 3 classes
dat_fri_c5 <- poly_fri %>% as.data.frame %>%
  tabyl(class5) %>%  
  filter(is.na(class5) == F) %>% 
  arrange(desc(class5)) %>%
  mutate(prop = n / sum(.$n)*100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop) %>%
  mutate(lbl = round(prop))

# plot
p1 <- ggplot(dat_grm_c5, aes(x = "", y = prop, fill = class5)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(y = ypos, label = str_c(lbl, "%")), size = 15) +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=30)) +
  labs(fill = "5 Species Classification") +
  ggtitle("Imputed 5 Species Classification \nof GRM Segmented Forest Stands")

p2 <- ggplot(dat_fri_c5, aes(x = "", y = prop, fill = class5)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(y = ypos, label = str_c(lbl, "%")), size = 15) +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=30)) +
  labs(fill = "5 Species Classification") +
  ggtitle("5 Species Classification \nof FRI Forest Stands")

grid.arrange(p1, p2, ncol = 2)
```

<img src="05-Imputation_files/figure-html/imputationp4f-1.png" width="1728" />
