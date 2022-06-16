# this code extracts final polygon summary stats from the best segmentation outputs

# load packages
library(terra)
library(tidyverse)
library(janitor)
library(magrittr)
library(sf)
library(lwgeom)
library(exactextractr)

# create list of polygon files and names
file <- list.files(path = 'D:/ontario_inventory/segmentation/grm/shp', 
                   pattern = glob2rx('*.shp'),
                   full.names = T)

file_names <- list.files(path = 'D:/ontario_inventory/segmentation/grm/shp', 
                         pattern = glob2rx('*.shp'))

file_names <- gsub("\\..*","", file_names)

# create standard error function
se <- function(x) sd(x)/sqrt(length(x))

#######################################
### EXTRACT P95 AND LOR CORRELATION ###
#######################################

# load p95
ras <- rast('D:/ontario_inventory/segmentation/spl_stack_mask_wat_ucl_polytype.tif')
ras <- ras[[1]]

# load lor
lor <- rast('D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_lor_FOR.tif') %>%
  project(ras) %>% mask(ras)

# as numeric
ras %<>% as.vector
lor %<>% as.vector

# check correlation
cor(ras, lor, use = 'complete.obs')

###########################################
### EXTRACT FINAL POLYGON SUMMARY STATS ###
###########################################

# loop through data sets
for(i in 1:length(file)){
  
  # load file
  p <- vect(file[i])
  
  # convert to sf
  p_sf <- st_as_sf(p)
  
  # calculate perimeter
  p$perim <- st_perimeter(p_sf) %>% as.numeric
  
  # calculate area
  p$area <- st_area(p_sf) %>% as.numeric
  
  # calculate perimeter to area ratio
  p$p_to_a <- p$perim/p$area
  p$p_to_a <- round(p$p_to_a, 3)
  
  # subset non masked WAT and UCL polygons
  p2_sf <- p[is.na(p$POLYTYPE)] %>% st_as_sf
  p2 <- p[is.na(p$POLYTYPE)] %>% as.data.frame
  
  # subset reliable forested polygons
  # p3 <- p %>% filter(dom_for == 'Yes', is.na(POLYTYPE), p95 >= 5)
  
  # calculate msi
  p2$msi <- p2$perim/sqrt(pi * p2$area)
  
  # load original raster input file
  ras <- rast('D:/ontario_inventory/segmentation/spl_stack_mask_wat_ucl_polytype.tif')
  
  # rename bands
  names(ras) <- c('p95', 'cc', 'cv')
  
  # extract pixel values
  pvals <- exact_extract(ras, p2_sf)
  
  # calculate SSE
  sse <- sapply(pvals, FUN = function(x){
    p95_mean <- mean(x$p95, na.rm = T)
    cc_mean <- mean(x$cc, na.rm = T)
    cv_mean <- mean(x$cv, na.rm = T)
    
    return(c(sum((x$p95 - p95_mean)^2, na.rm = T),
             sum((x$cc - cc_mean)^2, na.rm = T),
             sum((x$cv - cv_mean)^2, na.rm = T)))
  })
  
  # transpose
  sse <- t(sse)
  
  # calculate final sums
  sse <- colSums(sse)
  
  # unlist values
  pvals2 <- do.call(rbind, pvals)
  
  # calculate global mean values
  p95_mean <- mean(pvals2$p95, na.rm = T)
  cc_mean <- mean(pvals2$cc, na.rm = T)
  cv_mean <- mean(pvals2$cv, na.rm = T)
  
  rm(pvals2)
  
  # calculate SST
  sst <- sapply(pvals, FUN = function(x){
    return(c(sum((x$p95 - p95_mean)^2, na.rm = T),
             sum((x$cc - cc_mean)^2, na.rm = T),
             sum((x$cv - cv_mean)^2, na.rm = T)))
  })
  
  # transpose
  sst <- t(sst)
  
  # calculate final sums
  sst <- colSums(sst)
  
  # calculate r2 values
  r2_p95 <- 1 - (sse[1]/sst[1]) %>% round(3)
  r2_cc <- 1 - (sse[2]/sst[2]) %>% round(3)
  r2_cv <- 1 - (sse[3]/sst[3]) %>% round(3)
  r2_all <- (sum(r2_p95, r2_cc, r2_cv) / 3) %>% round(3)
  
  # create dataframe with values wanted
  ms_df <- data.frame(file = file_names[i],
                      min_pix = (min(p2$nbPixels)),
                      max_pix = (max(p2$nbPixels)),
                      mean_pix = (mean(p2$nbPixels)),
                      med_pix = (median(p2$nbPixels)),
                      num_poly = NROW(p2),
                      mean_area = mean(p2$area),
                      se_area = se(p2$area),
                      sd_area = sd(p2$area),
                      mean_perim = mean(p2$perim),
                      se_perim = se(p2$perim),
                      sd_perim = sd(p2$perim),
                      mean_p_a = mean(p2$p_to_a),
                      se_p_a = se(p2$p_to_a),
                      sd_p_a = sd(p2$p_to_a),
                      mean_msi = mean(p2$msi),
                      se_msi = se(p2$msi),
                      sd_msi = sd(p2$msi),
                      r2_p95 <- r2_p95,
                      r2_cc <- r2_cc,
                      r2_cv <- r2_cv,
                      r2_all <- r2_all)
  
  # round numeric columns
  ms_df %<>% 
    mutate_at(c('min_pix',
                'max_pix', 
                'mean_pix', 
                'med_pix', 
                'mean_area', 
                'se_area', 
                'sd_area',
                'mean_perim',
                'se_perim',
                'sd_perim'), 
              function(x) round(x, 2)) %>%
    mutate_at(c('mean_p_a', 
                'se_p_a', 
                'sd_p_a', 
                'mean_msi', 
                'se_msi', 
                'sd_msi'), 
              function(x) round(x, 4))
  
  # merge datasets
  if(i == 1){
    df <- ms_df
  }else{
    df <- rbind(df, ms_df)
  }
}

#####################
### ADD FRI STATS ###
#####################

# load interpreter derived polygons to extract statistics
p <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

# convert to sf
p_sf <- st_as_sf(p)

# calculate perimeter
p$perim <- st_perimeter(p_sf) %>% as.numeric

# calculate area
p$area <- st_area(p_sf) %>% as.numeric

# calculate perimeter to area ratio
p$p_to_a <- p$perim/p$area
p$p_to_a <- round(p$p_to_a, 3)

# subset all non water/ucl polygons
p2_sf <- p[!(p$POLYTYPE %in% c('WAT', 'UCL'))] %>% st_as_sf
p2 <- p[!(p$POLYTYPE %in% c('WAT', 'UCL'))] %>% as.data.frame

# calculate msi
p2$msi <- p2$perim/sqrt(pi * p2$area)

# load original raster input file
ras <- rast('D:/ontario_inventory/segmentation/spl_stack_mask_wat_ucl_polytype.tif')

# rename bands
names(ras) <- c('p95', 'cc', 'cv')

# extract pixel values
pvals <- exact_extract(ras, p2_sf)

# calculate SSE
sse <- sapply(pvals, FUN = function(x){
  
  # subset values based on coverage fraction
  x %<>% filter(coverage_fraction >= 0.5)
  
  p95_mean <- mean(x$p95, na.rm = T)
  cc_mean <- mean(x$cc, na.rm = T)
  cv_mean <- mean(x$cv, na.rm = T)
  
  return(c(sum((x$p95 - p95_mean)^2, na.rm = T),
           sum((x$cc - cc_mean)^2, na.rm = T),
           sum((x$cv - cv_mean)^2, na.rm = T)))
})

# transpose
sse <- t(sse)

# calculate final sums
sse <- colSums(sse)

# unlist values
pvals2 <- do.call(rbind, pvals)

# subset values based on coverage fraction
pvals2 %<>% filter(coverage_fraction >= 0.5)

# calculate global mean values
p95_mean <- mean(pvals2$p95, na.rm = T)
cc_mean <- mean(pvals2$cc, na.rm = T)
cv_mean <- mean(pvals2$cv, na.rm = T)

rm(pvals2)

# calculate SST
sst <- sapply(pvals, FUN = function(x){
  # subset values based on coverage fraction
  x %<>% filter(coverage_fraction >= 0.5)
  
  return(c(sum((x$p95 - p95_mean)^2, na.rm = T),
           sum((x$cc - cc_mean)^2, na.rm = T),
           sum((x$cv - cv_mean)^2, na.rm = T)))
})

# transpose
sst <- t(sst)

# calculate final sums
sst <- colSums(sst)

# calculate r2 values
r2_p95 <- 1 - (sse[1]/sst[1]) %>% round(3)
r2_cc <- 1 - (sse[2]/sst[2]) %>% round(3)
r2_cv <- 1 - (sse[3]/sst[3]) %>% round(3)
r2_all <- (sum(r2_p95, r2_cc, r2_cv) / 3) %>% round(3)

# create dataframe with values wanted
ms_df <- data.frame(file = 'FRI',
                    min_pix = (min(p2$area / 400)),
                    max_pix = (max(p2$area / 400)),
                    mean_pix = (mean(p2$area / 400)),
                    med_pix = (median(p2$area / 400)),
                    num_poly = NROW(p2),
                    mean_area = mean(p2$area),
                    se_area = se(p2$area),
                    sd_area = sd(p2$area),
                    mean_perim = mean(p2$perim),
                    se_perim = se(p2$perim),
                    sd_perim = sd(p2$perim),
                    mean_p_a = mean(p2$p_to_a),
                    se_p_a = se(p2$p_to_a),
                    sd_p_a = sd(p2$p_to_a),
                    mean_msi = mean(p2$msi),
                    se_msi = se(p2$msi),
                    sd_msi = sd(p2$msi),
                    r2_p95 <- r2_p95,
                    r2_cc <- r2_cc,
                    r2_cv <- r2_cv,
                    r2_all <- r2_all)

# round numeric columns
ms_df %<>% 
  mutate_at(c('min_pix',
              'max_pix', 
              'mean_pix', 
              'med_pix'),
            function(x) round(x)) %>%
  mutate_at(c('mean_area', 
              'se_area', 
              'sd_area',
              'mean_perim',
              'se_perim',
              'sd_perim'), 
            function(x) round(x, 2)) %>%
  mutate_at(c('mean_p_a', 
              'se_p_a', 
              'sd_p_a', 
              'mean_msi', 
              'se_msi', 
              'sd_msi'), 
            function(x) round(x, 4))

# bind df
df <- rbind(df, ms_df)

# write df as csv
write.csv(df,
          file = 'D:/ontario_inventory/segmentation/results_paper/summary_stats_grm_2022_04_28.csv',
          row.names = F)
