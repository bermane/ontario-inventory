# this code extracts final polygon summary stats from the best segmentation outputs

# load packages
library(terra)
library(tidyverse)
library(janitor)
library(magrittr)
library(sf)
library(lwgeom)
library(exactextractr)

# create list of polygon files, names and parameters
# file <- c('D:/ontario_inventory/segmentation/grm/shp/grm_10_01_05.shp',
#           'D:/ontario_inventory/segmentation/grm/shp/grm_5_01_05.shp',
#           'D:/ontario_inventory/segmentation/grm/shp/grm_15_02_05.shp',
#           'D:/ontario_inventory/segmentation/python/shp/segments_slic_1000k_01_short_add_pt.shp',
#           'D:/ontario_inventory/segmentation/python/shp/segments_slic0_1000k_01_short_add_pt.shp',
#           'D:/ontario_inventory/segmentation/python/shp/segments_slic_400k_01_add_pt.shp',
#           'D:/ontario_inventory/segmentation/mask_wat_ucl/ms_10_10_100_add_pt.shp',
#           'D:/ontario_inventory/segmentation/mask_wat_ucl/ms_5_15_50_add_pt.shp')
# alg <- c('GRM',
#          'GRM',
#          'GRM',
#          'SLIC',
#          'SLIC',
#          'SLIC',
#          'MS',
#          'MS')
# param <- c('10_0.1_0.5',
#            '5_0.1_0.5',
#            '15_0.2_0.5',
#            '1000k_01',
#            '1000k_01_zero',
#            '400k_01',
#            '10_10_100',
#            '5_15_50')
# name <- c('10_01_05',
#           '5_01_05',
#           '15_02_05',
#           '1000k_01',
#           '1000k_01_zero',
#           '400k_01',
#           '10_10_100',
#           '5_15_50')

# create list of polygon files, names and parameters
file <- c('D:/ontario_inventory/segmentation/grm/shp/grm_10_01_05.shp',
          'D:/ontario_inventory/segmentation/python/shp/segments_slic_1000k_01_short_add_pt.shp',
          'D:/ontario_inventory/segmentation/mask_wat_ucl/ms_10_10_100_add_pt.shp')
alg <- c('GRM',
         'SLIC',
         'MS')
param <- c('10_0.1_0.5',
           '1000k_01',
           '10_10_100')
name <- c('10_01_05',
          '1000k_01',
          '10_10_100')

# create standard error function
se <- function(x) sd(x)/sqrt(length(x))

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
  ms_df <- data.frame(alg = alg[i],
                      param = param[i],
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
  
  #############################
  ### CALCULATE FINAL PLOTS ###
  #############################
  
  # plot density
  ggplot(data.frame(nbPixels = p2$nbPixels), aes(x = nbPixels)) +
    geom_density(fill = 'grey') +
    xlim(c(0,500)) +
    ylim(c(0, 0.015)) +
    geom_vline(aes(xintercept = median(nbPixels)), 
               linetype = "dashed", size = 0.6) +
    theme_bw() +
    xlab('Number of Pixels') +
    ylab('Density') +
    ggtitle(str_c(alg[i])) +
    theme(text = element_text(size = 20))
  
  ggsave(str_c('D:/ontario_inventory/segmentation/results_paper/plots/', alg[i], '_final_pix_dens.png'),
         width = 2100, height = 2100, units = 'px')
  
  # # plot perimeter to area ratio
  # ggplot(data.frame(p_to_a = as.numeric(p2$p_to_a)), aes(x = p_to_a)) +
  #   geom_density() +
  #   xlim(c(0,0.2)) +
  #   ylim(c(0, 100)) +
  #   geom_vline(aes(xintercept = median(p_to_a)),
  #              linetype = "dashed", size = 0.6) +
  #   theme_bw() +
  #   xlab('Perimeter to Area Ratio') +
  #   ylab('Density') +
  #   ggtitle(str_c(alg[i], ' ', param[i])) +
  #   theme(text = element_text(size = 20))
  # 
  # ggsave(str_c('D:/ontario_inventory/segmentation/results_paper/plots/', alg[i], '_', name[i], '_ptoa_dens.png'),
  #        width = 2100, height = 2100, units = 'px')
  
  # plot shape index
  ggplot(data.frame(msi = as.numeric(p2$msi)), aes(x = msi)) +
    geom_density(fill = 'grey') +
    xlim(c(0,8)) +
    ylim(c(0, 1.5)) +
    geom_vline(aes(xintercept = median(msi)),
               linetype = "dashed", size = 0.6) +
    theme_bw() +
    xlab('Shape Index') +
    ylab('Density') +
    ggtitle(str_c(alg[i])) +
    theme(text = element_text(size = 20))
  
  ggsave(str_c('D:/ontario_inventory/segmentation/results_paper/plots/', alg[i], '_final_shape_dens.png'),
         width = 2100, height = 2100, units = 'px')
  
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
ms_df <- data.frame(alg = 'FRI',
                    param = 'FRI',
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

# # write df as csv
# write.csv(df,
#           file = 'D:/ontario_inventory/segmentation/results_paper/summary_stats_2022_04_28.csv',
#           row.names = F)

#############################
### CALCULATE FINAL PLOTS ###
#############################

# plot density
ggplot(data.frame(nbPixels = p2$area / 400), aes(x = nbPixels)) +
  geom_density(fill = 'grey') +
  xlim(c(0,500)) +
  ylim(c(0, 0.015)) +
  geom_vline(aes(xintercept = median(nbPixels)), 
             linetype = "dashed", size = 0.6) +
  theme_bw() +
  xlab('Number of Pixels') +
  ylab('Density') +
  ggtitle('FRI') +
  theme(text = element_text(size = 20))

ggsave(str_c('D:/ontario_inventory/segmentation/results_paper/plots/fri_final_pix_dens.png'),
       width = 2100, height = 2100, units = 'px')

# # plot perimeter to area ratio
# ggplot(data.frame(p_to_a = as.numeric(p2$p_to_a)), aes(x = p_to_a)) +
#   geom_density() +
#   xlim(c(0,0.2)) +
#   ylim(c(0, 100)) +
#   geom_vline(aes(xintercept = median(p_to_a)),
#              linetype = "dashed", size = 0.6) +
#   theme_bw() +
#   xlab('Perimeter to Area Ratio') +
#   ylab('Density') +
#   ggtitle('FRI') +
#   theme(text = element_text(size = 20))
# 
# ggsave(str_c('D:/ontario_inventory/segmentation/results_paper/plots/fri_ptoa_dens.png'),
#        width = 2100, height = 2100, units = 'px')

# plot shape index
ggplot(data.frame(msi = as.numeric(p2$msi)), aes(x = msi)) +
  geom_density(fill = 'grey') +
  xlim(c(0,8)) +
  ylim(c(0, 1.5)) +
  geom_vline(aes(xintercept = median(msi)),
             linetype = "dashed", size = 0.6) +
  theme_bw() +
  xlab('Shape Index') +
  ylab('Density') +
  ggtitle('FRI') +
  theme(text = element_text(size = 20))

ggsave(str_c('D:/ontario_inventory/segmentation/results_paper/plots/fri_final_shape_dens.png'),
       width = 2100, height = 2100, units = 'px')
