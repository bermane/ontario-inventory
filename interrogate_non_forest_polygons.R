# this code looks at the attributes of non-forest polygons in the FRI datasets
# to better understand height, canopy, cover, and other important variables.
# this will assist in deciding what to do with these polygons throughout
# the process

# There is also a comparison to landcover classes from the VLCE
# To see how well they match and also assess an updated LC dataset
# since the FRI is not updated regularly

# load packages
library(terra)
library(tidyverse)
library(sf)
library(exactextractr)

##############################
###INTERROGATE FRI POLYGONS###
##############################

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

# load values extracted without edge pixels
load('D:/ontario_inventory/dat/dat_fri_100.RData')

# move dat so can load values extracted without edge pixels
dat_100 <- dat
rm(dat)

# load values extracted without edge pixels 95%
load('D:/ontario_inventory/dat/dat_fri_95.RData')

# move dat so can load values extracted without edge pixels
dat_95 <- dat
rm(dat)

# Calculate summary stats by POLYTYPE compare 100%/95%
# EFI variables are not useful outside of forest polytype bc modeled
# from forest plots
dat_100 %>% group_by(POLYTYPE) %>% dplyr::summarize(mean_HT = mean(HT, na.rm = T),
                                                mean_p95 = mean(p95, na.rm = T),
                                                mean_CC = mean(CC, na.rm = T),
                                                mean_perc2m = mean(cc, na.rm = T),
                                                mean_BA = mean(BA, na.rm = T),
                                                mean_ba_efi = mean(ba, na.rm = T),
                                                mean_age = mean(AGE, na.rm = T),
                                                min_p95 = min(p95, na.rm = T))

dat_95 %>% group_by(POLYTYPE) %>% dplyr::summarize(mean_HT = mean(HT, na.rm = T),
                                                    mean_p95 = mean(p95, na.rm = T),
                                                    mean_CC = mean(CC, na.rm = T),
                                                    mean_perc2m = mean(cc, na.rm = T),
                                                    mean_BA = mean(BA, na.rm = T),
                                                    mean_ba_efi = mean(ba, na.rm = T),
                                                    mean_age = mean(AGE, na.rm = T),
                                                    min_p95 = min(p95, na.rm = T))

summary <- dat_100 %>% group_by(POLYTYPE) %>% dplyr::summarize(median_HT = median(HT, na.rm = T),
                                                median_p95 = median(p95, na.rm = T),
                                                median_CC = median(CC, na.rm = T),
                                                median_perc2m = median(cc, na.rm = T),
                                                median_BA = median(BA, na.rm = T),
                                                median_ba_efi = median(ba, na.rm = T),
                                                median_age = median(AGE, na.rm = T),
                                                min_p95 = min(p95, na.rm = T))

#write to disk
write.csv(summary, 
          file = 'D:/ontario_inventory/imputation/distributions/polytype_stats/polytype_summary_stats.csv',
          row.names = F)

dat_95 %>% group_by(POLYTYPE) %>% dplyr::summarize(median_HT = median(HT, na.rm = T),
                                                median_p95 = median(p95, na.rm = T),
                                                median_CC = median(CC, na.rm = T),
                                                median_perc2m = median(cc, na.rm = T),
                                                median_BA = median(BA, na.rm = T),
                                                median_ba_efi = median(ba, na.rm = T),
                                                median_age = median(AGE, na.rm = T),
                                                min_p95 = min(p95, na.rm = T))

# create plots of stats by polytype
ggplot(dat_100, aes(x = HT, color = POLYTYPE)) +
  geom_density() + lims(y = c(0, 0.5)) +
  theme_bw() + ggtitle("Height by Polytype Density Plot") +
  xlab("Height (m)") + ylab("Density")
ggsave(filename = 'D:/ontario_inventory/imputation/distributions/polytype_stats/ht_plot.png')

ggplot(dat_100, aes(x = p95, color = POLYTYPE)) +
  geom_density() + lims(y = c(0, 0.5)) +
  theme_bw() + ggtitle("P95 by Polytype Density Plot") +
  xlab("Height (m)") + ylab("Density")
ggsave(filename = 'D:/ontario_inventory/imputation/distributions/polytype_stats/p95_plot.png')

ggplot(dat_100, aes(x = CC, color = POLYTYPE)) +
  geom_density() + lims(y = c(0, 0.1)) +
  theme_bw() + ggtitle("Interp Canopy Cover by Polytype Density Plot") +
  xlab("Canopy Cover (%)") + ylab("Density")
ggsave(filename = 'D:/ontario_inventory/imputation/distributions/polytype_stats/cc_plot.png')

ggplot(dat_100, aes(x = cc, color = POLYTYPE)) +
  geom_density() + lims(y = c(0, 0.1)) +
  theme_bw() + ggtitle("% 2m Returns by Polytype Density Plot") +
  xlab("2 Meter Returns (%)") + ylab("Density")
ggsave(filename = 'D:/ontario_inventory/imputation/distributions/polytype_stats/2mrtn_plot.png')

ggplot(dat_100, aes(x = BA, color = POLYTYPE)) +
  geom_density() + lims(y = c(0, 0.2)) +
  theme_bw() + ggtitle("Interp Basal Area by Polytype Density Plot") +
  xlab(expression("Basal Area (m"^2*"/ha)")) + ylab("Density")
ggsave(filename = 'D:/ontario_inventory/imputation/distributions/polytype_stats/ba_plot.png')

ggplot(dat_100, aes(x = ba, color = POLYTYPE)) +
  geom_density() + lims(y = c(0, 0.2)) +
  theme_bw() + ggtitle("EFI Basal Area by Polytype Density Plot") +
  xlab(expression("Basal Area (m"^2*"/ha)")) + ylab("Density")
ggsave(filename = 'D:/ontario_inventory/imputation/distributions/polytype_stats/ba_efi_plot.png')

ggplot(dat_100, aes(x = AREA, color = POLYTYPE)) +
  geom_density() + lims(x = c(0,200000)) +
  theme_bw() + ggtitle("Polygon Area by Polytype Density Plot") +
  xlab(expression("Area (m"^2*")")) + ylab("Density")
ggsave(filename = 'D:/ontario_inventory/imputation/distributions/polytype_stats/poly_area_plot.png')

###################################
###INTERROGATE LANDCOVER CLASSES###
###################################

# clear workspace
rm(list=ls())

# load landcover data
lc <- rast('D:/ontario_inventory/VLCE/CA_forest_VLCE2_2018_CLIPPED.tif')

# change forests to one class
lc[lc == 220] <- 210
lc[lc == 230] <- 210

# convert raster to polygon
poly_lc <- as.polygons(lc)

# disagg polygons
# poly_lc <- disagg(poly_lc)

# convert to data frame
dat <- as.data.frame(poly_lc)

# load LiDAR datasets we need to extract over polygons
# create named vector with variable names and data links
# lidar <- c('cc' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif',
#            'avg' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_avg.tif',
#            'max' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_max.tif',
#            'p95' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_p95.tif',
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

lidar <- c('cc' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif',
           'p95' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_p95.tif',
           'lor' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_lor.tif',
           'ba' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_ba_ha.tif')

# load LiDAR rasters as raster stack
lidar_ras <- rast(lidar)

# project poly to crs of raster
poly_ras <- project(poly_lc, lidar_ras)

# convert to sf
poly_ras <- st_as_sf(poly_ras)

#extract median values
vec <- exact_extract(lidar_ras, poly_ras, function(values, coverage_fraction){
  values <- values[coverage_fraction == 1,]
  apply(values, 2, function(x) median(x, na.rm = T))
})

# transpose matrix and make data.frame
vec <- t(vec) %>% as.data.frame

# change column names
colnames(vec) <- names(lidar)

# add new column into dat
dat <- cbind(dat, vec)

# change first col name
colnames(dat)[1] <- 'lc'

# change values of lc to classes
dat$lc <- as.character(dat$lc)
dat$lc <- plyr::revalue(dat$lc, c('NaN' = NA,
                                  '100' = 'Herbs',
                                  '210' = 'Forest',
                                  '81' = 'Wetland-Treed',
                                  '20' = 'Water',
                                  '80' = 'Wetland',
                                  '50' = 'Shrubland',
                                  '33' = 'Exposed/Barren Land',
                                  '32' = 'Rock/Rubble',
                                  '40' = 'Bryoids'))

# save extracted dataframe for fast rebooting
save(dat, file = 'D:/ontario_inventory/dat/dat_lc_95.RData')

dat$lc <- as.factor(dat$lc)

# check summary stats
# just look at median doesn't matter only one value for each (but median calc)
summary <- dat %>% group_by(lc) %>% dplyr::summarize(median_p95 = median(p95, na.rm = T),
                                          median_perc2m = median(cc, na.rm = T),
                                          median_ba_efi = median(ba, na.rm = T),
                                          median_lor = median(lor, na.rm = T))

# write to disk
write.csv(summary, 
          file = 'D:/ontario_inventory/imputation/distributions/polytype_stats/lc_summary_stats.csv',
          row.names = F)

# create plots of stats by lc type
# no plots at the moment because just took median across whole polygon
