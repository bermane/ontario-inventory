# This code runs an imputation algorithm on LiDAR attributes connecting them to interpreter
# derived response attributes resulting in a model that can be applied to forest stands
# segmented from LiDAR rasters

# The purpose is to compare the distributions of various attributes (HT, CC, BA, AGE, POLYTYPE, WG)
# from FRI polygons, LiDAR attributes in FRI polygons, imputed vars into FRI polys, imputed vars into LiDAR polys

# this new version only runs forested pixels since we are only really interested in modelling
# and imputing attributes in forests
# it uses data screened for 10 and 90th percentile of residuals

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

# SPL vars only

# change all non-numeric variables to factor
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], 
                                         as.factor)

# build dataframe of reference variables
ref_vars <- c('cc', 'p95', 'avg', 'qav', 'cv', 'kur', 'max', 'ske')
x <- dat[, ref_vars]

# build dataframe of target variables
tar_vars <- c('HT', 'WG')
y <- dat[, tar_vars]

# build dataframe of ancillary data
anci <- dat[, !(colnames(dat) %in% ref_vars)]

# run imputation algorithm
rf <- yai(x = x, y = y, method = 'randomForest', ntree = 100*NCOL(y))

# plot output
plot(rf, vars = yvars(rf))

####################################################################
###LOAD LIDAR DERIVED POLYGONS AND POPULATE WITH LIDAR ATTRIBUTES###
####################################################################

# load lidar derived polygons
poly_lidar <- vect('D:/ontario_inventory/segmentation/ms_10_10_100_for_only_agg_na.shp')

# load LiDAR attributes
load('D:/ontario_inventory/imputation/seg_df_ms_10_10_100_for_only_agg_na.RData')

# we need to set rownames of lidar polygons so they don't overlap reference data
# find final rowname from reference and add 1
start_r <- as.numeric(rownames(dat)[NROW(dat)]) + 1
end_r <- start_r + NROW(dat_lidar) - 1

# change row names so don't overlap with reference dataset
rownames(dat_lidar) <- start_r:end_r

# run imputation over new polygons
imp_lidar <- predict(object = rf,
                     newdata = dat_lidar,
                     ancillaryData = anci,
                     observed = F)

# find any rows that were removed
rm_rows <- setdiff(rownames(dat_lidar), rownames(imp_lidar))

# re add any rows that were removed
if(length(rm_rows) != 0){
  add_rows <- imp_lidar[1:length(rm_rows),]
  add_rows[] <- NA
  rownames(add_rows) <- rm_rows
  imp_lidar <- rbind(imp_lidar, add_rows)
}

# re order the rows of imp_lidar to match dat_lidar
imp_lidar <- imp_lidar[order(as.numeric(row.names(imp_lidar))),]

# remove columns from imputed dataset if they exist from lidar
imp_lidar <- imp_lidar[, !(colnames(imp_lidar) %in% colnames(dat_lidar))]

# add additional columns from lidar dataset that were not imputed
imp_lidar <- imp_lidar %>% add_column(dat_lidar)

# change back to original row names
rownames(imp_lidar) <- 1:NROW(imp_lidar)

# repopulate polygons with data
values(poly_lidar) <- as.data.frame(imp_lidar)

# write poly_lidar without big missing NA polygon
poly_lidar2 <- poly_lidar[is.nan(poly_lidar$meanB0) == F,]

# write polygons with imputed data
writeVector(poly_lidar2, filename = 'D:/ontario_inventory/imputation/distributions_for_only/vector/ms_10_10_100_for_only_imp_10_perc.shp',
            overwrite = T)

# write dataframe with imputed data
write.csv(imp_lidar, file = 'D:/ontario_inventory/imputation/distributions_for_only/vector/ms_10_10_100_for_only_imp_10_perc.csv')

#################################
###COMPARISON WITH FRI DATASET###
#################################

# apply the imputation algorithm to the training/forest data
# we are not comparing to the other data because it doesn't 
# meet our criteria for clean forested polygons

# run imputation over new polygons
# we want to keep observed cols to see relationship
imp_full <- impute(object = rf,
              ancillaryData = anci,
              observed = T)

# add reference data
imp_full <- imp_full %>% add_column(x)

# create new polygon dataset for imp data
poly_imp <- poly

# load as df
df <- as.data.frame(poly_imp)

# subset polygons
poly_imp <- poly_imp[rownames(df) %in% rownames(imp_full)]

# fill with data
values(poly_imp) <- as.data.frame(imp_full)
rm(df)

# write polygons with imputed data
writeVector(poly_imp, filename = 'D:/ontario_inventory/imputation/distributions_for_only/vector/fri_polygons_for_only_imp_10_perc.shp',
            overwrite = T)

# write dataframe with imputed data
write.csv(imp_full, file = 'D:/ontario_inventory/imputation/distributions_for_only/vector/fri_polygons_for_only_imp_10_perc.csv')

#########################################
###COMPARE DISTRIBUTIONS OF ATTRIBUTES###
#########################################

#need to choose a color palette colorblind friendly
#https://personal.sron.nl/~pault/#sec:qualitative
#'#4477AA', '#000000', '#228833', '#CCBB44', '#EE6677', '#AA3377'

colscale_fri <- scale_color_manual(name = "Dataset", values = c('FRI Polygons' = '#4477AA', 'LiDAR Imp Polygons' = '#EE6677',
                                                                'FRI Imp Polygons' = '#228833'))

colscale_fri_full <- scale_color_manual(name = "Dataset", values = c('FRI Full Polygons' = '#CCBB44', 
                                                                     'FRI Screened Polygons' = '#4477AA', 
                                                                     'FRI Imp Polygons' = '#228833',
                                                                     'LiDAR Imp Polygons' = '#EE6677'))


# set line size
sz <- 1.2

# Interp Height and Lorey's Height
ggplot() +
  geom_density(data = imp_full, mapping = aes(x = lor, color = 'FRI (Lor Height)'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = lor, color = 'LiDAR (Lor Height)'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = HT.o, color = 'FRI (HT)'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = HT, color = 'LiDAR Imp (HT)'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = HT, color = 'FRI Imp (HT)'), size = sz) +
  theme_bw() + 
  ggtitle("Height Density Plot") +
  xlab("Height (m)") + ylab("Density") + 
  scale_color_manual(name = "Dataset", values = c('FRI (Lor Height)' = '#CCBB44', 'LiDAR (Lor Height)' = '#AA3377',
                                                  'FRI (HT)' = '#4477AA', 'LiDAR Imp (HT)' = '#EE6677',
                                                  'FRI Imp (HT)' = '#228833')) +
  theme(text = element_text(size = 20))
ggsave(filename = 'D:/ontario_inventory/imputation/distributions_for_only/height_for_only_10perc.png')

# Interp Height Only
ggplot() +
  geom_density(data = imp_full, mapping = aes(x = HT.o, color = 'FRI Polygons'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = HT, color = 'LiDAR Imp Polygons'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = HT, color = 'FRI Imp Polygons'), size = sz) +
  theme_bw() + 
  ggtitle("Height Density Plot") +
  xlab("Height (m)") + ylab("Density") + 
  colscale_fri +
  theme(text = element_text(size = 20))
ggsave(filename = 'D:/ontario_inventory/imputation/distributions_for_only/interp_height_for_only_10perc.png')

# Interp Height Only with full FRI also
ggplot() +
  geom_density(data = dat_orig, mapping = aes(x = HT, color = 'FRI Full Polygons'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = HT.o, color = 'FRI Screened Polygons'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = HT, color = 'LiDAR Imp Polygons'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = HT, color = 'FRI Imp Polygons'), size = sz) +
  theme_bw() + 
  ggtitle("Height Density Plot") +
  xlab("Height (m)") + ylab("Density") + 
  colscale_fri_full +
  theme(text = element_text(size = 20))
ggsave(filename = 'D:/ontario_inventory/imputation/distributions_for_only/interp_height_for_only_10perc_full.png')

# Canopy Cover two metrics
ggplot() +
  geom_density(data = imp_full, mapping = aes(x = cc, color = 'FRI (% 2m Rtn)'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = cc, color = 'LiDAR (% 2m Rtn)'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = CC.o, color = 'FRI (CC)'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = CC, color = 'LiDAR Imp (CC)'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = CC, color = 'FRI Imp (CC)'), size = sz) +
  theme_bw() + 
  ggtitle("Canopy Cover Density Plot") +
  xlab("Cover (%)") + ylab("Density") + 
  scale_color_manual(name = "Dataset", values = c('FRI (% 2m Rtn)' = '#CCBB44', 'LiDAR (% 2m Rtn)' = '#AA3377',
                                                  'FRI (CC)' = '#4477AA', 'LiDAR Imp (CC)' = '#EE6677',
                                                  'FRI Imp (CC)' = '#228833')) +
  theme(text = element_text(size = 20))
ggsave(filename = 'D:/ontario_inventory/imputation/distributions_for_only/canopy_cover_for_only_10perc.png')

# Interp Canopy Cover Only
ggplot() +
  geom_density(data = imp_full, mapping = aes(x = CC.o, color = 'FRI Polygons'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = CC, color = 'LiDAR Imp Polygons'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = CC, color = 'FRI Imp Polygons'), size = sz) +
  theme_bw() + 
  ggtitle("Canopy Cover Density Plot") +
  xlab("Cover (%)") + ylab("Density") + 
  colscale_fri +
  theme(text = element_text(size = 20))
ggsave(filename = 'D:/ontario_inventory/imputation/distributions_for_only/interp_canopy_cover_for_only_10perc.png')

# Interp Canopy Cover Only with Full FRI
ggplot() +
  geom_density(data = dat_orig, mapping = aes(x = CC, color = 'FRI Full Polygons'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = CC.o, color = 'FRI Screened Polygons'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = CC, color = 'LiDAR Imp Polygons'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = CC, color = 'FRI Imp Polygons'), size = sz) +
  theme_bw() + 
  ggtitle("Canopy Cover Density Plot") +
  xlab("Cover (%)") + ylab("Density") + 
  colscale_fri_full +
  theme(text = element_text(size = 20))
ggsave(filename = 'D:/ontario_inventory/imputation/distributions_for_only/interp_canopy_cover_for_only_10perc_full.png')

# Basal Area
ggplot() +
  geom_density(data = imp_full, mapping = aes(x = ba, color = 'FRI (EFI gen)'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = ba, color = 'LiDAR (EFI gen)'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = BA.o, color = 'FRI (BA)'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = BA, color = 'LiDAR Imp (BA)'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = BA, color = 'FRI Imp (BA)'), size = sz) +
  theme_bw() + 
  ggtitle("Basal Area Density Plot") +
  xlab(expression("Basal Area (m"^2*"/ha)")) + ylab("Density") + 
  scale_color_manual(name = "Dataset", values = c('FRI (EFI gen)' = '#CCBB44', 'LiDAR (EFI gen)' = '#AA3377',
                                                  'FRI (BA)' = '#4477AA', 'LiDAR Imp (BA)' = '#EE6677',
                                                  'FRI Imp (BA)' = '#228833')) +
  theme(text = element_text(size = 20))
ggsave(filename = 'D:/ontario_inventory/imputation/distributions_for_only/basal_area_for_only_10perc.png')

# Interp Basal Area Only
ggplot() +
  geom_density(data = imp_full, mapping = aes(x = BA.o, color = 'FRI Polygons'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = BA, color = 'LiDAR Imp Polygons'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = BA, color = 'FRI Imp Polygons'), size = sz) +
  theme_bw() + 
  ggtitle("Basal Area Density Plot") +
  xlab(expression("Basal Area (m"^2*"/ha)")) + ylab("Density") + 
  colscale_fri +
  theme(text = element_text(size = 20))
ggsave(filename = 'D:/ontario_inventory/imputation/distributions_for_only/interp_basal_area_for_only_10perc.png')

# Interp Basal Area Only with FUll FRI
ggplot() +
  geom_density(data = dat_orig, mapping = aes(x = BA, color = 'FRI Full Polygons'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = BA.o, color = 'FRI Screened Polygons'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = BA, color = 'LiDAR Imp Polygons'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = BA, color = 'FRI Imp Polygons'), size = sz) +
  theme_bw() + 
  ggtitle("Basal Area Density Plot") +
  xlab(expression("Basal Area (m"^2*"/ha)")) + ylab("Density") + 
  colscale_fri_full +
  theme(text = element_text(size = 20))
ggsave(filename = 'D:/ontario_inventory/imputation/distributions_for_only/interp_basal_area_for_only_10perc_full.png')

# Age
ggplot() +
  geom_density(data = imp_full, mapping = aes(x = AGE.o, color = 'FRI Polygons'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = AGE, color = 'LiDAR Imp Polygons'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = AGE, color = 'FRI Imp Polygons'), size = sz) +
  theme_bw() + 
  ggtitle("Interpreter Derived Age Density Plot") +
  xlab("Age (years)") + ylab("Density") + 
  colscale_fri +
  theme(text = element_text(size = 20))
ggsave(filename = 'D:/ontario_inventory/imputation/distributions_for_only/age_for_only_10perc.png')

# Age with Full FRI
ggplot() +
  geom_density(data = dat_orig, mapping = aes(x = AGE, color = 'FRI Full Polygons'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = AGE.o, color = 'FRI Screened Polygons'), size = sz) +
  geom_density(data = imp_lidar, mapping = aes(x = AGE, color = 'LiDAR Imp Polygons'), size = sz) +
  geom_density(data = imp_full, mapping = aes(x = AGE, color = 'FRI Imp Polygons'), size = sz) +
  theme_bw() + 
  ggtitle("Interpreter Derived Age Density Plot") +
  xlab("Age (years)") + ylab("Density") + 
  colscale_fri_full +
  theme(text = element_text(size = 20))
ggsave(filename = 'D:/ontario_inventory/imputation/distributions_for_only/age_for_only_10perc_full.png')

# Working group
# better as freq table
# load full data categories
# again don't include water polygons since skews the distribution
wg_freq <- plyr::count(imp_full, 'WG.o')

# change col names
wg_freq <- rename(wg_freq, WG = WG.o, fri_polygons = freq)

# add other polys
add1 <- plyr::count(imp_lidar, 'WG')
add1 <- drop_na(add1)
add2 <- plyr::count(imp_full, 'WG')
add2 <- drop_na(add2)

wg_freq$add1 <- 0
wg_freq$add2 <- 0

wg_freq$add1[wg_freq$WG %in% add1$WG] <- add1$freq
wg_freq$add2[wg_freq$WG %in% add2$WG] <- add2$freq

# change col names
wg_freq <- rename(wg_freq, lidar_imp_polygons = add1, fri_imp_polygons = add2)

# add percentages
wg_freq <- add_column(wg_freq, fri_poly_perc = wg_freq$fri_polygons/sum(wg_freq$fri_polygons)*100,
                      lidar_imp_poly_perc = wg_freq$lidar_imp_polygons/sum(wg_freq$lidar_imp_polygons)*100,
                      fri_imp_poly_perc = wg_freq$fri_imp_polygons/sum(wg_freq$fri_imp_polygons)*100)

# round
wg_freq[,(NCOL(wg_freq)-2):NCOL(wg_freq)] <- round(wg_freq[,(NCOL(wg_freq)-2):NCOL(wg_freq)], 2)

# export as csv
write.csv(wg_freq, 
          file = 'D:/ontario_inventory/imputation/distributions_for_only/working_group_10_perc.csv', 
          row.names = F)

##############################################################
###CALCULATE PERFORMANCE METRICS OF FRI OBSERVED VS IMPUTED###
##############################################################

# calculate rmse different variables
rmse <- function(obs, est) sqrt(mean((obs - est)^2))

rmse_df <- data.frame(AGE = c(rmse(imp_full$AGE.o, imp_full$AGE),
                              rmse(imp_full$AGE.o, imp_full$AGE)/sd(imp_full$AGE.o)),
                      HT = c(rmse(imp_full$HT.o, imp_full$HT),
                             rmse(imp_full$HT.o, imp_full$HT)/sd(imp_full$HT.o)),
                      BA = c(rmse(imp_full$BA.o, imp_full$BA),
                             rmse(imp_full$BA.o, imp_full$BA)/sd(imp_full$BA.o)),
                      CC = c(rmse(imp_full$CC.o, imp_full$CC),
                             rmse(imp_full$CC.o, imp_full$CC)/sd(imp_full$CC.o)))

rownames(rmse_df) <- c('raw_rmsd', 'scaled_rmsd')

# write table
write.csv(rmse_df, file = 'D:/ontario_inventory/imputation/distributions_for_only/rmsd_10_perc.csv')

# load function to calculate kappa coefficient
kappa <- function(m) {
  N <- sum(m)
  No <- sum(diag(m))
  Ne <- 1 / N * sum(colSums(m) * rowSums(m))
  return( (No - Ne) / (N - Ne) )
}

# create df of WG
wg <- data.frame(obs = imp_full$WG.o,
                 est = imp_full$WG)

# make estimate levels match obs levels
levels(wg$est) <- c(levels(wg$est), levels(wg$obs)[!(levels(wg$obs) %in% levels(wg$est))])
wg$est <- factor(wg$est, levels = levels(wg$obs))

# create column of match or not
wg$match <- wg$obs == wg$est

# total percent of matching WG
NROW(wg[wg$match == T,])/NROW(wg)

# check count of different working groups
plyr::count(wg, 'obs')
plyr::count(wg, 'est')

# build accuracy table
accmat_wg <- table("pred" = wg$est, "ref" = wg$obs)

# UA
ua_wg <- diag(accmat_wg) / rowSums(accmat_wg) * 100

# PA
pa_wg <- diag(accmat_wg) / colSums(accmat_wg) * 100

# OA
oa_wg <- sum(diag(accmat_wg)) / sum(accmat_wg) * 100

# build confusion matrix
accmat_wg_ext <- addmargins(accmat_wg)
accmat_wg_ext <- rbind(accmat_wg_ext, "Users" = c(pa_wg, NA))
accmat_wg_ext <- cbind(accmat_wg_ext, "Producers" = c(ua_wg, NA, oa_wg))
#colnames(accmat_wg_ext) <- c(levels(as.factor(shp.train$classes)), "Sum", "PA")
#rownames(accmat_wg_ext) <- c(levels(as.factor(shp.train$classes)), "Sum", "UA")
accmat_wg_ext <- round(accmat_wg_ext, digits = 1)
dimnames(accmat_wg_ext) <- list("Prediction" = colnames(accmat_wg_ext),
                                "Reference" = rownames(accmat_wg_ext))
class(accmat_wg_ext) <- "table"
accmat_wg_ext

# calculate kappa coefficient
kap <- kappa(accmat_wg) %>% round(2)

# add kappa to matrix
accmat_wg_ext <- rbind(accmat_wg_ext, c(kap, rep('', NCOL(accmat_wg_ext) - 1)))
rownames(accmat_wg_ext)[NROW(accmat_wg_ext)] <- 'Kappa'

# write to csv
write.csv(accmat_wg_ext, file = 'D:/ontario_inventory/imputation/distributions_for_only/accmat_wg_10_perc.csv')
