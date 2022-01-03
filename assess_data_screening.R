# This code loads photo interpreted polygons from the Romeo Malette Forest in Ontario
# and applies a number of data screens to curate an accurate and reliable dataset

# load packages
library(terra)
library(tidyverse)

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

# convert to df
dat <- as.data.frame(poly)

###########################
###DATA SCREENING PART 1###
###########################

# A. Polygon area < 50 ha (500000 m^2)

# check number of polygons with area < 50 ha
#NROW(dat[dat$AREA < 500000,])

# subset data
#dat_1a <- dat[dat$AREA < 500000,]

# look AT HISTOGRAM OF VERY LARGE POLYGONS AND GET RID OF THEM and ONLY FOR POLYGONS!!!

# let's skip this skip for now.

# B. Within-polygon height and canopy cover coefficient of variation < 0.5

# load loreys height and canopy cover from LiDAR
cc <- rast('D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif')
lor <- rast('D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_lor.tif')

# project poly to crs of rasters
poly_cc <- project(poly, cc)
poly_lor <- project(poly, lor)

# extract coeff of variation -- only using cells whose center is within polygon
# chen did not use edge pixels -- here edge pixels are still included so can consider
# writing code to remove edge pixels from the calculation
cc_cv <- terra::extract(cc, poly_cc, fun = function(x){sd(x, na.rm = T)/mean(x, na.rm = T)})
lor_cv <- terra::extract(lor, poly_lor, fun = function(x){sd(x, na.rm = T)/mean(x, na.rm = T)})

# remove first column
cc_cv <- cc_cv[,2]
lor_cv <- lor_cv[,2]

# cbind coeff of variation into dat
dat <- cbind(dat, cc_cv, lor_cv)

# remove layers
rm(cc, lor, cc_cv, lor_cv, poly_cc, poly_lor)

# check number of polygons with cv < 0.5
NROW(dat[dat$cc_cv < 0.5,])
NROW(dat[dat$lor_cv < 0.5,])

# subset dat based on cv
dat_1b <- dat[dat$cc_cv < 0.5 & dat$lor_cv < 0.5,]

# C. number of land cover types within a polygon < 2

# load VLCE 2.0 landcover dataset from 2018
lc <- rast('D:/ontario_inventory/VLCE/CA_forest_VLCE2_2018_CLIPPED.tif')

# project poly to crs of raster
poly_lc <- project(poly, lc)

# extract lc types within each polygon
lc_poly <- terra::extract(lc, poly_lc)

# change col names
colnames(lc_poly) <- c('ID', 'COVER')

# define treed classes
#treed <- c('Mixed Wood', 'Broadleaf', 'Coniferous', 'Wetland-Treed')
treed <- c(230, 220, 210, 81)

# calculate number of unique lc types in each polygon and % of polygon covered by treed classes
lc_poly <- lc_poly %>% group_by(ID) %>%
  mutate(num_uniq = length(unique(COVER))) %>%
  mutate(perc_for = sum(COVER %in% treed)/length(COVER))

# drop duplicate rows
lc_poly <- lc_poly[duplicated(lc_poly$ID)==F,]

# add new columns into dat
dat <- dat %>% add_column(num_uniq = lc_poly$num_uniq,
                          perc_for = lc_poly$perc_for)

# rm variables
rm(lc, lc_poly, poly_lc, treed)

# check number of polygons with LC types < 2
# perhaps Chen meant LC types <= 2??
NROW(dat[dat$num_uniq < 2,])
NROW(dat[dat$num_uniq <= 2,])

# check number of polygons with forested classes > 50%
NROW(dat[dat$perc_for > .5,])

# subset dat based on lc info
dat_1c <- dat[dat$num_uniq <=2 & dat$perc_for > .5,]

###########################
###DATA SCREENING PART 2###
###########################

# first load all data, take median values, and enter into dat
# Load Lorey's Height, cc, and ba from LiDAR
lor <- rast('D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_lor.tif')
cc <- rast('D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif')
ba <- rast('D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_ba_ha.tif')

# project poly to crs of rasters
poly_lor <- project(poly, lor)
poly_cc <- project(poly, cc)
poly_ba <- project(poly, ba)

# extract median within each polygon
# chen did not use edge pixels -- here edge pixels are still included so can consider
# writing code to remove edge pixels from the calculation
lor_med <- terra::extract(lor, poly_lor, fun = function(x){median(x, na.rm = T)})
cc_med <- terra::extract(cc, poly_cc, fun = function(x){median(x, na.rm = T)})
ba_med <- terra::extract(ba, poly_ba, fun = function(x){median(x, na.rm = T)})

# add new column into dat
dat <- dat %>% add_column(lor_med = lor_med[,2],
                          cc_med = cc_med[,2],
                          ba_med = ba_med[,2])

# remove layers
rm(lor, lor_med, poly_lor, cc, cc_med, poly_cc, ba, ba_med, poly_ba)

# A. Linear regression of Height vs. Lorey's Height

# check values by POLYTYPE
table(dat$HT[dat$POLYTYPE == 'FOR'])
table(dat$HT[dat$POLYTYPE == 'WAT'])
table(dat$HT[dat$POLYTYPE == 'DAL'])
table(dat$HT[dat$POLYTYPE == 'GRS'])
table(dat$HT[dat$POLYTYPE == 'ISL'])
table(dat$HT[dat$POLYTYPE == 'UCL'])
table(dat$HT[dat$POLYTYPE == 'BSH'])
table(dat$HT[dat$POLYTYPE == 'OMS'])
table(dat$HT[dat$POLYTYPE == 'RCK'])
table(dat$HT[dat$POLYTYPE == 'TMS'])

# ALL HEIGHTS ARE ZERO OUTSIDE OF POLYTYPE == 'FOR'

# remove rows with missing values
dat_2a <- dat[is.na(dat$HT) == F & is.na(dat$lor_med) == F,]

# only keep land polygons
dat_2a <- dat_2a[dat_2a$POLYTYPE != 'WAT',]

# plot variables against each other
plot(dat_2a$lor_med, dat_2a$HT)
plot(dat_2a$lor_med[dat_2a$POLYTYPE == 'FOR'], dat_2a$HT[dat_2a$POLYTYPE == 'FOR'])

# run simple linear model
lm_ht <- lm(HT ~ lor_med, data = dat_2a)
summary(lm_ht)

# bind residuals to data
dat_2a <- dat_2a %>% add_column(lm_ht_resid = resid(lm_ht))

# find 10th and 90th percentile
lm_ht_perc <- quantile(dat_2a$lm_ht_resid, probs = c(.10, .90))

# remove lower 10 and upper 10 from dat
dat_2a_10_90 <- dat_2a[dat_2a$lm_ht_resid > lm_ht_perc[1] & dat_2a$lm_ht_resid < lm_ht_perc[2],]

# find 200th and 80th percentile
lm_ht_perc <- quantile(dat_2a$lm_ht_resid, probs = c(.20, .80))

# remove lower 10 and upper 10 from dat
dat_2a_20_80 <- dat_2a[dat_2a$lm_ht_resid > lm_ht_perc[1] & dat_2a$lm_ht_resid < lm_ht_perc[2],]

# find 30th and 70th percentile
lm_ht_perc <- quantile(dat_2a$lm_ht_resid, probs = c(.30, .70))

# remove lower 10 and upper 10 from dat
dat_2a_30_70 <- dat_2a[dat_2a$lm_ht_resid > lm_ht_perc[1] & dat_2a$lm_ht_resid < lm_ht_perc[2],]

# find 25th and 75th percentile
lm_ht_perc <- quantile(dat_2a$lm_ht_resid, probs = c(.25, .75))

# remove lower 10 and upper 10 from dat
dat_2a_25_75 <- dat_2a[dat_2a$lm_ht_resid > lm_ht_perc[1] & dat_2a$lm_ht_resid < lm_ht_perc[2],]

# B. Linear regression of Canopy Closure vs. LiDAR Canopy Cover

# check values by POLYTYPE
table(dat$CC[dat$POLYTYPE == 'FOR'])
table(dat$CC[dat$POLYTYPE == 'WAT'])
table(dat$CC[dat$POLYTYPE == 'DAL'])
table(dat$CC[dat$POLYTYPE == 'GRS'])
table(dat$CC[dat$POLYTYPE == 'ISL'])
table(dat$CC[dat$POLYTYPE == 'UCL'])
table(dat$CC[dat$POLYTYPE == 'BSH'])
table(dat$CC[dat$POLYTYPE == 'OMS'])
table(dat$CC[dat$POLYTYPE == 'RCK'])
table(dat$CC[dat$POLYTYPE == 'TMS'])

# Not zero in all...need to think about whether to leave some in or not

# remove rows with missing values
dat_2b <- dat[is.na(dat$CC) == F & is.na(dat$cc_med) == F,]

# only keep land polygons
dat_2b <- dat_2b[dat_2b$POLYTYPE != 'WAT',]

# plot variables against each other
plot(dat_2b$cc_med, dat_2b$CC)
plot(dat_2b$cc_med[dat_2b$POLYTYPE == 'FOR'], dat_2b$CC[dat_2b$POLYTYPE == 'FOR'])

# run simple linear model
lm_cc <- lm(CC ~ cc_med, data = dat_2b)
summary(lm_cc)

# bind residuals to data
dat_2b <- dat_2b %>% add_column(lm_cc_resid = resid(lm_cc))

# find 10th and 90th percentile
lm_cc_perc <- quantile(dat_2b$lm_cc_resid, probs = c(.10, .90))

# remove lower 20 and upper 20 from dat
dat_2b_10_90 <- dat_2b[dat_2b$lm_cc_resid > lm_cc_perc[1] & dat_2b$lm_cc_resid < lm_cc_perc[2],]

# find 20th and 80th percentile
lm_cc_perc <- quantile(dat_2b$lm_cc_resid, probs = c(.20, .80))

# remove lower 20 and upper 20 from dat
dat_2b_20_80 <- dat_2b[dat_2b$lm_cc_resid > lm_cc_perc[1] & dat_2b$lm_cc_resid < lm_cc_perc[2],]

# find 30th and 70th percentile
lm_cc_perc <- quantile(dat_2b$lm_cc_resid, probs = c(.30, .70))

# remove lower 30 and upper 30 from dat
dat_2b_30_70 <- dat_2b[dat_2b$lm_cc_resid > lm_cc_perc[1] & dat_2b$lm_cc_resid < lm_cc_perc[2],]

# find 25th and 75th percentile
lm_cc_perc <- quantile(dat_2b$lm_cc_resid, probs = c(.25, .75))

# remove lower 20 and upper 20 from dat
dat_2b_25_75 <- dat_2b[dat_2b$lm_cc_resid > lm_cc_perc[1] & dat_2b$lm_cc_resid < lm_cc_perc[2],]

# re plot relationship
plot(dat_2b$cc_med, dat_2b$CC)

# check values by POLYTYPE
table(dat_2b$CC[dat_2b$POLYTYPE == 'FOR'])
table(dat_2b$CC[dat_2b$POLYTYPE == 'WAT'])
table(dat_2b$CC[dat_2b$POLYTYPE == 'DAL'])
table(dat_2b$CC[dat_2b$POLYTYPE == 'GRS'])
table(dat_2b$CC[dat_2b$POLYTYPE == 'ISL'])
table(dat_2b$CC[dat_2b$POLYTYPE == 'UCL'])
table(dat_2b$CC[dat_2b$POLYTYPE == 'BSH'])
table(dat_2b$CC[dat_2b$POLYTYPE == 'OMS'])
table(dat_2b$CC[dat_2b$POLYTYPE == 'RCK'])
table(dat_2b$CC[dat_2b$POLYTYPE == 'TMS'])


# C. Linear regression of Basal Area vs. LiDAR Basal Area

# check values by POLYTYPE
table(dat$BA[dat$POLYTYPE == 'FOR'])
table(dat$BA[dat$POLYTYPE == 'WAT'])
table(dat$BA[dat$POLYTYPE == 'DAL'])
table(dat$BA[dat$POLYTYPE == 'GRS'])
table(dat$BA[dat$POLYTYPE == 'ISL'])
table(dat$BA[dat$POLYTYPE == 'UCL'])
table(dat$BA[dat$POLYTYPE == 'BSH'])
table(dat$BA[dat$POLYTYPE == 'OMS'])
table(dat$BA[dat$POLYTYPE == 'RCK'])
table(dat$BA[dat$POLYTYPE == 'TMS'])

# check same values for LiDAR calculated attribute
table(dat$ba_med[dat$POLYTYPE == 'FOR'])
table(dat$ba_med[dat$POLYTYPE == 'WAT'])
table(dat$ba_med[dat$POLYTYPE == 'DAL'])
table(dat$ba_med[dat$POLYTYPE == 'GRS'])
table(dat$ba_med[dat$POLYTYPE == 'ISL'])
table(dat$ba_med[dat$POLYTYPE == 'UCL'])
table(dat$ba_med[dat$POLYTYPE == 'BSH'])
table(dat$ba_med[dat$POLYTYPE == 'OMS'])
table(dat$ba_med[dat$POLYTYPE == 'RCK'])
table(dat$ba_med[dat$POLYTYPE == 'TMS'])

# some interpreter polygons not zero but still issue with matching...

# remove rows with missing values
dat_2c <- dat[is.na(dat$BA) == F & is.na(dat$ba_med) == F,]

# only keep land polygons
dat_2c <- dat_2c[dat_2c$POLYTYPE != 'WAT',]

# plot variables against each other
plot(dat_2c$ba_med, dat_2c$BA)
plot(dat_2c$ba_med[dat_2c$POLYTYPE == 'FOR'], dat_2c$BA[dat_2c$POLYTYPE == 'FOR'])

# run simple linear model
lm_ba <- lm(BA ~ ba_med, data = dat_2c)
summary(lm_ba)

# bind residuals to data
dat_2c <- dat_2c %>% add_column(lm_ba_resid = resid(lm_ba))

# find 10th and 90th percentile
lm_ba_perc <- quantile(dat_2c$lm_ba_resid, probs = c(.10, .90))

# remove lower 10 and upper 10 from dat
dat_2c_10_90 <- dat_2c[dat_2c$lm_ba_resid > lm_ba_perc[1] & dat_2c$lm_ba_resid < lm_ba_perc[2],]

# find 20th and 80th percentile
lm_ba_perc <- quantile(dat_2c$lm_ba_resid, probs = c(.20, .80))

# remove lower 20 and upper 20 from dat
dat_2c_20_80 <- dat_2c[dat_2c$lm_ba_resid > lm_ba_perc[1] & dat_2c$lm_ba_resid < lm_ba_perc[2],]

# find 30th and 70th percentile
lm_ba_perc <- quantile(dat_2c$lm_ba_resid, probs = c(.30, .70))

# remove lower 30 and upper 30 from dat
dat_2c_30_70 <- dat_2c[dat_2c$lm_ba_resid > lm_ba_perc[1] & dat_2c$lm_ba_resid < lm_ba_perc[2],]

# find 25th and 75th percentile
lm_ba_perc <- quantile(dat_2c$lm_ba_resid, probs = c(.25, .75))

# remove lower 25 and upper 25 from dat
dat_2c_25_75 <- dat_2c[dat_2c$lm_ba_resid > lm_ba_perc[1] & dat_2c$lm_ba_resid < lm_ba_perc[2],]

# re plot relationship
plot(dat_2c$ba_med, dat_2c$BA)

table(dat_2c$BA[dat_2c$POLYTYPE == 'FOR'])
table(dat_2c$BA[dat_2c$POLYTYPE == 'WAT'])
table(dat_2c$BA[dat_2c$POLYTYPE == 'DAL'])
table(dat_2c$BA[dat_2c$POLYTYPE == 'GRS'])
table(dat_2c$BA[dat_2c$POLYTYPE == 'ISL'])
table(dat_2c$BA[dat_2c$POLYTYPE == 'UCL'])
table(dat_2c$BA[dat_2c$POLYTYPE == 'BSH'])
table(dat_2c$BA[dat_2c$POLYTYPE == 'OMS'])
table(dat_2c$BA[dat_2c$POLYTYPE == 'RCK'])
table(dat_2c$BA[dat_2c$POLYTYPE == 'TMS'])

####################################
### COMBINE INTERSECTION OF DATA ###
####################################

# # combine results from part 1 (POLYID column only)
# dat_1 <- intersect(dat_1b %>% subset(select = POLYID), dat_1c %>% subset(select = POLYID))
# 
# # combine results from part 2 (POLYID column only)
# dat_2 <- intersect(dat_2a %>% subset(select = POLYID), dat_2b %>% subset(select = POLYID)) %>% 
#   intersect(., dat_2c %>% subset(select = POLYID))
# 
# # combine all
# dat_screen <- intersect(dat_1, dat_2)
# 
# # reload attributes
# dat_screen <- dat[dat$POLYID %in% dat_screen$POLYID,]

# combine results from part 1 (POLYID column only)
dat_1 <- dat_1b %>% subset(select = POLYID)

# combine results from part 2 (POLYID column only)
# four different thresholds on residuals
dat_2_10_90 <- intersect(dat_2a_10_90 %>% subset(select = POLYID), dat_2b_10_90 %>% subset(select = POLYID)) %>%
  intersect(., dat_2c_10_90 %>% subset(select = POLYID))

dat_2_20_80 <- intersect(dat_2a_20_80 %>% subset(select = POLYID), dat_2b_20_80 %>% subset(select = POLYID)) %>%
  intersect(., dat_2c_20_80 %>% subset(select = POLYID))

dat_2_30_70 <- intersect(dat_2a_30_70 %>% subset(select = POLYID), dat_2b_30_70 %>% subset(select = POLYID)) %>%
  intersect(., dat_2c_30_70 %>% subset(select = POLYID))

dat_2_25_75 <- intersect(dat_2a_25_75 %>% subset(select = POLYID), dat_2b_25_75 %>% subset(select = POLYID)) %>%
  intersect(., dat_2c_25_75 %>% subset(select = POLYID))

# combine all
dat_screen_10_90 <- intersect(dat_1, dat_2_10_90)
dat_screen_20_80 <- intersect(dat_1, dat_2_20_80)
dat_screen_30_70 <- intersect(dat_1, dat_2_30_70)
dat_screen_25_75 <- intersect(dat_1, dat_2_25_75)

# reload attributes
dat_screen_10_90 <- dat[dat$POLYID %in% dat_screen_10_90$POLYID,]
dat_screen_20_80 <- dat[dat$POLYID %in% dat_screen_20_80$POLYID,]
dat_screen_30_70 <- dat[dat$POLYID %in% dat_screen_30_70$POLYID,]
dat_screen_25_75 <- dat[dat$POLYID %in% dat_screen_25_75$POLYID,]

# check some values/distributions
table(dat_screen_10_90$POLYTYPE)
table(dat_screen_20_80$POLYTYPE)
table(dat_screen_25_75$POLYTYPE)
table(dat_screen_30_70$POLYTYPE)
table(dat$POLYTYPE)
# we lost a lot of the variability in polygon type probably because of the linear regression models not
# matching in those other poly types

# write to disk
write.csv(dat_screen_10_90, file = 'D:/ontario_inventory/imputation/dat_screen_1b_2_10_90.csv', row.names = F)
write.csv(dat_screen_20_80, file = 'D:/ontario_inventory/imputation/dat_screen_1b_2_20_80.csv', row.names = F)
write.csv(dat_screen_30_70, file = 'D:/ontario_inventory/imputation/dat_screen_1b_2_30_70.csv', row.names = F)
write.csv(dat_screen_25_75, file = 'D:/ontario_inventory/imputation/dat_screen_1b_2_25_75.csv', row.names = F)

###############################################
### EXPORT USING PART 2 DATA SCREENING ONLY ###
###############################################

# # combine results from part 1 (POLYID column only)
# dat_1 <- intersect(dat_1b %>% subset(select = POLYID), dat_1c %>% subset(select = POLYID))
# 
# # combine results from part 2 (POLYID column only)
# dat_2 <- intersect(dat_2a %>% subset(select = POLYID), dat_2b %>% subset(select = POLYID)) %>% 
#   intersect(., dat_2c %>% subset(select = POLYID))
# 
# # combine all
# dat_screen <- intersect(dat_1, dat_2)
# 
# # reload attributes
# dat_screen <- dat[dat$POLYID %in% dat_screen$POLYID,]

# combine results from part 1 (POLYID column only)
dat_1 <- dat_1b %>% subset(select = POLYID)

# combine results from part 2 (POLYID column only)
# four different thresholds on residuals
dat_2_10_90 <- intersect(dat_2a_10_90 %>% subset(select = POLYID), dat_2b_10_90 %>% subset(select = POLYID)) %>%
  intersect(., dat_2c_10_90 %>% subset(select = POLYID))

dat_2_20_80 <- intersect(dat_2a_20_80 %>% subset(select = POLYID), dat_2b_20_80 %>% subset(select = POLYID)) %>%
  intersect(., dat_2c_20_80 %>% subset(select = POLYID))

dat_2_30_70 <- intersect(dat_2a_30_70 %>% subset(select = POLYID), dat_2b_30_70 %>% subset(select = POLYID)) %>%
  intersect(., dat_2c_30_70 %>% subset(select = POLYID))

dat_2_25_75 <- intersect(dat_2a_25_75 %>% subset(select = POLYID), dat_2b_25_75 %>% subset(select = POLYID)) %>%
  intersect(., dat_2c_25_75 %>% subset(select = POLYID))

# combine all
dat_screen_10_90 <- dat_2_10_90
dat_screen_20_80 <- dat_2_20_80
dat_screen_30_70 <- dat_2_30_70
dat_screen_25_75 <- dat_2_25_75

# reload attributes
dat_screen_10_90 <- dat[dat$POLYID %in% dat_screen_10_90$POLYID,]
dat_screen_20_80 <- dat[dat$POLYID %in% dat_screen_20_80$POLYID,]
dat_screen_30_70 <- dat[dat$POLYID %in% dat_screen_30_70$POLYID,]
dat_screen_25_75 <- dat[dat$POLYID %in% dat_screen_25_75$POLYID,]

# check some values/distributions
table(dat_screen_10_90$POLYTYPE)
table(dat_screen_20_80$POLYTYPE)
table(dat_screen_25_75$POLYTYPE)
table(dat_screen_30_70$POLYTYPE)
table(dat$POLYTYPE)
# we lost a lot of the variability in polygon type probably because of the linear regression models not
# matching in those other poly types

# write to disk
write.csv(dat_screen_10_90, file = 'D:/ontario_inventory/imputation/dat_screen_2_10_90.csv', row.names = F)
write.csv(dat_screen_20_80, file = 'D:/ontario_inventory/imputation/dat_screen_2_20_80.csv', row.names = F)
write.csv(dat_screen_30_70, file = 'D:/ontario_inventory/imputation/dat_screen_2_30_70.csv', row.names = F)
write.csv(dat_screen_25_75, file = 'D:/ontario_inventory/imputation/dat_screen_2_25_75.csv', row.names = F)


# save working image to speed up future changes
# save.image('D:/ontario_inventory/imputation/photo_data_screening.RData')
