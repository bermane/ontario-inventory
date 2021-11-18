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

# find 20th and 80th percentile
lm_ht_perc <- quantile(dat_2a$lm_ht_resid, probs = c(.10, .90))

# remove lower 20 and upper 20 from dat
dat_2a <- dat_2a[dat_2a$lm_ht_resid > lm_ht_perc[1] & dat_2a$lm_ht_resid < lm_ht_perc[2],]

# re plot relationship
plot(dat_2a$lor_med, dat_2a$HT)

# run simple linear model on remaining data
lm_ht2 <- lm(HT ~ lor_med, data = dat_2a)
summary(lm_ht2)

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

# find 20th and 80th percentile
lm_cc_perc <- quantile(dat_2b$lm_cc_resid, probs = c(.10, .90))

# remove lower 20 and upper 20 from dat
dat_2b <- dat_2b[dat_2b$lm_cc_resid > lm_cc_perc[1] & dat_2b$lm_cc_resid < lm_cc_perc[2],]

# re plot relationship
plot(dat_2b$cc_med, dat_2b$CC)

# run simple linear model on remaining data
lm_cc2 <- lm(CC ~ cc_med, data = dat_2b)
summary(lm_cc2)

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

# find 20th and 80th percentile
lm_ba_perc <- quantile(dat_2c$lm_ba_resid, probs = c(.10, .90))

# remove lower 20 and upper 20 from dat
dat_2c <- dat_2c[dat_2c$lm_ba_resid > lm_ba_perc[1] & dat_2c$lm_ba_resid < lm_ba_perc[2],]

# re plot relationship
plot(dat_2c$ba_med, dat_2c$BA)

# run simple linear model on remaining data
lm_ba2 <- lm(BA ~ ba_med, data = dat_2c)
summary(lm_ba2)

###COMBINE IN

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

# combine results from pa rt 2 (POLYID column only)
dat_2 <- intersect(dat_2a %>% subset(select = POLYID), dat_2b %>% subset(select = POLYID)) %>%
  intersect(., dat_2c %>% subset(select = POLYID))

# combine all
dat_screen <- intersect(dat_1, dat_2)

# reload attributes
dat_screen <- dat[dat$POLYID %in% dat_screen$POLYID,]

# check some values/distributions
table(dat_screen$POLYTYPE)
table(dat$POLYTYPE)
# we lost a lot of the variability in polygon type probably because of the linear regression models not
# matching in those other poly types

# write to disk
write.csv(dat_screen, file = 'D:/ontario_inventory/imputation/dat_screen_1b_2_10perc.csv', row.names = F)

# save working image to speed up future changes
save.image('D:/ontario_inventory/imputation/photo_data_screening.RData')
