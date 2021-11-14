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
NROW(dat[dat$AREA < 500000,])

# subset data
dat_1a <- dat[dat$AREA < 500000,]

# look AT HISTOGRAM OF VERY LARGE POLYGONS AND GET RID OF THEM and ONLY FOR POLYGONS!!!

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

# load VLCE landcover dataset (2010 at the moment)
lc <- rast('D:/ontario_inventory/VLCE/LC_Class_HMM_17S_2010.dat')

# project poly to crs of raster
poly_lc <- project(poly, lc)

# extract lc types within each polygon
lc_poly <- terra::extract(lc, poly_lc)

# define treed classes
treed <- c('Mixed Wood', 'Broadleaf', 'Coniferous', 'Wetland-Treed')

# calculate number of unique lc types in each polygon and % of polygon covered by treed classes
lc_poly <- lc_poly %>% group_by(ID) %>%
  mutate(num_uniq = length(unique(category))) %>%
  mutate(perc_for = sum(category %in% treed)/length(category))

# drop duplicate rows
lc_poly <- lc_poly[duplicated(lc_poly$ID)==F,]

# add new columns into dat
dat <- dat %>% add_column(num_uniq = lc_poly$num_uniq,
                          perc_for = lc_poly$perc_for)

# rm variables
rm(lc, lc_poly, poly_lc, treed)

# check number of polygons with LC types < 2
# perhaps Chen meant LC types <= 2??
NROW(lc_poly[lc_poly$num_uniq < 2,])
NROW(lc_poly[lc_poly$num_uniq <= 2,])

# check number of polygons with forested classes > 50%
NROW(lc_poly[lc_poly$perc_for > .5,])

# subset dat based on lc info
dat_1c <- dat[dat$num_uniq <=2 & dat$perc_for > .5,]

###########################
###DATA SCREENING PART 2###
###########################

# A. Linear regression of Height vs. Lorey's Height

# WE PROBABLY WANT TO RUN THIS WITHOUT WATER POLYGONS!!
# PRODUCTIVE FOREST ONLY? OR ALL LAND POLYGONS? OR PROD and NONPROD FOREST ONLY?

# Load Lorey's Height from LiDAR
lor <- rast('D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_lor.tif')

# project poly to crs of rasters
poly_lor <- project(poly, lor)

# extract median within each polygon
# chen did not use edge pixels -- here edge pixels are still included so can consider
# writing code to remove edge pixels from the calculation
lor_med <- terra::extract(lor, poly_lor, fun = function(x){median(x, na.rm = T)})

# add new column into dat
dat <- dat %>% add_column(lor_med = lor_med[,2])

# remove layers
rm(lor, lor_med, poly_lor)

# remove rows with missing values
dat_2a <- dat[is.na(dat$HT) == F & is.na(dat$lor_med) == F,]

# plot variables against each other
plot(dat$lor_med, dat$HT)
plot(dat$lor_med[dat$POLYTYPE == 'FOR'], dat$HT[dat$POLYTYPE == 'FOR'])

# run simple linear model
lm_ht <- lm(HT ~ lor_med, data = dat_2a)
summary(lm_ht)

# bind residuals to data
dat_2a <- dat_2a %>% add_column(lm_ht_resid = resid(lm_ht))

# find 20th and 80th percentile
lm_ht_perc <- quantile(dat_2a$lm_ht_resid, probs = c(.20, .80))

# remove lower 20 and upper 20 from dat
dat_2a <- dat_2a[dat_2a$lm_ht_resid > lm_ht_perc[1] & dat_2a$lm_ht_resid < lm_ht_perc[2],]

# re plot relationship
plot(dat_2a$lor_med, dat_2a$HT)

# B. Linear regression of Canopy Closure vs. LiDAR Canopy Cover

# WE PROBABLY WANT TO RUN THIS WITHOUT WATER POLYGONS!!
# PRODUCTIVE FOREST ONLY? OR ALL LAND POLYGONS? OR PROD and NONPROD FOREST ONLY?

# Load canopy cover from LiDAR
cc <- rast('D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif')

# project poly to crs of rasters
poly_cc <- project(poly, cc)

# extract median within each polygon
# chen did not use edge pixels -- here edge pixels are still included so can consider
# writing code to remove edge pixels from the calculation
cc_med <- terra::extract(cc, poly_cc, fun = function(x){median(x, na.rm = T)})

# add new column into dat
dat <- dat %>% add_column(cc_med = cc_med[,2])

# remove layers
rm(cc, cc_med, poly_cc)

# remove rows with missing values
dat_2b <- dat[is.na(dat$CC) == F & is.na(dat$cc_med) == F,]

# plot variables against each other
plot(dat$cc_med, dat$CC)
plot(dat$cc_med[dat$POLYTYPE == 'FOR'], dat$CC[dat$POLYTYPE == 'FOR'])

# run simple linear model
lm_cc <- lm(CC ~ cc_med, data = dat_2b)
summary(lm_cc)

# bind residuals to data
dat_2b <- dat_2b %>% add_column(lm_cc_resid = resid(lm_cc))

# find 20th and 80th percentile
lm_cc_perc <- quantile(dat_2b$lm_cc_resid, probs = c(.20, .80))

# remove lower 20 and upper 20 from dat
dat_2b <- dat_2b[dat_2b$lm_cc_resid > lm_cc_perc[1] & dat_2b$lm_cc_resid < lm_cc_perc[2],]

# re plot relationship
plot(dat_2b$cc_med, dat_2b$CC)
