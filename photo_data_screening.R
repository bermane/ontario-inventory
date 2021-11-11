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

# EXTRA: Some sort of criteria based on year the polygon was updated? 
# i.e. how closely it matches the timing of the LiDAR acquisition

# do we only want polygons that are listed as FOREST? probably...