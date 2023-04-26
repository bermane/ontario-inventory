# This code loads photo interpreted polygons from the FSF in Ontario
# and applies a number of data screens to curate an accurate 
# and reliable dataset

# load packages
library(terra)
library(tidyverse)
library(exactextractr)
library(sf)
library(magrittr)

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/FSF/FRI/FSF_opi_polygon_CSRS_NAD83_17.shp')

# convert to df
dat <- as.data.frame(poly)

#############################
### DATA SCREENING PART 1 ###
#############################

# remove all non-forested polygons
dat <- filter(dat, POLYTYPE == 'FOR')

# create smaller polygon set only FOR polytypes
poly_dat <- poly[poly$POLYTYPE == 'FOR']

# export forested polygons only for maps
# writeVector(poly_dat, filename = 'D:/ontario_inventory/imputation/distributions_for_only/vector/fri_polygons_for_only.shp')

#############################
### DATA SCREENING PART 2 ###
#############################

# polygon landcover > 50% forested

# load VLCE 2.0 landcover dataset from 2018
lc <- rast('D:/ontario_inventory/VLCE/CA_forest_VLCE2_2018.tif')

# project poly to crs of raster
poly_lc <- project(poly_dat, lc)

# convert to sf
poly_lcsf <- st_as_sf(poly_lc)

# extract landcover values
lc_poly <- exact_extract(lc, poly_lcsf)

# set landcover class key with single forested class
lc_key_for <- c(`0` = 'NA',
                `20` = 'Water',
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

# add new columns into dat
dat <- dat %>% add_column(dom_for = lc_dom_for)

# check number of rows with forested lc > 50%
NROW(dat[dat$dom_for == 'Yes',])

# subset dat
dat_domfor <- dat[dat$dom_for == 'Yes',]

#############################
### DATA SCREENING PART 3 ###
#############################

# p95 > 5 m
# cc > 50%

# Load cc and p95
lidar <- c('p95' = 'D:/ontario_inventory/FSF/ALS/zq95.img', 
           'cc' = 'D:/ontario_inventory/FSF/ALS/cov_2m.img')

# load LiDAR rasters as raster stack
lidar_ras <- rast(lidar)

# change layer names
names(lidar_ras) <- names(lidar)

# project poly to crs of raster
poly_ras <- project(poly_dat, lidar_ras)

# convert to sf
poly_ras <- st_as_sf(poly_ras)

# extract median values
vec <- exact_extract(lidar_ras, poly_ras, 'median')

# change column names
colnames(vec) <- names(lidar)

# add new column into dat
dat <- cbind(dat, vec)

# How many rows have p95 > 5, which is definition of forest
NROW(dat[dat$p95 >= 5,])

# require p95 > 5
dat_p95 <- dat[dat$p95 >= 5,]

# How many rows have cc > 50
NROW(dat[dat$cc >= 50,])

# require cc > 50
dat_cc <- dat[dat$cc >= 50,]

####################################
### COMBINE INTERSECTION OF DATA ###
####################################

# combine all
dat_screen <- intersect(dat_domfor %>% subset(select = POLYID), dat_p95 %>% subset(select = POLYID)) %>%
  intersect(., dat_cc %>% subset(select = POLYID))

# reload attributes
dat_screen <- dat[dat$POLYID %in% dat_screen$POLYID,]

# write to disk
write.csv(dat_screen, file = 'D:/ontario_inventory/imputation/fsf/dat_screen_fri_fsf.csv', row.names = F)
