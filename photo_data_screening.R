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

# 1. Polygon area < 50 ha (500000 m^2)

# check number of polygons with area < 50 ha
NROW(dat[dat$AREA < 500000,])

# subset data
dat <- dat[dat$AREA < 500000,]

# 2. Number of micro stands in polygon < 15
# skip since we don't have micro stands

# 3. Within-polygon height and canopy cover coefficient of variation < 0.5