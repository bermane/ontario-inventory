# This code adds 3 and 5 species classification to the Romeo
# FRI Polygons to be used for visuals

# load packages
library(terra)
library(tidyverse)
library(magrittr)
library(janitor)

############################
### LOAD POLYGON DATASET ###
############################

# load photo interpreted polygons
poly <-
  vect(
    'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp'
  )

# make as dataframe
dat <- as.data.frame(poly)

###################
### ADD SP COMP ###
###################

# load sp_group data
sp_group <- read.csv(
  'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest_SPGROUP.shp'
)

# change POLYID to numeric
dat %<>% mutate(POLYID = as.numeric(POLYID))
sp_group %<>% mutate(POLYID = as.numeric(POLYID))

# join to dat
dat <- left_join(dat, 
                 sp_group %>% select(POLYID, SpeciesGroup2, SpeciesGroup3), 
                 by = 'POLYID')

# rename species groups
dat %<>% rename(species_5_class = SpeciesGroup2,
                species_3_class = SpeciesGroup3)

# add back to polygons
values(poly) <- dat

# write shapefile
writeVector(poly, filename = 'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest_SPClass.shp')
