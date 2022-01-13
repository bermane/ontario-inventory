# this code compares attributes from FRI and LiDAR generated and imputed forest stands
# across a grid of randomly sampled points

# load packages
library(terra)
library(tidyverse)
library(sp)
library(sf)

# load two polygon datasets
fri <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')
lid <- vect('D:/ontario_inventory/imputation/ms_10_10_100_agg_na_IMP.shp')

# remove most of the data for easy conversion to sp
fri_sp <- fri
values(fri_sp) <- values(fri)[,1]

# convert polygons to sp objects
fri_sp <- as(fri_sp, 'Spatial')

# convert to same CRS
fri <- project(fri, lid)

# create vector of extent
extent <- as.polygons(ext(lid))
crs(extent) <- crs(lid)

# convert extent to sp
extent_sp <- as(extent, 'Spatial')

# Create a grid of points
grid <- spsample(extent_sp, n = 500, type = 'regular')

# convert grid to spatvector and change crs
grid_sv <- vect(grid)

# get polygon values at points, fri first
int_fri <- terra::intersect(project(grid_sv, fri), fri)

# put values into df
df_fri <- values(int_fri)

# create new set of points where there was data in fri
grid_fri <- int_fri
values(grid_fri) <- NULL

#get polygon values at points of lid
int_lid <- terra::intersect(project(grid_fri, lid), lid)

# put values into df
df_lid <- values(int_lid)

# create df for comparison
poly_df <- data.frame(ht_fri = df_fri$HT,
                      ht_lid = df_lid$HT,
                      cc_fri = df_fri$CC,
                      cc_lid = df_lid$CC,
                      ba_fri = df_fri$BA,
                      ba_lid = df_lid$BA,
                      age_fri = df_fri$AGE,
                      age_lid = df_lid$AGE,
                      spcomp_fri = df_fri$SPCOMP,
                      spcomp_lid = df_lid$SPCOMP)

# plot relationships
plot(poly_df$ht_fri, poly_df$ht_lid)
plot(poly_df$cc_fri, poly_df$cc_lid)
plot(poly_df$ba_fri, poly_df$ba_lid)
plot(poly_df$age_fri, poly_df$age_lid)
