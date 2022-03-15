# this code loads polygons segmented from lidar to assess distribution
# of dominant landcover and forested polygons in different
# FRI polytypes

# load packages
library(terra)
library(tidyverse)
library(sf)
library(janitor)

##################################################
### ANALYZE FORESTED POLYGONS NON-FOR POLYTYPE ###
##################################################

# load polygons
poly <- vect('D:/ontario_inventory/segmentation/mask_wat_ucl/ms_5_20_50_add_pt.shp')

# convert to df
p <- as.data.frame(poly)

# pre allocate lists
dom_lc <- list()
dom_for <- list()

# loop through unique POLYTYPES to create table of dom_lc and dom_for
for(i in unique(p$POLYTYPE)){
  
  # subset by polytype
  psub <- p %>% filter(POLYTYPE == i)
  
  # add dom_lc and dom_for tables to list
  dom_lc[[i]] <- psub %>% tabyl(dom_lc) %>% arrange(desc(n)) %>% mutate(percent = round(percent*100,2))
  dom_for[[i]] <- psub %>% tabyl(dom_for) %>% arrange(desc(n)) %>% mutate(percent = round(percent*100,2))
  
}
