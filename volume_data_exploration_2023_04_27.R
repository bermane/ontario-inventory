# this code explores real volume data from cutblocks in the RMF over
# the last few years. GreenFirst has compared volume to FRI estimates
# and we can compare to EFI estimates?

# load packages
library(tidyverse)
library(terra)
library(magrittr)
library(exactextractr)
library(sf)
library(MISTR)

# load FRI
fri <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

# load GRM polygons and project to FRI CRS
grm <- vect('D:/ontario_inventory/SEGMENTATION OUTPUTS/RMF/grm_imp.shp') %>% 
  project(fri)

# load cut block shape file and project to FRI CRS
cuts <- vect('D:/ontario_inventory/Volume Data from Grant/RMF_Harvest/RMF_FMP2019_Harvest.shp') %>% 
  project(fri)

# load volume from EFI
vol <- rast('D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_V_ha.tif') %>%
  project(fri)

# convert to m3
vol <- vol * 0.04

# load merchantable volume from EFI
vol_merch <- rast('D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_Vmerch_ha.tif') %>%
  project(fri)

# convert to m3
vol_merch <- vol_merch * 0.04

####################################
### EXTRACT EFI VOLUME INTO CUTS ###
####################################

# we want the sum of volume based on % coverage

# convert cuts to sf
cuts_sf <- st_as_sf(cuts)

#extract median values of volume and merchantable volume
v <- exact_extract(vol, cuts_sf, 'sum')
v_merch <- exact_extract(vol_merch, cuts_sf, 'sum')

# add new column into cuts
cuts$v_efi <- v
cuts$vmerch_efi <- v_merch

######################################
### GET EFI VOLUME SUM BY BLOCK ID ###
######################################

# get cuts df
cuts_df <- as.data.frame(cuts)

# get efi_vol by BLOCKID
vol_efi_block <- cuts_df %>% 
  group_by(BLOCKID) %>% 
  mutate(v_efi_block = sum(v_efi),
         vmerch_efi_block = sum(vmerch_efi),
         area = sum(Shape_Area),
         n_ha = sum(Shape_Area)/10000) %>% 
  select(BLOCKID, v_efi_block, vmerch_efi_block, area, n_ha) %>% 
  distinct

####################################
### EXTRACT FRI VOLUME INTO CUTS ###
####################################

# convert fri to sf
fri_sf <- st_as_sf(fri)

# check CRS
st_crs(fri_sf)
st_crs(cuts_sf)

# run the intersect function, converting the output to a tibble in the process
int <- as_tibble(st_intersection(fri_sf, cuts_sf))

# calculate gross total volume
int$gtv <- NA
for(i in 1:NROW(int)){
  # skip NA and 'BOG'
  if(is.na(int$SFU[i]) == F){
    if(int$SFU[i] != 'BOG'){
      int$gtv[i] <- calc_GTV(int$BA[i], int$HT[i], int$SFU[i])
    }
  }
}

# add in an area count column to the tibble (area of each arable poly in intersect layer)
int$area_cuts <- st_area(int$geometry)

# subset based on cuts OBJECTID
int %<>% filter(OBJECTID %in% cuts_df$OBJECTID)

# calculate area sum based on OBJECTID
int_area <- int %>% group_by(OBJECTID) %>%
  summarise(area = sum(area_cuts))

#group data by county area and calculate the total arable land area per county
#output as new tibble
tb_ArableByCounty <- int %>%
  group_by(County_UA) %>%
  summarise(areaArable = sum(areaArable))

#change data type of areaArable field to numeric (from 'S3: units' with m^2 suffix)
tb_ArableByCounty$areaArable <- as.numeric(tb_ArableByCounty$areaArable)

#join tibble to original county polygon shapefile and export as new shapefile
shp_out <- st_write(merge(counties, tb_ArableByCounty, by = 'County_UA'), "ArablebyCounty.shp")

########################################
### EXTRACT VOLUME INTO GRM POLYGONS ###
########################################