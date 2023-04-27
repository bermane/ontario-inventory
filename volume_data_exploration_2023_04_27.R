# this code explores real volume data from cutblocks in the RMF over
# the last few years. GreenFirst has compared volume to FRI estimates
# and we can compare to EFI estimates?

# load packages
library(tidyverse)
library(terra)
library(magrittr)
library(exactextractr)
library(sf)

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

####################################
### EXTRACT EFI VOLUME INTO CUTS ###
####################################

# we want the sum of volume based on % coverage

# convert cuts to sf
cuts_sf <- st_as_sf(cuts)

#extract median values
vol <- exact_extract(vol, cuts_sf, 'sum')

# add new column into cuts
cuts$vol_efi <- vol

######################################
### GET EFI VOLUME SUM BY BLOCK ID ###
######################################

# get cuts df
cuts_df <- as.data.frame(cuts)

# get efi_vol by BLOCKID
vol_efi_block <- cuts_df %>% 
  group_by(BLOCKID) %>% 
  mutate(vol_efi_block = sum(vol_efi)) %>% 
  select(BLOCKID, vol_efi_block) %>% 
  distinct

########################################
### EXTRACT VOLUME INTO GRM POLYGONS ###
########################################