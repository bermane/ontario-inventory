library(terra)
library(lidR)
library(sf)
library(stringr)
library(dplyr)
library(tidyr)
library(doParallel)
library(foreach)
library(ggplot2)
library(readr)

assign_common_name <- function(sp_abbrev) {
  sp_abbrev  <- toupper(sp_abbrev)
  
  dict <- data.frame(SB = "black spruce", 
                     LA = "eastern larch", 
                     BW = "white birch", 
                     BF = "balsam fir", 
                     CE = "cedar", 
                     SW = "white spruce", 
                     PT = "trembling aspen", 
                     PJ = "jack pine", 
                     PO = "poplar", 
                     PB = "balsam poplar", 
                     PR = "red pine", 
                     PW = "white pine", 
                     SX = "spruce", 
                     MR = "red maple", 
                     AB = "black ash", 
                     BY = "yellow birch") %>%
    pivot_longer(everything(), names_to = "abb", values_to = "common")
  
  dict$common[match(sp_abbrev, dict$abb)]
  
}

assign_type <- function(sp_common) {
  sp_common  <- tolower(sp_common)
  ifelse(stringr::str_detect(sp_common, pattern = "pine|spruce|fir|cedar|larch"), "Coniferous", "Deciduous")
}

VRI <- st_read('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

VRI_FOR <- VRI %>%
  st_drop_geometry() %>%
  filter(POLYTYPE == "FOR") %>%
  select(POLYID, POLYTYPE, DEVSTAGE, WG, SPCOMP, YRORG, YRDEP, YRUPD) %>%
  mutate(new_SP = str_match_all(SPCOMP, "[A-Z]{2}[ ]+[0-9]+")) %>%
  unnest(new_SP) %>%
  mutate(new_SP = as.character(new_SP)) %>%
  separate(new_SP, into = c("SP", "PROP")) %>%
  mutate(PROP = as.numeric(PROP), 
         Common = assign_common_name(SP), 
         sp_type = assign_type(Common))

# Polygon-level species groups ---- 


poly_dom_type <- VRI_FOR %>%
  group_by(POLYID) %>%
  summarize(per_conif = sum(PROP[sp_type == "Coniferous"]), 
            per_decid = sum(PROP[sp_type == "Deciduous"]))

poly_dom_sp <- VRI_FOR %>%
  group_by(POLYID) %>%
  slice_max(PROP, n = 1, with_ties = FALSE)


poly_dom_sp_group <- inner_join(poly_dom_type, poly_dom_sp, by = "POLYID")

poly_dom_sp_group <- poly_dom_sp_group %>%
  mutate(SpeciesGroup1 = ifelse(PROP >= 70, Common, 
                                ifelse(PROP < 70 & per_conif >= 70, "Mixed Coniferous",
                                       ifelse(PROP < 70 & per_decid >= 70, "Mixed Deciduous", "Mixedwoods"))), 
         SpeciesGroup2 = ifelse(Common == "jack pine" & PROP >= 50 & per_conif >= 70, "Jack Pine Dominated", ifelse(
           Common == "black spruce" & PROP >= 50 & per_conif >= 70, "Black Spruce Dominated", ifelse(
             per_decid >= 70, "Hardwood", ifelse(
               per_decid >= 30 & per_decid <= 70 & per_conif >= 30 & per_conif <= 70, "Mixedwood", "Mixed Conifers"
             )))), 
         SpeciesGroup3 = ifelse(per_conif >= 70, "Softwood", 
                                ifelse(per_decid >= 70, "Hardwood", "Mixedwood"))) %>%
  mutate(across(.cols = starts_with("SpeciesGroup"), .fns = as.factor))

abbreviation_dict <- data.frame(
  SpeciesGroup2 = c("Black Spruce Dominated", 
                    "Jack Pine Dominated", 
                    "Hardwood", 
                    "Mixedwood", 
                    "Mixed Conifers"), 
  SpeciesGroup2_short = c("BS", "JP", "HW", "MW", "MC")
)

poly_dom_sp_group <- inner_join(poly_dom_sp_group, 
                                abbreviation_dict)

write.csv(poly_dom_sp_group,
          'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest_SPGROUP.shp',
          row.names = F)
